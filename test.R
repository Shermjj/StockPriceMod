library(quadprog)
library(tidyverse)
library(readr)
set.seed(123)

kalman <- function(y,phi,Sigmav,Sigmaw,m0,Sigma0,n=2){

  T <- length(y)

  #initialization
  mu.p <- rep(NA,T+n)
  Sigma.p <- rep(NA,T+n)
  mu.f <- rep(NA,T)
  Sigma.f <- rep(NA,T)



  #forward recursion time1
  mu.p[1] <- m0
  Sigma.p[1] <- Sigma0
  mu.f[1] <- m0 + (y[1]-m0)*(Sigma0/(Sigma0+Sigmaw))
  Sigma.f[1] <- Sigma0-(Sigma0^2/(Sigma0+Sigmaw))

  #forward recursion time 2:T
  for (t in 2:T){

    #prediction
    mu.p[t] <- phi*mu.f[t-1]
    Sigma.p[t] <- phi^2 * Sigma.f[t-1] + Sigmaw

    #update
    deno <- Sigmaw + Sigma.p[t]
    mu.f[t] <- Sigmaw*mu.p[t]/deno + Sigma.p[t]*y[t]/deno
    Sigma.f[t] <- Sigmaw*Sigma.p[t]/deno
  }
  #predict for T+1:T+n
  for (t in (T+1):(T+n)){
    if (t == T+1){
      mu.p[t] <- phi*mu.f[t-1]
      Sigma.p[t] <- phi^2 * Sigma.f[t-1] + Sigmaw
    }
    else{
      mu.p[t] <- phi*mu.p[t-1]
      Sigma.p[t] <- phi^2 * Sigma.p[t-1] + Sigmaw
    }
  }
  return (list(mu.f=mu.f,Sigma.f=Sigma.f,mu.p=mu.p,Sigma.p=Sigma.p))
}

kf.gp <- function(y,gamma,Sigmaw,m0=0,Sigma0=1,n=2){
  T=length(y)
  #update Sigmav and phi
  phi <- exp(-1/gamma)
  Sigmav <- 1-exp(-2/gamma)
  result <- kalman(y,phi=phi,Sigmav=Sigmav,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0,n)

  return (list(mu.p=result$mu.p, Sigma.p=result$Sigma.p,
               mu.f=result$mu.f, Sigma.f=result$Sigma.f))

}

kf.loglikelihood1 <- function(y,mu.p,Sigma.p,Sigmaw,m0=0,Sigma0=1){
  T <- length(y)
  likelihood <- rep(NA,T)

  #at time 1
  likelihood[1] <- log(dnorm(y[1],mean=m0,sd = sqrt(Sigma0 + Sigmaw)))

  #time 2:T
  for (t in 2:T){
    likelihood[t] <- log(dnorm(y[t],mean=mu.p[t],sd=sqrt(Sigmaw+Sigma.p[t])))
  }
  return (sum(likelihood))
}

kf.loglikelihood <- function(y,gamma,Sigmaw,m0=0,Sigma0=1){
  o <- kf.gp(y=y,gamma=gamma,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0)
  mu.p <- o$mu.p
  Sigma.p <- o$Sigma.p
  result <- kf.loglikelihood1(y=y,mu.p=mu.p,Sigma.p=Sigma.p,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0)
  return (result)
}

optim_parm <- function(y){
  opt_param <- optim(par = c(5,0.5),
                     fn = function(parm) -1*kf.loglikelihood(y,parm[1],  parm[2]))
  return(list(gamma = opt_param$par[1],
              Sigmaw=opt_param$par[2]))
}

# Load the data from a CSV file
stock_prices <- read_csv("./data/sp500_stock_data.csv")

# Convert Date to a Date object if it's not already
stock_prices$Date <- as.Date(stock_prices$Date)

# Using all of 2020
start_date <- as.Date("2020-01-02")
end_date <- as.Date("2020-12-31")

# Filter for specific date range and select only Date, Close, and Ticker
stock_prices_filtered <- stock_prices %>%
  filter(Date >=  start_date & Date <= end_date) %>%
  select(Date, Close, Ticker)

# Reshape data to a wide format
wide_data <- stock_prices_filtered %>%
  pivot_wider(names_from = Ticker, values_from = Close)

wide_data <- wide_data[, !apply(is.na(wide_data), 2, any)]

# Calculate logarithmic returns
log_returns <- wide_data %>%
  mutate(across(-Date, ~log(. / lag(.)))) %>%
  select(-Date) %>%  # Remove the Date column
  na.omit()  # Remove NAs resulting from lag calculation

#Looking at 10 days
start_day <- 200
end_day <- 210

log_returns_monthly <- log_returns[start_day:end_day, ]

# Perform PCA on logarithmic returns
pca_result <- prcomp(log_returns_monthly, scale = TRUE, center= TRUE)


explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_variance <- cumsum(explained_variance)

# Find the number of components that explain at least the threshold variance
num_components <- which(cumulative_variance >= 0.7)[1]
print(num_components)

# Return the principal components
y <- pca_result$x[, 1:num_components]

num_pred_days <- 2

# Initialize an empty list to store the results
mu_p_list <- list()

# Initialize an empty dataframe to store the time series data
time_series_df <- data.frame(matrix(ncol = num_components, nrow = num_pred_days))
colnames(time_series_df) <- paste("Component", 1:num_components, sep = "")

for (m in 1:num_components) {
  # Obtain optimized hyperparameters for the m-th component
  optim_hyperparam <- optim_parm(y[, m])

  # Run the Gaussian Process with Kalman filter
  gp_result <- kf.gp(y[, m], optim_hyperparam$gamma, optim_hyperparam$Sigmaw, n = num_pred_days)
  #  print(gp_result$Sigma.p)
  # Store the last num_pred_days of the time series mu.p from gp_result
  mu_p_list[[m]] <- tail(gp_result$mu.p, num_pred_days)

  # Add the time series to the dataframe
  time_series_df[, m] <- mu_p_list[[m]]
}

agg_exp_returns = apply(time_series_df, MARGIN = 2, FUN = sum)

covar_mat <- diag(pca_result$sdev[1:num_components] ** 2)
Amat <- matrix(0, num_components, num_components)
Amat[,1] <- 1
bvec <- rep(0,num_components)
bvec[1] <- 1
lambda <- 1

QP_result <- solve.QP(2*covar_mat, lambda * agg_exp_returns, matrix(1, nrow=num_components,ncol=1), c(1), meq=1)
QP_result$solution
