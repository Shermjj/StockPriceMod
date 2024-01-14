# Load necessary libraries
library(tidyverse)
library(readr)

# Load the data from a CSV file
stock_prices <- read_csv("sp500_stock_data.csv")

# Convert Date to a Date object if it's not already
stock_prices$Date <- as.Date(stock_prices$Date)

# Filter for specific date range and select only Date, Close, and Ticker
stock_prices_filtered <- stock_prices %>%
  filter(Date >= as.Date("2020-01-02") & Date <= as.Date("2020-01-10")) %>%
  select(Date, Close, Ticker)

# Reshape data to a wide format
wide_data <- stock_prices_filtered %>%
  pivot_wider(names_from = Ticker, values_from = Close)

# Calculate logarithmic returns
log_returns <- wide_data %>%
  mutate(across(-Date, ~log(. / lag(.)))) %>%
  select(-Date) %>%  # Remove the Date column
  na.omit()  # Remove NAs resulting from lag calculation

# Perform PCA on logarithmic returns
pca_result <- prcomp(log_returns, scale = TRUE)

# View summary of PCA results
summary(pca_result)

# Extract and display loadings for the first three PCs
loadings <- as.data.frame(pca_result$rotation[, 1:3])
print(loadings)


# Extract variance explained
var_explained <- pca_result$sdev ^ 2 / sum(pca_result$sdev ^ 2)
cum_var_explained <- cumsum(var_explained)

# Create a dataframe for plotting
pc_var_explained <- data.frame(PC = 1:length(cum_var_explained), 
                               CumulativeVariance = cum_var_explained)

# Plot using ggplot2
ggplot(pc_var_explained, aes(x = PC, y = CumulativeVariance)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:length(cum_var_explained)) +
  labs(title = "Cumulative Variance Explained by Each Principal Component",
       x = "Principal Component",
       y = "Cumulative Variance Explained") +
  theme_minimal()

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

# Project the data onto the first principal component
pc1_scores <- as.data.frame(pca_result$x[,1])

# Add Date back to the projected data for plotting
pc1_scores$Date <- wide_data$Date[2:length(wide_data$Date)]

# Plot the first principal component over time
ggplot(pc1_scores, aes(x = Date, y = pca_result$x[,1])) +  # V1 is the first principal component
  geom_line() +
  labs(title = "Time Series Plot of the First Principal Component",
       x = "Date",
       y = "PC1 Score") +
  theme_minimal()

y <- pc1_scores[,1]
z <- optim_parm(y)
kf.gp(y, z$gamma, z$Sigmaw )
