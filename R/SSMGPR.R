#' Kalman Filter Function
#'
#' Implements a Kalman filter for time series data.
#'
#' @param y Numeric vector of observations.
#' @param phi Numeric, the parameter phi in the Kalman filter.
#' @param Sigmav Numeric, the parameter Sigmav in the Kalman filter.
#' @param Sigmaw Numeric, the parameter Sigmaw in the Kalman filter.
#' @param m0 Numeric, initial value of m.
#' @param Sigma0 Numeric, initial value of Sigma.
#' @param n Numeric, number of future time points to predict.
#' @return A list containing updated and predicted means and variances.
#' @export
kalman <- function(y, phi, Sigmav, Sigmaw, m0 = 0, Sigma0 = 1, n = 2) {
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


#' Gaussian Process Regression with Kalman Filter
#'
#' Applies Gaussian Process Regression (GPR) using a Kalman filter approach.
#'
#' @param y Numeric vector, the time series data to model.
#' @param gamma Numeric, the gamma parameter in GPR.
#' @param Sigmaw Numeric, the noise variance parameter.
#' @param m0 Numeric, initial mean.
#' @param Sigma0 Numeric, initial variance.
#' @param n Numeric, number of steps for prediction.
#' @return A list with updated and predicted means and variances.
#' @export
kf.gp <- function(y, gamma, Sigmaw, m0 = 0, Sigma0 = 1, n = 2) {
  T=length(y)
  #update Sigmav and phi
  phi <- exp(-1/gamma)
  Sigmav <- 1-exp(-2/gamma)
  result <- kalman(y,phi=phi,Sigmav=Sigmav,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0,n)

  return (list(mu.p=result$mu.p, Sigma.p=result$Sigma.p,
               mu.f=result$mu.f, Sigma.f=result$Sigma.f))
}

#' Log Likelihood for Gaussian Process Regression
#'
#' Calculates the log likelihood for a given set of parameters in GPR.
#'
#' @param y Numeric vector, the time series data.
#' @param mu.p Numeric vector, the predicted means.
#' @param Sigma.p Numeric vector, the predicted variances.
#' @param Sigmaw Numeric, the noise variance parameter.
#' @param m0 Numeric, initial mean.
#' @param Sigma0 Numeric, initial variance.
#' @return Numeric, the log likelihood value.
#' @export
kf.loglikelihood1 <- function(y, mu.p, Sigma.p, Sigmaw, m0 = 0, Sigma0 = 1) {
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

#' Wrapper for Log Likelihood Calculation in GPR
#'
#' A wrapper function to calculate the log likelihood in GPR.
#'
#' @param y Numeric vector, the time series data.
#' @param gamma Numeric, the gamma parameter in GPR.
#' @param Sigmaw Numeric, the noise variance parameter.
#' @param m0 Numeric, initial mean.
#' @param Sigma0 Numeric, initial variance.
#' @return Numeric, the log likelihood value.
#' @export
kf.loglikelihood <- function(y, gamma, Sigmaw, m0 = 0, Sigma0 = 1) {
  o <- kf.gp(y=y,gamma=gamma,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0)
  mu.p <- o$mu.p
  Sigma.p <- o$Sigma.p
  result <- kf.loglikelihood1(y=y,mu.p=mu.p,Sigma.p=Sigma.p,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0)
  return (result)
}

#' Optimize Parameters for Gaussian Process Regression
#'
#' Finds optimal parameters for GPR using optimization techniques.
#'
#' @param y Numeric vector, the time series data.
#' @return A list containing optimal values for gamma and Sigmaw.
#' @export
optim_parm <- function(y) {
  opt_param <- optim(par = c(5,0.5),
                     fn = function(parm) -1*kf.loglikelihood(y,parm[1],  parm[2]))
  return(list(gamma = opt_param$par[1],
              Sigmaw=opt_param$par[2]))
}
