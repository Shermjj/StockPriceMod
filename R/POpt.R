#' Optimal Portfolio Weights Calculation
#'
#' This function calculates the optimal portfolio weights based on a given covariance matrix, expected returns, and the number of components.
#' It uses a quadratic programming approach to find the solution.
#'
#' @param covar_mat A square, symmetric matrix representing the covariance matrix of asset returns.
#' @param returns A numeric vector of expected returns for each asset.
#' @param num_components An integer indicating the number of assets in the portfolio.
#' @param lambda A numeric value representing the risk aversion parameter (default is 1).
#'
#' @return A numeric vector of portfolio weights that minimizes the risk for a given level of return.
#'
#' @examples
#' covar_mat <- matrix(c(0.1, 0.2, 0.2, 0.3), ncol=2)
#' returns <- c(0.1, 0.15)
#' num_components <- length(returns)
#' port_opt(covar_mat, returns, num_components)
#'
#' @export
#'
port_opt <- function(covar_mat, returns, num_components, lambda = 1) {
  Amat <- matrix(0, num_components, num_components)
  Amat[,1] <- 1
  bvec <- rep(0,num_components)
  bvec[1] <- 1

  QP_result <- solve.QP(2*covar_mat, lambda * returns, matrix(1, nrow=num_components,ncol=1), c(1), meq=1)
  return(QP_result$solution)
}
