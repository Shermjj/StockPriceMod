% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSMGPR.R
\name{kf.loglikelihood}
\alias{kf.loglikelihood}
\title{Wrapper for Log Likelihood Calculation in GPR}
\usage{
kf.loglikelihood(y, gamma, Sigmaw, m0 = 0, Sigma0 = 1)
}
\arguments{
\item{y}{Numeric vector, the time series data.}

\item{gamma}{Numeric, the gamma parameter in GPR.}

\item{Sigmaw}{Numeric, the noise variance parameter.}

\item{m0}{Numeric, initial mean.}

\item{Sigma0}{Numeric, initial variance.}
}
\value{
Numeric, the log likelihood value.
}
\description{
A wrapper function to calculate the log likelihood in GPR.
}
