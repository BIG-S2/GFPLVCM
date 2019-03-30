\name{parameter_estimate_boot}
\alias{parameter_estimate_boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypothesis testing via bootstrap
}
\description{
Test the nullity of both the scalar and the functional covariates.
}
\usage{
parameter_estimate_boot(y, x, z, lambda_s, lambda_u, bd, order_1, order_2, breaks, boot_R, gamma_est_or, b_est_or)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a list of the response with the elements \code{y_ID}, \code{y_time_point} and \code{y_value}.}
  \item{x}{a list of the functional covariate with the elements \code{x_ID}, \code{x_time_point} and \code{x_value}.}
  \item{z}{a list of the scalar covariate with the elements \code{z_ID}, \code{z_time_point} and \code{z_value}.}
  \item{lambda_s}{the tuning parameter in the \code{s} direction.}
  \item{lambda_u}{the tuning parameter in the \code{u} direction.}
  \item{bd}{the bandwidth of the kernel function.}
  \item{order_1}{order of the B-splines in the \code{s} direction.}
  \item{order_2}{order of the B-splines in the \code{u} direction.}
  \item{breaks}{knots of the B-splines.}
  \item{boot_R}{boostrap times}
  \item{gamma_est_or}{estimate of the scalar parameters using the original estimating equations.}
  \item{b_est_or}{estimate of the functional parameter using the original estimating equations.}
}
\value{
\item{boot_gamma}{a vector of statistics of the bootstrap samples when testing the scalar parameters.}
\item{boot_b}{a vector of statistics of the bootstrap samples when testing the functional parameter.}
\item{gamma_judge}{TRUE or FALSE indicating the significance of the scalar covariates.}
\item{b_judge}{TRUE or FALSE indicating the significance of the functional covariate.}
\item{gamma_pvalue}{p-value of testing the nullity of the scalar parameters.}
\item{b_pvalue}{p-value of testing the nullity of the functional parameters.}
}
\references{
see the paper "Generalized functional partial varying-coefficient model".
}
