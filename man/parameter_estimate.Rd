\name{parameter_estimate}
\alias{parameter_estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the parameters
}
\description{
Esimate both the functiona and scalar parameterzs.
}
\usage{
parameter_estimate(y, x, z, lambda_s, lambda_u, bd, order_1, order_2, breaks,  pois_par, grid)
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
  \item{pois_par}{intensity of the Possion distribution when generating the time points of the functional parameter.}
  \item{grid}{observation points in the \code{u} direction.}
}

\value{
\item{gamma_est}{estimates of the salcar parameters.}
\item{b_est}{estimates of the B-splines coefficients.}
\item{beta_est}{estimate of the functional parameter.}
\item{gamma_stat}{Statistic with respect to testing the nullity of the scalar parameters.}
\item{b_stat}{Statistic with respect to testing the nullity of the functional parameter.}
}
\references{
see the paper "Generalized functional partial varying-coefficient model".
}
