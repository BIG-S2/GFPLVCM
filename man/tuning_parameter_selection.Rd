\name{tuning_parameter_selection}
\alias{tuning_parameter_selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select the tuning parameters
}
\description{
Select the tuning parameters via the roughly estimated mean square errors and do estimation, prediction as well as hypothesis testing with the selected tuning parameters.
}
\usage{
tuning_parameter_selection( lambda_range, h_range, y, x, z, order_1, order_2, breaks, pois_par, grid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda_range}{a vector of the tuning parameter lambda.}
  \item{h_range}{a vector of the bandwidth.}
  \item{y}{a list of the response with the elements \code{y_ID}, \code{y_time_point} and \code{y_value}.}
  \item{x}{a list of the functional covariate with the elements \code{x_ID}, \code{x_time_point} and \code{x_value}.}
  \item{z}{a list of the scalar covariate with the elements \code{z_ID}, \code{z_time_point} and \code{z_value}.}
  \item{order_1}{order of the B-splines in the \code{s} direction.}
  \item{order_2}{order of the B-splines in the \code{u} direction.}
  \item{breaks}{knots of the B-splines.}
  \item{pois_par}{intensity of the Possion distribution when generating the time points of the functional parameter.}
  \item{grid}{observation points in the \code{u} direction.}
}

\value{
\item{final_lambda}{The seleted roughness penality tuning parameter.}
\item{final_bd}{The selected bandwidth.}
}
\references{
see the paper "Generalized functional partial varying-coefficient model".
}

