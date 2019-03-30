\name{Gen_x_z}
\alias{Gen_x_z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the covariates
}
\description{
Generate the scalar and the functional covariates for one subject at given time points in the simulations.
}
\usage{
Gen_x_z(time_number, time_point, u_time, len_k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time_number}{the number of the observation time points for the subject.}
  \item{time_point}{a vector of the time points.}
  \item{u_time}{length of the observation time in the \code{u} direction.}
  \item{len_k}{number of the basis function used when generating the functional parameter.}
}

\value{
\item{z_temp_value}{values of the scalar covariates at the observation time points.}
\item{x_temp_value}{values of the functional covariate at the observation time points.}
\item{beta_temp_value}{values of the functional parameter at the observation time points.}
}
\references{
see the paper "Generalized functional partial varying-coefficient model".
}
