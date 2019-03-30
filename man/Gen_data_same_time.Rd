\name{Gen_data_same_time}
\alias{Gen_data_same_time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the data at the same time points
}
\description{
Generate the response and the covariates at the same time points.
}
\usage{
Gen_data_same_time(n, pois_par, u_time, gamma, len_k, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{n}{The number of sample size to be generated.}
  \item{pois_par}{The number of time points in the \code{t} direction follows from a Possion(\code{pois_par}).}
  \item{u_time}{length of the observation time in the u direction.}
  \item{gamma}{true values of the scalar parameters}
  \item{len_k}{number of the basis function used when generating the functional parameter.}
  \item{B}{signal of the functional parameter.}
}

\value{
\item{y}{The response in the simulations.}
\item{x}{The functional variable in the simulations with the same observation times with the response.}
\item{z}{The scalar variables in the simulations with the same observation times with the response.}
}

\references{
see the paper "Generalized functional partial varying-coefficient model".
}
