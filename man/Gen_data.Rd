\name{Gen_data}
\alias{Gen_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the data in the simulations
}
\description{
Generate the data used in the simulations.
}
\usage{
Gen_data(n, pois_par, u_time, gamma, len_k, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of sample size to be generated.}
  \item{pois_par}{number of time points in the \code{t} direction follows from a Possion(\code{pois_par}).}
  \item{u_time}{length of the observation time in the u direction.}
  \item{gamma}{true values of the scalar parameters}
  \item{len_k}{number of the basis function used when generating the functional parameter.}
  \item{B}{true signal of the functional parameter.}
}

\value{
\item{y}{the response in the simulations.}
\item{x}{the functional variable in the simulations.}
\item{z}{the scalar variables in the simulations.}
}

\references{
see the paper "Generalized functional partial varying-coefficient model".
}

