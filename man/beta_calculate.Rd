\name{beta_calculate}
\alias{beta_calculate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the true value of the functional parameter
}
\description{
Calculate the true value of the functional parameter beta when given the observation points of the two directions in the simulations.
}
\usage{
beta_calculate(s_time_point, u_time, len_k, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s_time_point}{a vector of the observation times in the s direction.}
  \item{u_time}{length of the observation time in the \code{u} direction.}
  \item{len_k}{number of the basis function used when generating the functional parameter.}
  \item{B}{signal of the functional parameter.}
}

\value{
\code{beta_calculate()} returns the true value of the functional parameter given the observation times.
}
\references{
see the paper "Generalized functional partial varying-coefficient model".
}
