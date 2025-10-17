\name{bup}
\alias{bup}
\alias{bup,unmarkedRanef-method}

\title{Extract Best Unbiased Predictors (BUPs) of latent variables from \code{ranef} output}

\description{
  Extract Best Unbiased Predictors (BUPs) of latent variables such as 
  occupancy or abundance from a \code{unmarkedRanef} object produced by the
  \code{ranef} function. Either the posterior mean or median can be requested.
}

\usage{
\S4method{bup}{unmarkedRanef}(object, stat = c("mean", "mode"), ...)
}

\arguments{
  \item{object}{A \code{unmarkedRanef} object, the output of a call to \code{ranef}.}
  \item{stat}{The statistic to calculate from the posterior, either mean or mode.}
  \item{...}{Other arguments, currently ignored} 
}

\value{
  A vector of length M for single-season models, where M is the number of sites.
  A matrix MxT for multi-season models, where T is the number of primary periods.
  Values in the vector or matrix are the BUPs (occupancy or abundance states).
}

\seealso{
\code{\link{ranef}}, \code{\link{unmarkedRanef-class}}
}
