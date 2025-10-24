\name{getStarts}
\alias{getStarts}
\title{Get order of starting values for an unmarked model}
\description{
  This function provides a named vector of starting values for a given
  unmarked model, to serve as a guide for providing your own starting values.
  To use it, write the code to fit an unmarked model and then wrap the code
  in \code{getStarts()} (see example below). Note that due to environment/scoping
  issues, this function may not work if you call it from inside other functions.
}
\usage{
getStarts(code)
}
\arguments{
  \item{code}{A line of R code creating an \code{unmarkedFit} object}
}
\value{A named vector showing the correct order of starting values}
\examples{
# Create unmarkedFrame
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))

# Get correct sequence of starting values for model
getStarts(occu(~ obsvar1 ~ 1, pferUMF))
}
