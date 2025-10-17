\name{[-methods}
\docType{methods}
\alias{[-methods}
\alias{[,unmarkedEstimateList,ANY,ANY,ANY-method}
\alias{[,unmarkedFit,ANY,ANY,ANY-method}
\alias{[,unmarkedFrame,numericOrLogical,missing,missing-method}
\alias{[,unmarkedFrame,missing,numericOrLogical,missing-method}
\alias{[,unmarkedFrame,numericOrLogical,numericOrLogical,missing-method}
\title{Methods for bracket extraction [ in Package `unmarked'}
\description{
Methods for bracket extraction \code{[} in Package `unmarked'
}
\usage{
\S4method{[}{unmarkedEstimateList,ANY,ANY,ANY}(x, i, j, drop)
\S4method{[}{unmarkedFit,ANY,ANY,ANY}(x, i, j, drop)
\S4method{[}{unmarkedFrame,numericOrLogical,numericOrLogical,missing}(x, i, j)
}
\arguments{
	\item{x}{Object of appropriate S4 class}
	\item{i}{Site numbers or a vector of TRUE/FALSE}
	\item{j}{Observation numbers (occasions or primary periods) or a vector of TRUE/FALSE}
	\item{drop}{Not currently used}
}
\section{Methods}{
\describe{

\item{x = "unmarkedEstimateList", i = "ANY", j = "ANY", drop = "ANY"}{Extract a
	unmarkedEstimate object from an unmarkedEstimateList by name (either 'det'
	or 'state')}

\item{x = "unmarkedFit", i = "ANY", j = "ANY", drop = "ANY"}{Extract a
	unmarkedEstimate object from an unmarkedFit by name (either 'det'
	or 'state')}

\item{x = "unmarkedFrame", i = "numericOrLogical", j = "numericOrLogical", drop = "missing"}{
	Extract sites and/or observations from an unmarkedFrame }
}}


\examples{

data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
	obsCovs = mallard.obs)
summary(mallardUMF)

mallardUMF[1:5,]
mallardUMF[,1:2]
mallardUMF[1:5, 1:2]

}
\keyword{methods}

