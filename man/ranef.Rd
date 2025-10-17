\name{ranef}
\alias{ranef}
\alias{ranef,unmarkedFit-method}

\title{Estimate posterior distributions of latent occupancy or abundance}
\description{
Estimate posterior distributions of the random variables (latent
occupancy or abundance) using empirical Bayes methods. These methods
return an object storing the posterior distributions of the latent
variables at each site, and for each year (primary period) in the case
of open population models. See \link{unmarkedRanef-class} and \link{posteriorSamples}
for methods used to manipulate the returned object.
}

\usage{
\S4method{ranef}{unmarkedFit}(object, ...)
}

\arguments{
  \item{object}{A \code{unmarkedFit} object}
  \item{...}{Other arguments. For some abundance models you can specify \code{K}; 
  for multi-species models you can specify \code{species}.}
}

\note{
  From Carlin and Louis (1996): \dQuote{... the Bayesian approach to
  inference depends on a prior distribution for the model
  parameters. This prior can depend on unknown parameters which in turn
  may follow some second-stage prior. This sequence of parameters and
  priors consitutes a hierarchical model. The hierarchy must stop at
  some point, with all remaining prior parameters assumed known. Rather
  than make this assumption, the basic empirical Bayes approach uses the
  observed data to estimate these final stage parameters (or to estimate
  the Bayes rule), and proceeds as in a standard Bayesian analysis.}
}


\section{Warning}{
  Empirical Bayes methods can underestimate the variance of the
  posterior distribution because they do not account for uncertainty in
  the hyperparameters (lambda or psi). Eventually, we hope to add
  methods to account for the uncertainty of the hyperparameters.

  Note also that the posterior mode appears to exhibit some bias as an
  estimator or abundance. Consider using the posterior mean instead,
  even though it will not be an integer in general. More
  simulation studies are needed to evaluate the performance of empirical
  Bayes methods for these models.
}

\author{Richard Chandler \email{rbchan@uga.edu}}

\seealso{
\code{\link{unmarkedRanef-class}}, \code{\link{bup}}
}

\references{
  Laird, N.M. and T.A. Louis. 1987. Empirical Bayes confidence intervals
  based on bootstrap samples. Journal of the American Statistical
  Association 82:739--750.

  Carlin, B.P and T.A Louis. 1996. Bayes and Empirical Bayes Methods for
  Data Analysis. Chapman and Hall/CRC.

  Royle, J.A and R.M. Dorazio. 2008. Hierarchical Modeling and Inference
  in Ecology. Academic Press.
}

\examples{
# Simulate data under N-mixture model
set.seed(4564)
R <- 20
J <- 5
N <- rpois(R, 10)
y <- matrix(NA, R, J)
y[] <- rbinom(R*J, N, 0.5)

# Fit model
umf <- unmarkedFramePCount(y=y)
fm <- pcount(~1 ~1, umf, K=50)

# Estimates of conditional abundance distribution at each site
(re <- ranef(fm))
# Best Unbiased Predictors
bup(re, stat="mean")           # Posterior mean
bup(re, stat="mode")           # Posterior mode
confint(re, level=0.9) # 90\% CI

# Plots
plot(re, subset=site \%in\% c(1:10), layout=c(5, 2), xlim=c(-1,20))

# Compare estimates to truth
sum(N)
sum(bup(re))

# Extract all values in convenient formats
post.df <- as(re, "data.frame")
head(post.df)
post.arr <- as(re, "array")

#Generate posterior predictive distribution for a function
#of random variables using predict()

#First, create a function that operates on a vector of 
#length M (if you fit a single-season model) or a matrix of 
#dimensions MxT (if a dynamic model), where
#M = nsites and T = n primary periods
#Our function will generate mean abundance for sites 1-10 and sites 11-20
myfunc <- function(x){ #x will be length 20 since M=20
  
  #Mean of first 10 sites
  group1 <- mean(x[1:10])
  #Mean of sites 11-20
  group2 <- mean(x[11:20])
  
  #Naming elements of the output is optional but helpful
  return(c(group1=group1, group2=group2))

}

#Get 100 samples of the values calculated in your function
(pr <- predict(re, func=myfunc, nsims=100))

#Summarize posterior
data.frame(mean=rowMeans(pr),
           se=apply(pr, 1, stats::sd),
           lower=apply(pr, 1, stats::quantile, 0.025),
           upper=apply(pr, 1, stats::quantile, 0.975))

#Alternatively, you can return the posterior predictive distribution
#and run operations on it separately
(ppd <- posteriorSamples(re, nsims=100))

}
