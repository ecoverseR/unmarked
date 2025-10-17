\name{occuComm}
\alias{occuComm}

\title{Fit a community occupancy model with species-level random effects}

\usage{
  occuComm(formula, data, starts, method="BFGS", se=TRUE, ...)
}

\arguments{
  \item{formula}{Double right-hand side formula describing covariates of
    detection and occupancy in that order. It is also possible to include
    random intercepts using \code{lme4}-type syntax.
  }
  \item{data}{An \code{\link{unmarkedFrameOccuComm}} object}
  \item{starts}{Vector of parameter starting values.}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{Logical specifying whether or not to compute standard errors.}
  \item{\dots}{Additional arguments to optim, such as lower and upper bounds}
}

\description{This function fits a multi-species, community occupancy model with 
  random intercepts and slopes by species.}

\details{

In a community occupancy model, detection-nondetection data from multiple species 
are analyzed together. Intercepts and slopes in the model are species-specific and 
come from common distributions, allowing for information sharing across species. 
This structure also allows estimation of site richness. For example, suppose you 
have sites indexed \eqn{i}, occasions \eqn{j} and species \eqn{s}. 
The true occupancy state at a site is

\deqn{z_{is} \sim \mathrm{Bernoulli}(\psi_{is})}

with detection data \eqn{y_{ijs}} modeled as

\deqn{y_{ijs} \sim \mathrm{Bernoulli}(p_{ijs} \cdot z_{is})}

Occupancy probability \eqn{\psi_{is}} can be modeled as a function of covariates 
with species-specific random intercepts and slopes coming from common distributions:

\deqn{\psi_{is} = \mathrm{logit}(\beta_{0,s} + \beta_{i,s} \cdot x_i)}
\deqn{\beta_{0,s} \sim \mathrm{Normal}(\mu_{\beta_0}, \sigma_{\beta_0})}
\deqn{\beta_{1,s} \sim \mathrm{Normal}(\mu_{\beta_1}, \sigma_{\beta_1})}

A similar structure can be implemented for detection probability \eqn{p}.

Note there is a variety of this model that incorporates hypothetical completely 
unobserved species using data augmentation, but that is not a model that 
\code{unmarked} is able to fit at this time.

See the vignette for \code{occuComm} and Kery and Royle (2016), Section 11.6
for more details.
}

\references{
Kery, Marc, and J. Andrew Royle. 2016. \emph{Applied Hierarchical Modeling in
  Ecology}, Volume 1. Academic Press.
}

\value{unmarkedFitOccuComm object describing the model fit.}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarkedFrameOccuComm}}}

\examples{

# Simulate multispecies data
M <- 300
J <- 5
S <- 30

set.seed(123)
x <- rnorm(M)

mu_0 <- 0
sd_0 <- 0.4
beta0 <- rnorm(S, mu_0, sd_0)

mu_x <- 1
sd_x <- 0.3
beta_x <- rnorm(S, mu_x, sd_x)

mu_a <- 0
sd_a <- 0.2
alpha0 <- rnorm(S, mu_a, sd_a)

ylist <- list()
z <- matrix(NA, M, S)
for (s in 1:S){
  psi <- plogis(beta0[s] + beta_x[s] * x)
  z[,s] <- rbinom(M, 1, psi)

  p <- plogis(alpha0[s])

  y <- matrix(NA, M, J)
  for (m in 1:M){
    y[m,] <- rbinom(J, 1, p * z[m,s])
  }
  ylist[[s]] <- y
}
names(ylist) <- paste0("sp", sprintf("\%02d", 1:S))

# Create unmarkedFrame
sc <- data.frame(x=x, a=factor(sample(letters[1:5], M, replace=TRUE)))
umf <- unmarkedFrameOccuComm(ylist, siteCovs=sc)

\donttest{
# Fit model with one covariate on occupancy
(fit <- occuComm(~1~x, umf))

# Look at species-specific random intercepts and slopes
rt <- randomTerms(fit, addMean = TRUE)
head(rt)

# Compare true and estimated random occupancy intercepts
beta0_est <- subset(rt, Model == "psi" & Name == "(Intercept)")$Estimate 
plot(beta0, beta0_est, xlab="truth", ylab="estimate", pch=19, 
     main="Random occ intercepts")
abline(a=0, b=1, col='red')

# Estimate richness for each site
r <- richness(fit)

# Compare true and estimated richness
plot(r, apply(z, 1, sum), xlab="estimate", ylab="truth", main="Richness")
abline(a=0, b=1)
}

}
