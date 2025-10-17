\name{powerAnalysis}
\alias{powerAnalysis}
\alias{powerAnalysis,unmarkedFrame-method}
\alias{powerAnalysis,list-method}
\alias{powerAnalysis,unmarkedFit-method}

\title{Conduct a power analysis for an \code{unmarked} model}

\description{
This function uses a simulation-based approach to estimate power for parameters 
in unmarked models. At a minimum, users must provide an \code{unmarkedFrame}
object describing the experimental design and a list of effect sizes for each 
parameter in the model. See the \code{unmarkedPower} vignette for more details and
examples.
}

\usage{
  \S4method{powerAnalysis}{unmarkedFrame}(object, model = NULL, effects = NULL,
    alpha=0.05, nsim = 100, parallel = FALSE, nulls=NULL, ...)
  \S4method{powerAnalysis}{list}(object, model = NULL, effects = NULL,
    alpha=0.05, nsim = length(object), parallel = FALSE, nulls=NULL, ...)
}

\arguments{
  \item{object}{An \code{unmarkedFrame} object representing the desired study
    design. The values in the response (\code{y}) don't matter and can be missing.
    Alternatively, you can provide a list of such objects with the 
    response data already simulated (such as the output from \code{simulate}).
  }
  \item{model}{The model to use when the \code{unmarkedFrame} type is used for
    multiple model types. For example, if the \code{object} is an
    \code{unmarkedFrameOccu}, model should be set to \code{occu} or \code{occuRN}.
  }
  \item{effects}{A list containing the desired effect sizes/parameter values
    for which you want to estimate power. This list must follow a specific format. 
    There is one named entry in the list per submodel (e.g., occupancy, detection). 
    Each list element should be a numeric vector with length equal to the number
    of parameters in that submodel. Parameter values are on the inverse link
    scale. You can leave \code{effects=NULL}, which
    will generate an error message with a template that you can fill in.
  }
  \item{alpha}{Desired Type I error rate.}
  \item{nsim}{Number of simulations to conduct.}
  \item{parallel}{Logical; run simulations in parallel?}
  \item{nulls}{If provided, a list matching the structure of \code{effects} which
    defines the null hypothesis value for each parameter. By default the null
    is 0 for all parameters.
  }
  \item{...}{Arguments to send to the fitting function for the model. 
    Most importantly this will include formula argument(s), but could also include 
    distributions, key functions, etc. For example, for simulating occupancy data, 
    you must also supply the argument \code{formula = ~1~1} for a no-covariate model, 
  \code{formula=~1~x} for a covariate effect of \code{x} on occupancy, etc.}
} 

\value{\code{unmarkedPower} object containing the results of the power analysis.
      For information on interpretation of the output, see the power analysis
      vignette.}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{unmarkedPowerList}}
}

\examples{

\dontrun{

# Create experimental design
M <- 50
J <- 3
y <- matrix(NA, M, J)
sc <- data.frame(x=rnorm(M))
umf <- unmarkedFrameOccu(y, siteCovs=sc)

# Power analysis
p <- powerAnalysis(umf, model=occu, formula=~1~x, 
                   effects = list(state = c(-0.2, 0.3), det = 0))

p
summary(p, alpha=0.3)
plot(p, ylim=c(-3, 3))
plot(p, ylim=c(-3, 3))

# Simulate your own datasets first and pass to power analysis
cf <- list(state=c(0,1), det=0)
s <- simulate(umf, model = occu, formula=~1~x, coefs=cf, nsim = 100)
p2 <- powerAnalysis(s, model=occu, formula=~1~x, effects=cf)
p2

}
}
