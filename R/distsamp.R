
distsamp <- function(formula, data,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("density", "abund"), unitsOut=c("ha", "kmsq"), starts = NULL,
    method="BFGS", se = TRUE, engine = c("C", "R", "TMB"),
    rel.tol=0.001, ...){

  # Check arguments
  engine <- match.arg(engine)
  formulas <- split_formula(formula)
  if(any(sapply(formulas, has_random))) engine <- "TMB"
  keyfun <- match.arg(keyfun)
  output <- match.arg(output)
  unitsOut <- match.arg(unitsOut)

  #Generate design matrix
  dm <- getDesign(data, formulas)
  y <- dm$y
  M <- nrow(y)
  J <- ncol(y)

  ua <- getUA(data)
  a <- ua$a; u <- ua$u
  A <- get_ds_area(data, unitsOut)
  if(output == "abund") A <- rep(1, length(A))

  # Distance sampling design info
  db <- data@dist.breaks
  tlength <- data@tlength
  survey <- data@survey
  w <- diff(db)
  unitsIn <- data@unitsIn

  # Set up submodels
  state_name <- switch(output, abund = "Abundance", density = "Density")
  estimateList <- unmarkedEstimateList(list(
    state = unmarkedEstimate(state_name, short.name = "lam", invlink="exp")
  ))
  if(keyfun != "uniform") {
    estimateList@estimates$det <- 
      unmarkedEstimate(name = "Detection", short.name = "p", invlink = "exp")
  }
  if(keyfun == "hazard"){
    estimateList@estimates$scale <- 
      unmarkedEstimate(name = "Hazard-rate(scale)", short.name = "p", invlink = "exp")
  }

  # Set up parameter names and indices
  par_inds <- get_parameter_inds(estimateList, dm)

  if(engine=="R") {

    cp <- matrix(NA, M, J)
    switch(keyfun,
    halfnorm = {
        nll <- function(param) {
            sigma <- drop(exp(dm$X_det %*% param[par_inds$det] + dm$offset_det))
            lambda <- drop(exp(dm$X_state %*% param[par_inds$state] + dm$offset_state))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    f.0 <- 2 * dnorm(0, 0, sd=sigma[i])
                    int <- 2 * (pnorm(db[-1], 0, sd=sigma[i]) -
                        pnorm(db[-(J+1)], 0, sd=sigma[i]))
                    cp[i,] <- int / f.0 / w
                    },
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhn, db[j], db[j+1], sigma=sigma[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    exp = {
        nll <- function(param) {
            rate <- drop(exp(dm$X_det %*% param[par_inds$det] + dm$offset_det))
            lambda <- drop(exp(dm$X_state %*% param[par_inds$state] + dm$offset_state))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxexp, db[j], db[j+1], rate=rate[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value / w[j]
                        else {
                            cp[i, j] <- NA
                            }
                        }},
                point = {
                    for(j in 1:J) {
                        int <- integrate(grexp, db[j], db[j+1], rate=rate[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    hazard = {
        nll <- function(param) {
            shape <- drop(exp(dm$X_det %*% param[par_inds$det] + dm$offset_det))
            scale <- drop(exp(param[par_inds$scale]))
            lambda <- drop(exp(dm$X_state %*% param[par_inds$state] + dm$offset_state))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxhaz, db[j], db[j+1], shape=shape[i],
                            scale=scale, stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value / w[j]
                        else {
                            cp[i, j] <- NA
                            }
                        }},
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhaz, db[j], db[j+1], shape=shape[i],
                            scale=scale, stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }

                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    uniform = {
        nll <- function(param) {
            lambda <- drop(exp(dm$X_state %*% param + dm$offset_state))
            if(identical(output, "density"))
                lambda <- lambda * A
            ll <- dpois(y, lambda * u, log=TRUE)
            -sum(ll)
            }
        })
  } else if(engine=="C") {
        nll <- function(param) {
            beta.lam <- param[par_inds$state]
            if(identical(keyfun, "hazard")) {
                beta.sig <- param[par_inds$det]
                scale <- exp(param[par_inds$scale])
            } else {
                beta.sig <- param[par_inds$det]
                scale <- -99.0
            }
            lambda <- drop(exp(dm$X_state %*% beta.lam + dm$offset_state))
            if(identical(output, "density"))
                lambda <- lambda * A
            sigma <- drop(exp(dm$X_det %*% beta.sig + dm$offset_det))
            nll_distsamp(
                  y, lambda, sigma, scale,
                  a, u, w, db,
                  keyfun, survey
            )
        }
  }

  if(engine %in% c("C","R")){
    # This is overly complicated to maintain backwards compatability
    # Maybe we should just set the sigma init to log(median(db)) for all
    if(is.null(starts)){
      starts <- rep(0, max(unlist(par_inds)))
      if(keyfun == "halfnorm"){
        starts[par_inds$det[1]] <- log(max(db))  
      } else if(keyfun == "exp"){
        if(engine == "TMB"){
          starts[par_inds$det[1]] <- log(median(db))
        }
      } else if(keyfun == "hazard"){
        starts[par_inds$det[1]] <- log(median(db))
        starts[par_inds$scale] <- 1
      }
    }

    fit <- fit_optim(nll, starts, method, se, estimateList, par_inds, ...)

  } else if(engine == "TMB"){

    # Set up TMB input data
    survey_type <- switch(survey, line={0}, point={1})
    keyfun_type <- switch(keyfun, uniform={0}, halfnorm={1}, exp={2},
                          hazard={3})
    tmb_inputs <- get_TMB_inputs(formulas = formulas, dm = dm, par_inds = par_inds, umf = data, 
                                 survey_type = survey_type, keyfun_type = keyfun_type,
                                 A = A, db = db, a = a, w = w, u = u)
    if(is.null(starts)){
      if(keyfun != "uniform") tmb_inputs$pars$beta_det[1] <- log(median(db))
    }

    # Fit model with TMB
    fit <- fit_TMB("tmb_distsamp", starts, method, estimateList, par_inds,
                    tmb_inputs, data, ...)
  }

  new("unmarkedFitDS", fitType = "distsamp", call = match.call(),
      opt = fit$opt, formula = formula, formlist = formulas, data = data, keyfun=keyfun,
      sitesRemoved = dm$removed.sites, unitsOut=unitsOut,
      estimates = fit$estimate_list, AIC = fit$AIC, negLogLike = fit$opt$value,
      nllFun = fit$nll, output=output, TMB=fit$TMB)
}


# Detection functions

gxhn <- function(x, sigma) exp(-x^2/(2 * sigma^2))
gxexp <- function(x, rate) exp(-x / rate)
gxhaz <- function(x, shape, scale)  1 - exp(-(x/shape)^-scale)
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
grexp <- function(r, rate) exp(-r / rate) * r
grhaz <- function(r, shape, scale)  (1 - exp(-(r/shape)^-scale)) * r

dxhn <- function(x, sigma)
	gxhn(x=x, sigma=sigma) / integrate(gxhn, 0, Inf, sigma=sigma)$value
drhn <- function(r, sigma)
	grhn(r=r, sigma=sigma) / integrate(grhn, 0, Inf, sigma=sigma)$value
dxexp <- function(x, rate)
	gxexp(x=x, rate=rate) / integrate(gxexp, 0, Inf, rate=rate)$value
drexp <- function(r, rate)
	grexp(r=r, rate=rate) / integrate(grexp, 0, Inf, rate=rate)$value
dxhaz <- function(x, shape, scale)
	gxhaz(x=x, shape=shape, scale=scale) / integrate(gxhaz, 0, Inf,
		shape=shape, scale=scale)$value
drhaz <- function(r, shape, scale)
	grhaz(r=r, shape=shape, scale=scale) / integrate(grhaz, 0, Inf,
		shape=shape, scale=scale)$value


