
gdistsamp <- function(lambdaformula, phiformula, pformula, data,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("abund", "density"), unitsOut=c("ha", "kmsq"),
    mixture=c("P", "NB", 'ZIP'), K, starts = NULL, method = "BFGS", se = TRUE, 
    engine=c("C","R"), rel.tol=1e-4, threads=1, ...){
  
  if(!is(data, "unmarkedFrameGDS"))
      stop("Data is not of class unmarkedFrameGDS.")

  engine <- match.arg(engine, c("C", "R"))
  keyfun <- match.arg(keyfun)
  if(!keyfun %in% c("halfnorm", "exp", "hazard", "uniform"))
    stop("keyfun must be 'halfnorm', 'exp', 'hazard', or 'uniform'")
  output <- match.arg(output)
  unitsOut <- match.arg(unitsOut)
  db <- data@dist.breaks
  w <- diff(db)
  tlength <- data@tlength
  survey <- data@survey
  unitsIn <- data@unitsIn
  mixture <- match.arg(mixture)

  formulas <- list(lambda = lambdaformula, phi = phiformula, det = pformula)
  check_no_support(formulas)
  D <- getDesign(data, formulas)
  y <- D$y  # MxJT

  M <- nrow(y)
  T <- data@numPrimary
  R <- ncol(y)
  J <- R / T

  y <- array(y, c(M, J, T))
  y <- aperm(y, c(1,3,2))
  yt <- apply(y, 1:2, function(x) {
    if(all(is.na(x)))
        return(NA)
    else return(sum(x, na.rm=TRUE))
    })

  minK <- max(yt, na.rm=TRUE)
  if(missing(K) || is.null(K)){
    K <- minK + 100
  } else {
    if(K < minK){
      stop("K should be larger than the max observed abundance, ", minK, call.=FALSE)
    }
  }
  k <- 0:K
  lk <- length(k)

  ua <- getUA(data)
  a <- ua$a; u <- ua$u
  A <- get_ds_area(data, unitsOut)
  if(output == "abund") A <- rep(1, length(A))
  if(length(D$removed.sites) > 0){
    a <- a[-D$removed.sites,,drop=FALSE]
    u <- u[-D$removed.sites,,drop=FALSE]
    A <- A[-D$removed.sites]
  }

  # Set up submodels
  estimateList <- unmarkedEstimateList(list(
    lambda = unmarkedEstimate(name = "Abundance", short.name="lambda", invlink="exp")
  ))

  if(T>1){
    estimateList@estimates$phi <- 
      unmarkedEstimate(name = "Availability", short.name = "phi", invlink = "logistic")
  }

  if(keyfun != "uniform"){
    estimateList@estimates$det <- 
      unmarkedEstimate(name = "Detection", short.name = "p", invlink = "exp")
  }

  if(keyfun == "hazard"){
    estimateList@estimates$scale <- 
      unmarkedEstimate(name = "Hazard-rate(scale)", short.name = "scale", invlink = "exp")
  }

  if(mixture == "NB"){
    estimateList@estimates$alpha <- 
      unmarkedEstimate(name = "Dispersion", short.name = "alpha", invlink = "exp")
  } else if(mixture == "ZIP"){
    estimateList@estimates$psi <- 
      unmarkedEstimate(name="Zero-inflation", short.name = "psi", invlink = "logistic")
  }

  # Set up parameter names and indices
  par_inds <- get_parameter_inds(estimateList, D)
  
  cp <- array(as.numeric(NA), c(M, T, J+1))
  g <- matrix(as.numeric(NA), M, lk)

  lfac.k <- lgamma(k+1)
  kmyt <- array(NA, c(M, T, lk))
  lfac.kmyt <- array(0, c(M, T, lk))
  fin <- matrix(NA, M, lk)
  naflag <- array(NA, c(M, T, J))
  for(i in 1:M) {
    fin[i, ] <- k - max(yt[i,], na.rm=TRUE) >= 0
    for(t in 1:T) {
      naflag[i,t,] <- is.na(y[i,t,])
      if(!all(naflag[i,t,])) {
        kmyt[i,t,] <- k - yt[i,t]
        lfac.kmyt[i, t, fin[i,]] <- lgamma(kmyt[i, t, fin[i,]] + 1)
      }
    }
  }

  switch(keyfun,
    halfnorm = {
      altdetParms <- paste("sigma", colnames(D$X_det), sep="")
      if(missing(starts)) {
        starts <- rep(0, max(unlist(par_inds)))
        starts[par_inds$det[1]] <- log(max(db))
      }

      nll_R <- function(pars) {
        lambda <- exp(D$X_lambda %*% pars[par_inds$lambda] + D$offset_lambda)
        if(identical(output, "density"))
            lambda <- lambda * A

        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(D$X_phi %*% pars[par_inds$phi] + D$offset_phi)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        sigma <- exp(D$X_det %*% pars[par_inds$det]+D$offset_det)
        sigma <- matrix(sigma, M, T, byrow=TRUE)

        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[par_inds$alpha]))),
            ZIP = f <- sapply(k, function(x) dzip(rep(x, length(lambda)), lambda=lambda, psi=plogis(pars[par_inds$psi])))
        )
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(any(naflag[i,t,]))
                    next
                p <- rep(NA, J)
                switch(survey,
                line = {
                    f.0 <- 2 * dnorm(0, 0, sd=sigma[i, t])
                    int <- 2 * (pnorm(db[-1], 0, sd=sigma[i, t]) -
                        pnorm(db[-(J+1)], 0, sd=sigma[i, t]))
                    p <- int / f.0 / w
                    },
                point = {
                    for(j in 1:J) {
#                        int <- integrate(grhn, db[j], db[j+1],
#                            sigma=sigma[i, t], rel.tol=rel.tol,
#                            stop.on.error=FALSE, subdivisions=50)
#                        if(!identical(int$message, "OK"))
#                            int$value <- NA
                        int <- sigma[i,t]^2 *
                            (1-exp(-db[j+1]^2 / (2*sigma[i,t]^2))) -
                                sigma[i,t]^2 *
                                    (1-exp(-db[j]^2 / (2*sigma[i,t]^2)))
#                        p[j] <- int$value * 2 * pi / a[i,j]
                        p[j] <- int * 2 * pi / a[i,j]
                        }
                    })
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)
                if(!is.na(cp[J+1]) && cp[J+1]<0){ cp[J+1] <- 0 }

                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, ] *
                    log(cp[1:J])) +
                    kmyt[i, t,] * log(cp[J+1])
            }
            g[i,] <- exp(rowSums(mn))
        }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
        }
    },
    exp = {
      altdetParms <- paste("rate", colnames(D$X_det), sep="")
      if(missing(starts)) {
        starts <- rep(0, max(unlist(par_inds)))
        starts[par_inds$det[1]] <- log(max(db))
      }

      nll_R <- function(pars) {
        lambda <- exp(D$X_lambda %*% pars[par_inds$lambda] + D$offset_lambda)
        if(identical(output, "density"))
            lambda <- lambda * A
        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(D$X_phi %*% pars[par_inds$phi] + D$offset_phi)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        rate <- exp(D$X_det %*% pars[par_inds$det] + D$offset_det)
        rate <- matrix(rate, M, T, byrow=TRUE)

        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[par_inds$alpha]))),
            ZIP = f <- sapply(k, function(x) dzip(rep(x, length(lambda)), lambda=lambda, psi=plogis(pars[par_inds$psi])))
        )
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(any(naflag[i,t,]))
                    next
                p <- rep(NA, J)
                switch(survey,
                line = {
                    for(j in 1:J) {
#                        int <- integrate(gxexp, db[j], db[j+1],
#                             rate=rate[i,t], rel.tol=rel.tol,
#                             stop.on.error=FALSE, subdivisions=50)
#                        if(!identical(int$message, "OK"))
#                            int$value <- NA
                        int <- rate[i,t]*(1-exp(-db[j+1]/rate[i,t])) -
                            rate[i,t]*(1-exp(-db[j]/rate[i,t]))
#                        p[j] <- int$value / w[j]
                        p[j] <- int / w[j]
                        }
                    },
                point = {
                    for(j in 1:J) {
                        int <- integrate(grexp, db[j], db[j+1],
                            rate=rate[i, t], rel.tol=rel.tol,
                            stop.on.error=FALSE, subdivisions=50)
                        if(!identical(int$message, "OK"))
                            int$value <- NA
                        p[j] <- int$value * 2 * pi / a[i,j]
                        }
                    })
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)
                if(!is.na(cp[J+1]) && cp[J+1]<0){ cp[J+1] <- 0 }
                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, ] *
                    log(cp[1:J])) +
                    kmyt[i, t,] * log(cp[J+1])
            }
            g[i,] <- exp(rowSums(mn))
        }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
        }
    },
    hazard = {
      altdetParms <- paste("shape", colnames(D$X_det), sep="")
      if(missing(starts)) {
        starts <- rep(0, max(unlist(par_inds)))
      }
      nll_R <- function(pars) {
        lambda <- exp(D$X_lambda %*% pars[par_inds$lambda] + D$offset_lambda)
        if(identical(output, "density"))
            lambda <- lambda * A

        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(D$X_phi %*% pars[par_inds$phi] + D$offset_phi)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        shape <- exp(D$X_det %*% pars[par_inds$det]+D$offset_det)
        shape <- matrix(shape, M, T, byrow=TRUE)

        scale <- exp(pars[par_inds$scale])

        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[par_inds$alpha]))),
            ZIP = f <- sapply(k, function(x) dzip(rep(x, length(lambda)), lambda=lambda, psi=plogis(pars[par_inds$psi])))
        )
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(any(naflag[i,t,]))
                    next
                p <- rep(NA, J)
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxhaz, db[j], db[j+1],
                             shape=shape[i,t], scale=scale,
                             rel.tol=rel.tol,
                             stop.on.error=FALSE, subdivisions=50)
                        if(!identical(int$message, "OK"))
                            int$value <- NA
                        p[j] <- int$value / w[j]
                        }
                    },
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhaz, db[j], db[j+1],
                            shape=shape[i, t], scale=scale,
                            rel.tol=rel.tol,
                            stop.on.error=FALSE, subdivisions=50)
                        if(!identical(int$message, "OK"))
                            int$value <- NA
                        p[j] <- int$value * 2 * pi / a[i,j]
                        }
                    })
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)
                if(!is.na(cp[J+1]) && cp[J+1]<0){ cp[J+1] <- 0 }
                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, ] *
                    log(cp[1:J])) +
                    kmyt[i, t,] * log(cp[J+1])
                }
            g[i,] <- exp(rowSums(mn))
            }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
        }
    },
    uniform = {
      if(missing(starts)) {
        starts <- rep(0, max(unlist(par_inds)))
      }
      nll_R <- function(pars) {
        lambda <- exp(D$X_lambda %*% pars[par_inds$lambda] + D$offset_lambda)
        if(identical(output, "density"))
            lambda <- lambda * A
        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(D$X_phi %*% pars[par_inds$phi] + D$offset_phi)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        p <- 1
        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[par_inds$alpha]))),
            ZIP = f <- sapply(k, function(x) dzip(rep(x, length(lambda)), lambda=lambda, psi=plogis(pars[par_inds$psi])))
        )
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(any(naflag[i,t,]))
                    next
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)
                if(!is.na(cp[J+1]) && cp[J+1]<0){ cp[J+1] <- 0 }
                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, ] *
                    log(cp[1:J])) +
                    kmyt[i, t,] * log(cp[J+1])
            }
            g[i,] <- exp(rowSums(mn))
        }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
    }
  })

  if(engine =="C"){
    long_format <- function(x){
      out <- matrix(aperm(x,c(1,3,2)),nrow=nrow(x),ncol=dim(x)[2]*dim(x)[3])
      as.vector(t(out))
    }
    y_long <- long_format(y)
    # Vectorize these arrays as using arma::subcube sometimes crashes
    kmytC <- as.vector(aperm(kmyt, c(3,2,1)))
    lfac.kmytC <- as.vector(aperm(lfac.kmyt, c(3,2,1)))
    if(output!='density'){
      A <- rep(1, M)
    }
    mixture_code <- switch(mixture, P={1}, NB={2}, ZIP={3})
    # TODO: use more consistent pattern here
    n_param <- c(length(par_inds$lambda), length(par_inds$phi), length(par_inds$det),
                 length(par_inds$scale), length(par_inds$alpha)+length(par_inds$psi))
    Kmin <- apply(yt, 1, max, na.rm=TRUE)

    nll <- function(params){
      nll_gdistsamp(params, n_param, y_long, mixture_code, keyfun, survey,
                    D$X_lambda, D$offset_lambda, A, D$X_phi, D$offset_phi, D$X_det, D$offset_det,
                    db, a, t(u), w, k, lfac.k, lfac.kmytC, kmytC, Kmin, threads)
    }

  } else {
    nll <- nll_R
  }

  fit <- fit_optim(nll, starts, method, se, estimateList, par_inds, ...)
  
  new("unmarkedFitGDS", fitType = "gdistsamp", call = match.call(), 
      formlist = formulas, data = data, estimates = fit$estimate_list, 
      sitesRemoved = D$removed.sites, AIC = fit$AIC, opt = fit$opt, 
      negLogLike = fit$opt$value, nllFun = fit$nll,
      mixture=mixture, K=K, keyfun=keyfun, unitsOut=unitsOut, output=output)
}
