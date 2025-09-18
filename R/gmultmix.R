
# data will need to be an unmarkedMultFrame
gmultmix <- function(lambdaformula, phiformula, pformula, data,
    mixture=c("P", "NB", "ZIP"), K, starts = NULL, method = "BFGS", se = TRUE,
    engine=c("C","R"), threads=1, ...){
  
  if(!is(data, "unmarkedFrameGMM"))
    stop("Data is not of class unmarkedFrameGMM.")

  engine <- match.arg(engine, c("C", "R"))
  mixture <- match.arg(mixture)

  formulas <- list(lambda = lambdaformula, phi = phiformula, det = pformula)
  check_no_support(formulas)
  D <- getDesign(data, formulas)
  y <- D$y  # MxJT

  K <- check_K_multinomial(K, K_adjust = 100, y, data@numPrimary)
  k <- 0:K
  lk <- length(k)
  M <- nrow(y)
  T <- data@numPrimary
  R <- numY(data) / T
  J <- obsNum(data) / T

  y <- array(y, c(M, R, T))
  y <- aperm(y, c(1,3,2))
  yt <- apply(y, 1:2, function(x) {
    if(all(is.na(x)))
        return(NA)
    else return(sum(x, na.rm=TRUE))
  })

  piFun <- data@piFun

  # Set up submodels
  estimateList <- unmarkedEstimateList(list(
    lambda = unmarkedEstimate(name = "Abundance", short.name="lambda", invlink="exp")
  ))

  if(T>1){
    estimateList@estimates$phi <- 
      unmarkedEstimate(name = "Availability", short.name = "phi", invlink = "logistic")
  }

  estimateList@estimates$det <- 
    unmarkedEstimate(name = "Detection", short.name = "p", invlink = "logistic")

  if(mixture == "NB"){
    estimateList@estimates$alpha <- 
      unmarkedEstimate(name = "Dispersion", short.name = "alpha", invlink = "exp")
  } else if(mixture == "ZIP"){
    estimateList@estimates$psi <- 
      unmarkedEstimate(name="Zero-inflation", short.name = "psi", invlink = "logistic")
  }

  # Set up parameter names and indices
  par_inds <- get_parameter_inds(estimateList, D)
  

  lfac.k <- lgamma(k+1)
  kmyt <- array(NA, c(M, T, lk))
  lfac.kmyt <- array(0, c(M, T, lk))
  fin <- matrix(NA, M, lk)
  naflag <- array(NA, c(M, T, R))
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

  nll_R <- function(pars) {
    lambda <- exp(D$X_lambda %*% pars[par_inds$lambda] + D$offset_lambda)
    if(T==1)
        phi <- 1
    else if(T>1)
        phi <- drop(plogis(D$X_phi %*% pars[par_inds$phi] + D$offset_phi))
    p <- plogis(D$X_det %*% pars[par_inds$det] + D$offset_det)

    phi.mat <- matrix(phi, M, T, byrow=TRUE)
    phi <- as.numeric(phi.mat)

    p <- matrix(p, nrow=M, byrow=TRUE)
    p <- array(p, c(M, J, T))
    p <- aperm(p, c(1,3,2))
    cp <- array(as.numeric(NA), c(M, T, R+1))

    for(t in 1:T) cp[,t,1:R] <- do.call(piFun, list(p[,t,]))
    cp[,,1:R] <- cp[,,1:R] * phi
    cp[,, 1:R][is.na(y)]<- NA   # andy added 5/29
    cp[,,R+1] <- 1 - apply(cp[,,1:R,drop=FALSE], 1:2, sum, na.rm=TRUE)

    switch(mixture,
      P = f <- sapply(k, function(x) dpois(x, lambda)),
      NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[par_inds$alpha]))),
      ZIP = f <- sapply(k, function(x) dzip(rep(x, length(lambda)), lambda=lambda, psi=plogis(pars[par_inds$psi])))
    )
    g <- matrix(as.numeric(NA), M, lk)
    for(i in 1:M) {
        A <- matrix(0, lk, T)
        for(t in 1:T) {
            na <- naflag[i,t,]
            if(!all(na))
                A[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, !na] * log(cp[i, t, which(!na)])) +
                    kmyt[i, t,] * log(cp[i, t, R+1])
            }
        g[i,] <- exp(rowSums(A))
        }
    f[!fin] <- g[!fin] <- 0
    ll <- rowSums(f*g)
    -sum(log(ll))
  }


  if(engine=="R"){
    nll <- nll_R
  } else {
    long_format <- function(x){
      out <- matrix(aperm(x,c(1,3,2)),nrow=nrow(x),ncol=dim(x)[2]*dim(x)[3])
      as.vector(t(out))
    }
    y_long <- long_format(y)

    # Vectorize these arrays as using arma::subcube sometimes crashes
    kmytC <- as.vector(aperm(kmyt, c(3,2,1)))
    kmytC[which(is.na(kmytC))] <- 0
    lfac.kmytC <- as.vector(aperm(lfac.kmyt, c(3,2,1)))
 
    mixture_code <- switch(mixture, P={1}, NB={2}, ZIP={3})
    n_param <- c(length(par_inds$lambda), length(par_inds$phi),
                 length(par_inds$det), mixture%in%c("NB","ZIP"))
    Kmin <- apply(yt, 1, max, na.rm=TRUE)

    nll <- function(params) {
      nll_gmultmix(params, n_param, y_long, mixture_code, piFun, D$X_lambda, D$offset_lambda,
                   D$X_phi, D$offset_phi, D$X_det, D$offset_det, k, lfac.k, lfac.kmytC,
                   kmytC, Kmin, threads)
    }

    if(!piFun%in%c('doublePiFun','removalPiFun','depDoublePiFun')){
      warning("Custom pi functions are not supported by C engine. Using R engine instead.")
      nll <- nll_R
    }
  }

  fit <- fit_optim(nll, starts, method, se, estimateList, par_inds, ...)
 
  new("unmarkedFitGMM", fitType = "gmn", call = match.call(), 
      formlist = formulas, data = data, estimates = fit$estimate_list, 
      sitesRemoved = D$removed.sites, AIC = fit$AIC, opt = fit$opt, 
      negLogLike = fit$opt$value, nllFun = fit$nll, mixture=mixture, K=K)
}
