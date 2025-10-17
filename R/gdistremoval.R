

unmarkedFrameGDR <- function(yDistance, yRemoval, numPrimary=1,
                                     siteCovs=NULL, obsCovs=NULL,
                                     yearlySiteCovs=NULL, dist.breaks,
                                     unitsIn, period.lengths=NULL){

  if(is.null(period.lengths)){
    period.lengths <- rep(1, ncol(yRemoval)/numPrimary)
  }

  M <- nrow(yDistance)
  Jdist <- ncol(yDistance) / numPrimary
  Jrem <- ncol(yRemoval) / numPrimary

  if(length(dist.breaks) != Jdist+1){
    stop(paste("dist.breaks must have length",Jdist+1), call.=FALSE)
  }
  if(length(period.lengths) != Jrem){
    stop(paste("period.lengths must have length",Jrem), call.=FALSE)
  }

  dist_array <- array(as.vector(t(yDistance)), c(Jdist, numPrimary, M))
  dist_sums <- apply(dist_array, c(2,3), sum, na.rm=T)

  rem_array <- array(as.vector(t(yRemoval)), c(Jrem, numPrimary, M))
  rem_sums <- apply(rem_array, c(2,3), sum, na.rm=T)

  if(!all(dist_sums == rem_sums)){
    stop("Some sites/primary periods do not have the same number of distance and removal observations", call.=FALSE)
  }

  umf <- new("unmarkedFrameGDR", y=yRemoval, yDistance=yDistance,
             yRemoval=yRemoval, numPrimary=numPrimary, siteCovs=siteCovs,
             obsCovs=obsCovs, yearlySiteCovs=yearlySiteCovs, survey="point",
             dist.breaks=dist.breaks, unitsIn=unitsIn, period.lengths=period.lengths,
             obsToY=diag(ncol(yRemoval)))
  umf <- umf_to_factor(umf)
  umf
}

setAs("unmarkedFrameGDR", "data.frame", function(from){

  out <- callNextMethod(from, "data.frame")
  J <- obsNum(from)
  out <- out[,(J+1):ncol(out), drop=FALSE]

  yDistance <- from@yDistance
  colnames(yDistance) <- paste0("yDist.",1:ncol(yDistance))

  yRemoval <- from@yRemoval
  colnames(yRemoval) <- paste0("yRem.",1:ncol(yRemoval))

  data.frame(yDistance, yRemoval, out)
})

gdistremoval <- function(lambdaformula=~1, phiformula=~1, removalformula=~1,
  distanceformula=~1, data, keyfun=c("halfnorm", "exp", "hazard", "uniform"),
  output=c("abund", "density"), unitsOut=c("ha", "kmsq"), mixture=c('P', 'NB', 'ZIP'),
  K, starts, method = "BFGS", se = TRUE, engine=c("C","TMB"), threads=1, ...){

  keyfun <- match.arg(keyfun)
  output <- match.arg(output)
  unitsOut <- match.arg(unitsOut)
  mixture <- match.arg(mixture)
  engine <- match.arg(engine)

  formulas <- list(lambda = lambdaformula, phi = phiformula, dist = distanceformula,
                   rem = removalformula)
  if(any(sapply(formulas, has_random))) engine <- "TMB"

  M <- numSites(data)
  T <- data@numPrimary
  Rdist <- ncol(data@yDistance)
  Rrem <- ncol(data@yRemoval)
  mixture_code <- switch(mixture, P={1}, NB={2}, ZIP={3})

  gd <- getDesign(data, formulas)

  Jdist <- Rdist / T
  ysum <- array(t(gd$yDist), c(Jdist, T, M))
  ysum <- t(apply(ysum, c(2,3), sum, na.rm=T))

  Kmin = apply(ysum, 1, max, na.rm=T)

  # Parameters-----------------------------------------------------------------
  n_param <- c(ncol(gd$X_lambda), ifelse(mixture=="P",0,1),
              ifelse(T>1,ncol(gd$X_phi),0),
              ifelse(keyfun=="uniform", 0, ncol(gd$X_dist)),
              ifelse(keyfun=="hazard",1,0),
              ncol(gd$X_rem))
  nP <- sum(n_param)

  pnames <- colnames(gd$X_lambda)
  if(mixture!="P") pnames <- c(pnames, "alpha")
  if(data@numPrimary > 1) pnames <- c(pnames, colnames(gd$X_phi))
  if(keyfun!="uniform") pnames <- c(pnames, colnames(gd$X_dist))
  if(keyfun=="hazard") pnames <- c(pnames, "scale")
  pnames <- c(pnames, colnames(gd$X_rem))

  lam_ind <- 1:n_param[1]
  a_ind <- n_param[1]+1
  phi_ind <- (sum(n_param[1:2])+1):(sum(n_param[1:3]))
  dist_ind <- (sum(n_param[1:3])+1):(sum(n_param[1:4]))
  sc_ind <- (sum(n_param[1:4])+1)
  rem_ind <- (sum(n_param[1:5])+1):(sum(n_param[1:6]))

  # Distance info--------------------------------------------------------------
  db <- data@dist.breaks
  w <- diff(db)
  umf_new <- data
  umf_new@y <- umf_new@yDistance
  ua <- getUA(umf_new) #in utils.R
  u <- ua$u; a <- ua$a
  A <- rowSums(a)
  switch(data@unitsIn, m = A <- A / 1e6, km = A <- A)
  switch(unitsOut,ha = A <- A * 100, kmsq = A <- A)
  if(output=='abund') A <- rep(1, numSites(data))

  # Removal info---------------------------------------------------------------
  pl <- data@period.lengths

  # Get K----------------------------------------------------------------------
  if(missing(K) || is.null(K)) K <- max(Kmin, na.rm=TRUE) + 40

  # Using C++ engine-----------------------------------------------------------
  if(engine == "C"){

    if(missing(starts)){
      starts <- rep(0, nP)
      starts[sum(n_param[1:3])+1] <- log(median(db))
    } else if(length(starts)!=nP){
      stop(paste0("starts must be length ",sum(n_param)), call.=FALSE)
    }

    nll <- function(param){
      nll_gdistremoval(param, n_param, gd$yDist, gd$yRem, ysum, mixture_code, keyfun,
                      gd$X_lambda, A, gd$X_phi, gd$X_rem, gd$X_dist, db, a, t(u), w, pl,
                      K, Kmin, threads=threads)
    }

    opt <- optim(starts, nll, method=method, hessian=se, ...)

    covMat <- invertHessian(opt, nP, se)
    ests <- opt$par
    names(ests) <- pnames
    fmAIC <- 2 * opt$value + 2 * nP
    tmb_mod <- NULL

    #Organize fixed-effect estimates
    lambda_coef <- list(ests=ests[lam_ind], cov=as.matrix(covMat[lam_ind,lam_ind]))
    if(mixture != "P"){
      alpha_coef <- list(ests=ests[a_ind], cov=as.matrix(covMat[a_ind,a_ind]))
    }
    if(T > 1){
      phi_coef <- list(ests=ests[phi_ind], cov=as.matrix(covMat[phi_ind,phi_ind]))
    }
    if(keyfun != "uniform"){
      dist_coef <- list(ests=ests[dist_ind], cov=as.matrix(covMat[dist_ind,dist_ind]))
    }
    if(keyfun == "hazard"){
      scale_coef <- list(ests=ests[sc_ind], cov=as.matrix(covMat[sc_ind,sc_ind]))
    }
    rem_coef <- list(ests=ests[rem_ind], cov=as.matrix(covMat[rem_ind,rem_ind]))

    # No random effects in C engine
    lambda_rand_info <- phi_rand_info <- dist_rand_info <- rem_rand_info <- list()

  }

  # Using TMB engine-----------------------------------------------------------

  if(engine == "TMB"){
    if(missing(starts)) starts <- NULL

    dlist <- list(lambda=siteCovs(data), phi=yearlySiteCovs(data),
                  dist=siteCovs(data), rem=obsCovs(data))
    inps <- get_ranef_inputs(formulas, dlist,
                             gd[c("X_lambda","X_phi","X_dist","X_rem")],
                             gd[c("Z_lambda","Z_phi","Z_dist","Z_rem")])

    keyfun_type <- switch(keyfun, uniform={0}, halfnorm={1}, exp={2},
                          hazard={3})
    tmb_dat <- c(list(y_dist=gd$yDist, y_rem=gd$yRem, y_sum=ysum,
                      mixture=mixture_code, K=K, Kmin=Kmin, T=T, keyfun_type=keyfun_type,
                      A=A, db=db, a=a, w=w, u=u, per_len=pl), inps$data)

    tmb_param <- c(inps$pars, list(beta_alpha=rep(0,0), beta_scale=rep(0,0)))

    if(is.null(starts)){
      if(keyfun != "uniform") tmb_param$beta_dist[1] <- log(median(db))
    }

    if(keyfun == "hazard") tmb_param$beta_scale <- rep(0,1)
    if(keyfun == "uniform") tmb_param$beta_dist <- rep(0,0)
    if(mixture != "P") tmb_param$beta_alpha <- rep(0,1)
    if(T == 1){
      tmb_param$beta_phi <- tmb_param$b_phi <- tmb_param$lsigma_phi <- rep(0,0)
    }
    # Fit model in TMB
    tmb_out <- fit_TMB("tmb_gdistremoval", tmb_dat, tmb_param, inps$rand_ef,
                       starts=starts, method)
    tmb_mod <- tmb_out$TMB
    opt <- tmb_out$opt
    fmAIC <- tmb_out$AIC
    nll <- tmb_mod$fn

    # Organize estimates from TMB output
    lambda_coef <- get_coef_info(tmb_out$sdr, "lambda", pnames[lam_ind], lam_ind)
    lambda_rand_info <- get_randvar_info(tmb_out$sdr, "lambda", lambdaformula, siteCovs(data))
    if(mixture != "P"){
      alpha_coef <- get_coef_info(tmb_out$sdr, "alpha", pnames[a_ind], a_ind)
    }
    if(T > 1){
      phi_coef <- get_coef_info(tmb_out$sdr, "phi", pnames[phi_ind], phi_ind)
      phi_rand_info <- get_randvar_info(tmb_out$sdr, "phi", phiformula, yearlySiteCovs(data))
    }
    if(keyfun != "uniform"){
      dist_coef <- get_coef_info(tmb_out$sdr, "dist", pnames[dist_ind], dist_ind)
      dist_rand_info <- get_randvar_info(tmb_out$sdr, "dist", distanceformula, siteCovs(data))
    }
    if(keyfun == "hazard"){
      scale_coef <- get_coef_info(tmb_out$sdr, "scale", pnames[sc_ind], sc_ind)
    }
    rem_coef <- get_coef_info(tmb_out$sdr, "rem", pnames[rem_ind], rem_ind)
    rem_rand_info <- get_randvar_info(tmb_out$sdr, "rem", removalformula, obsCovs(data))

  }

  lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = lambda_coef$ests, covMat = lambda_coef$cov, fixed = 1:length(lam_ind),
    invlink = "exp", invlinkGrad = "exp", randomVarInfo=lambda_rand_info)

  estimateList <- unmarkedEstimateList(list(lambda=lamEstimates))

  if(mixture!="P"){
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = alpha_coef$ests, covMat = alpha_coef$cov,
        fixed = 1, invlink = "exp", invlinkGrad = "exp", randomVarInfo=list())
  }

  if(T>1){
    estimateList@estimates$phi <- unmarkedEstimate(name = "Availability",
        short.name = "phi", estimates = phi_coef$ests, covMat = phi_coef$cov,
        fixed = 1:length(phi_ind), invlink = "logistic", invlinkGrad = "logistic.grad",
        randomVarInfo=phi_rand_info)
  }

  if(keyfun!="uniform"){
    estimateList@estimates$dist <- unmarkedEstimate(name = "Distance",
        short.name = "dist", estimates = dist_coef$ests, covMat = dist_coef$cov,
        fixed = 1:length(dist_ind), invlink = "exp", invlinkGrad = "exp",
        randomVarInfo=dist_rand_info)
  }

  if(keyfun=="hazard"){
    estimateList@estimates$scale <- unmarkedEstimate(name = "Hazard-rate (scale)",
        short.name = "scale", estimates = scale_coef$ests, covMat = scale_coef$cov,
        fixed = 1, invlink = "exp", invlinkGrad = "exp", randomVarInfo=list())
  }

  estimateList@estimates$rem <- unmarkedEstimate(name = "Removal",
      short.name = "rem", estimates = rem_coef$ests, covMat = rem_coef$cov,
      fixed = 1:length(rem_ind), invlink = "logistic", invlinkGrad = "logistic.grad",
      randomVarInfo=rem_rand_info)

  new("unmarkedFitGDR", fitType = "gdistremoval",
    call = match.call(), formlist = formulas, data = data, 
    estimates = estimateList, sitesRemoved = numeric(0),
    AIC = fmAIC, opt = opt, negLogLike = opt$value, nllFun = nll,
    mixture=mixture, K=K, keyfun=keyfun, unitsOut=unitsOut, output=output, TMB=tmb_mod)

}

# Methods

setMethod("getP_internal", "unmarkedFitGDR", function(object){

  M <- numSites(object@data)
  T <- object@data@numPrimary
  Jrem <- ncol(object@data@yRemoval)/T
  Jdist <- ncol(object@data@yDistance)/T

  rem <- predict(object, "rem", level=NULL)$Predicted
  rem <- array(rem, c(Jrem, T, M))
  rem <- aperm(rem, c(3,1,2))

  pif <- array(NA, dim(rem))
  int_times <- object@data@period.lengths
  removalPiFun2 <- makeRemPiFun(int_times)
  for (t in 1:T){
    pif[,,t] <- removalPiFun2(rem[,,t])
  }

  phi <- rep(1, M*T)
  if(T>1) phi <- predict(object, "phi", level=NULL)$Predicted
  phi <- matrix(phi, M, T, byrow=TRUE)

  keyfun <- object@keyfun
  sig <- predict(object, "dist", level=NULL)$Predicted
  sig <- matrix(sig, M, T, byrow=TRUE)
  if(keyfun=="hazard") scale <- exp(coef(object, type="scale"))

  db <- object@data@dist.breaks
  a <- u <- rep(NA, Jdist)
  a[1] <- pi*db[2]^2
  for (j in 2:Jdist){
    a[j] <- pi*db[j+1]^2 - sum(a[1:(j-1)])
  }
  u <- a/sum(a)

  cp <- array(NA, c(M, Jdist, T))
  kf <- switch(keyfun, halfnorm=grhn, exp=grexp, hazard=grhaz,
               uniform=NULL)

  for (m in 1:M){
    for (t in 1:T){
      if(object@keyfun == "uniform"){
        cp[m,,t] <- u
      } else {
        for (j in 1:Jdist){
          cl <- call("integrate", f=kf, lower=db[j], upper=db[j+1], sigma=sig[m])
          names(cl)[5] <- switch(keyfun, halfnorm="sigma", exp="rate",
                                 hazard="shape")
          if(keyfun=="hazard") cl$scale=scale
          cp[m,j,t] <- eval(cl)$value * 2*pi / a[j] * u[j]
        }
      }
    }
  }

  #p_rem <- apply(pif, c(1,3), sum)
  #p_dist <- apply(cp, c(1,3), sum)

  out <- list(dist=cp, rem=pif)
  if(T > 1) out$phi <- phi
  out
})

# ranef

setMethod("ranef_internal", "unmarkedFitGDR", function(object, ...){

  M <- numSites(object@data)
  T <- object@data@numPrimary
  K <- object@K
  mixture <- object@mixture

  Jdist <- ncol(object@data@yDistance) / T
  ysum <- array(t(object@data@yDistance), c(Jdist, T, M))
  dist_has_na <- t(apply(ysum, c(2,3), function(x) any(is.na(x))))
  ysum <- t(apply(ysum, c(2,3), sum, na.rm=T))

  Jrem <- ncol(object@data@yRemoval) / T
  ysum_rem <- array(t(object@data@yRemoval), c(Jrem, T, M))
  rem_has_na <- t(apply(ysum_rem, c(2,3), function(x) any(is.na(x))))
  has_na <- dist_has_na | rem_has_na

  Kmin = apply(ysum, 1, max, na.rm=T)

  #loglam <- log(predict(object, "lambda", level=NULL)$Predicted)
  #lam <- exp(loglam)
  lam <- predict(object, "lambda", level=NULL)$Predicted
  if(object@output == "density"){
    ua <- getUA(object@data)
    A <- rowSums(ua$a)
    switch(object@data@unitsIn, m = A <- A / 1e6, km = A <- A)
    switch(object@unitsOut,ha = A <- A * 100, kmsq = A <- A)
    lam <- lam * A
  }

  if(object@mixture != "P"){
    alpha <- backTransform(object, "alpha")@estimate
  }

  dets <- getP(object)
  phi <- matrix(1, M, T)
  if(T > 1){
    phi <- dets$phi
  }
  cp <- dets$dist
  pif <- dets$rem

  pr <- apply(cp, c(1,3), sum)
  prRem <- apply(pif, c(1,3), sum)

  post <- array(0, c(M, K+1, 1))
  colnames(post) <- 0:K
  for (i in 1:M){
    if(mixture=="P"){
      f <- dpois(0:K, lam[i])
    } else if(mixture=="NB"){
      f <- dnbinom(0:K, mu=lam[i], size=alpha)
    } else if(mixture=="ZIP"){
      f <- dzip(0:K, lam[i], alpha)
    }

    # All sampling periods at site i have at least one missing value
    if(all(has_na[i,])){
      g <- rep(NA,K+1)
      next
    } else {
      # At least one sampling period wasn't missing
      g <- rep(1, K+1)
      for (t in 1:T){
        if(has_na[i,t]){
          next
        }
        for (k in 1:(K+1)){
          g[k] <- g[k] * dbinom(ysum[i,t], k-1, prob=pr[i,t]*prRem[i,t]*phi[i,t],
                                log=FALSE)
        }
      }
    }
    fg <- f*g
    post[i,,1] <- fg/sum(fg)
  }

  new("unmarkedRanef", post=post)
})

setMethod("simulate_internal", "unmarkedFitGDR", function(object, nsim){

  # Adjust log lambda when there is a random intercept
  #loglam <- log(predict(object, "lambda", level=NULL)$Predicted)
  #lam <- exp(loglam)
  lam <- predict(object, "lambda", level=NULL)$Predicted
  if(object@output == "density"){
    ua <- getUA(object@data)
    A <- rowSums(ua$a)
    switch(object@data@unitsIn, m = A <- A / 1e6, km = A <- A)
    switch(object@unitsOut,ha = A <- A * 100, kmsq = A <- A)
    lam <- lam * A
  }
  dets <- getP(object)

  if(object@mixture != "P"){
    alpha <- backTransform(object, "alpha")@estimate
  }

  M <- length(lam)
  T <- object@data@numPrimary

  if(T > 1){
    phi <- dets$phi
  } else {
    phi <- matrix(1, M, T)
  }

  Jrem <- dim(dets$rem)[2]
  Jdist <- dim(dets$dist)[2]

  p_dist <- apply(dets$dist, c(1,3), sum)
  p_rem <- apply(dets$rem, c(1,3), sum)

  dist_scaled <- array(NA, dim(dets$dist))
  rem_scaled <- array(NA, dim(dets$rem))
  for (t in 1:T){
    dist_scaled[,,t] <- dets$dist[,,t] / p_dist[,t]
    rem_scaled[,,t] <- dets$rem[,,t] / p_rem[,t]
  }

  p_total <- p_dist * p_rem * phi
  stopifnot(dim(p_total) == c(M, T))

  out <- vector("list", nsim)

  for (i in 1:nsim){

    switch(object@mixture,
      P = N <- rpois(M, lam),
      NB = N <- rnbinom(M, size=alpha, mu=lam),
      ZIP = N <- rzip(M, lam, alpha)
    )

    ydist <- matrix(NA, M, T*Jdist)
    yrem <- matrix(NA, M, T*Jrem)

    for (m in 1:M){
      ysum <- suppressWarnings(rbinom(T, N[m], p_total[m,]))

      ydist_m <- yrem_m <- c()

      for (t in 1:T){
        if(is.na(ysum[t])){
          yrem_m <- c(yrem_m, rep(NA, Jrem))
          ydist_m <- c(ydist_m, rep(NA, Jdist))
        } else {
          rem_class <- sample(1:Jrem, ysum[t], replace=TRUE, prob=rem_scaled[m,,t])
          rem_class <- factor(rem_class, levels=1:Jrem)
          yrem_m <- c(yrem_m, as.numeric(table(rem_class)))
          dist_class <- sample(1:Jdist, ysum[t], replace=TRUE, prob=dist_scaled[m,,t])
          dist_class <- factor(dist_class, levels=1:Jdist)
          ydist_m <- c(ydist_m, as.numeric(table(dist_class)))
        }
      }
      stopifnot(length(ydist_m)==ncol(ydist))
      stopifnot(length(yrem_m)==ncol(yrem))

      ydist[m,] <- ydist_m
      yrem[m,] <- yrem_m
    }
    out[[i]] <- list(yRemoval=yrem, yDistance=ydist)
  }
  out
})

setMethod("get_fitting_function", "unmarkedFrameGDR", 
          function(object, model, ...){
  gdistremoval
})

setMethod("y_to_zeros", "unmarkedFrameGDR", function(object, ...){
  object@yDistance[] <- 0
  object@yRemoval[] <- 0
  object
})

setMethod("rebuild_call", "unmarkedFitGDR", function(object){           
  cl <- object@call
  cl[["data"]] <- quote(object@data)
  cl[["lambdaformula"]] <- object@formlist$lambda
  cl[["phiformula"]] <- object@formlist$phi
  cl[["removalformula"]] <- object@formlist$rem
  cl[["distanceformula"]] <- object@formlist$dist
  cl[["mixture"]] <- object@mixture
  cl[["K"]] <- object@K
  cl[["keyfun"]] <- object@keyfun
  cl[["unitsOut"]] <- object@unitsOut
  cl[["output"]] <- object@output
  cl
})


setMethod("replaceY", "unmarkedFrameGDR",
          function(object, newY, replNA=TRUE, ...){

      ydist <- newY$yDistance
      stopifnot(dim(ydist)==dim(object@yDistance))
      yrem <- newY$yRemoval
      stopifnot(dim(yrem)==dim(object@yRemoval))

      if(replNA){
        ydist[is.na(object@yDistance)] <- NA
        yrem[is.na(object@yRemoval)] <- NA
      }

      object@yDistance <- ydist
      object@yRemoval <- yrem
      object
})


setMethod("SSE", "unmarkedFitGDR", function(fit, ...){
    r <- sapply(residuals(fit), function(x) sum(x^2, na.rm=T))
    return(c(SSE = sum(r)))
})


setMethod("residual_plot", "unmarkedFitGDR", function(x, ...)
{
    r <- residuals(x)
    e <- fitted(x)

    old_mfrow <- graphics::par("mfrow")
    on.exit(graphics::par(mfrow=old_mfrow))
    graphics::par(mfrow=c(2,1))

    plot(e[[1]], r[[1]], ylab="Residuals", xlab="Predicted values",
         main="Distance")
    abline(h = 0, lty = 3, col = "gray")

    plot(e[[2]], r[[2]], ylab="Residuals", xlab="Predicted values",
         main="Removal")
    abline(h = 0, lty = 3, col = "gray")
})

# Used with fitList
setMethod("fl_getY", "unmarkedFitGDR", function(fit, ...){
  getDesign(getData(fit), fit@formlist)$yDist
})
