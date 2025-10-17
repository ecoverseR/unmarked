setMethod("simulate", "unmarkedFrame",
  function(object, nsim = 1, seed = NULL, model = NULL, coefs = NULL,
           quiet = FALSE, ...){
  object <- y_to_zeros(object)
  fit <- get_fit(object, model, ...)
  coefs <- check_coefs(coefs, fit, quiet = quiet)
  lapply(1:nsim, function(i){
    coefs <- generate_random_effects(coefs, fit)
    fit <- replace_estimates(fit, coefs)
    # Cannot use nsim > 1 below as the sims would use the same set of random effects
    sim <- simulate(fit, nsim = 1)[[1]] 
    replaceY(object, sim)
  })
})

# Other fit-specific methods at the bottom of the file
setMethod("y_to_zeros", "unmarkedFrame", function(object, ...){
  object@y[] <- 0
  object
})

get_fit <- function(object, model, ...){
  fun <- get_fitting_function(object, model)
  fun(..., data = object, method = "SANN",
      control=list(maxit=0), se=FALSE)
}

# Get fitting function associated with an unmarkedFrame type
# Other fit-specific methods at the bottom of the file
setMethod("get_fitting_function", "unmarkedFrameOccu",
          function(object, model, ...){
 
  if(!(identical(model, occuRN) | identical(model, occu) | identical(model, occuPEN))){
    stop("model argument must be occu, occuRN, or occuPEN", call.=FALSE)
  }
  model
})

check_coefs <- function(coefs, fit, name = "coefs", quiet = FALSE){
  required_subs <- names(fit@estimates@estimates)
  required_coefs <- lapply(fit@estimates@estimates, function(x) names(x@estimates))

  formulas <- sapply(names(fit), function(x) get_formula(fit, x))

  # If there are random effects, adjust the expected coefficient names
  # to remove the b vector and add the grouping covariate name
  rand <- lapply(formulas, reformulas::findbars)
  if(!all(sapply(rand, is.null))){
    stopifnot(all(required_subs %in% names(formulas)))
    rvar <- lapply(rand, function(x) unlist(lapply(x, all.vars)))
    if(!all(sapply(rvar, length)<2)){
      stop("Only 1 random effect per parameter is supported", call.=FALSE)
    }
    for (i in required_subs){
      if(!is.null(rand[[i]][[1]])){
        signame <- rvar[[i]]
        old_coefs <- required_coefs[[i]]
        new_coefs <- old_coefs[!grepl("b_", old_coefs, fixed=TRUE)]
        new_coefs <- c(new_coefs, signame)
        required_coefs[[i]] <- new_coefs
      }
    }
  }
  required_lens <- lapply(required_coefs, length)

  dummy_coefs <- lapply(required_coefs, function(x){
                    out <- rep(0, length(x))
                    names(out) <- x
                    out
                  })

  if(is.null(coefs)){
    cat(name, "should be a named list of vectors, with the following structure
        (replace 0s with your values):\n\n")
    print(dummy_coefs)
    stop(paste("Specify", name, "argument as shown above"), call.=FALSE)
  }

  for (i in 1:length(required_subs)){
    if(!required_subs[i] %in% names(coefs)){
      stop(paste0("Missing required list element '",
                  required_subs[i], "' in ", name, " list"), call.=FALSE)
    }

    sub_coefs <- coefs[[required_subs[i]]]

    if(!quiet){
      message(paste0("Assumed parameter order for ", required_subs[i], ":\n",
                  paste(required_coefs[[i]], collapse=", ")))
    }
    
    if(length(sub_coefs) != required_lens[i]){
      stop(paste0("Entry '",required_subs[[i]], "' in ", name, " list must be length ",
                  required_lens[[i]]), call.=FALSE)
    }
    names(coefs[[required_subs[i]]]) <- required_coefs[[i]]

  }
  coefs[required_subs]
}

generate_random_effects <- function(coefs, fit){
  required_subs <- names(fit@estimates@estimates)
  formulas <- sapply(names(fit), function(x) get_formula(fit, x))
  rand <- lapply(formulas, reformulas::findbars)
  if(!all(sapply(rand, is.null))){
    rvar <- lapply(rand, function(x) unlist(lapply(x, all.vars)))
    for (i in required_subs){
      if(!is.null(rand[[i]][[1]])){
        signame <- rvar[[i]]
        old_coefs <- coefs[[i]]
        new_coefs <- old_coefs[names(old_coefs)!=signame]

        # Find levels of factor variable
        if(signame %in% names(siteCovs(fit@data))){
          lvldata <- siteCovs(fit@data)[[signame]]
        } else if(signame %in% names(obsCovs(fit@data))){
          lvldata <- obsCovs(fit@data)[[signame]]
        } else if(methods::.hasSlot(fit@data, "yearlySiteCovs") && signame %in% names(yearlySiteCovs(fit@data))){
          lvldata <- yearlySiteCovs(fit@data)[[signame]]
        } else {
          stop("Random effect covariate missing from data", call.=FALSE)
        }

        if(!is.factor(lvldata)){
          stop("Random effect covariates must be specified as factors", call.=FALSE)
        }
        sigma <- old_coefs[signame]
        if(sigma <= 0){
          stop("estimate for random effect represents sigma and must be positive",
               call.=FALSE)
        }
        b <- stats::rnorm(length(levels(lvldata)), 0, sigma)
        names(b) <- rep(paste0("b_",i), length(b))
        new_coefs <- c(new_coefs, b)
        coefs[[i]] <- new_coefs
      }
    }
  }
  coefs
}

replace_estimates <- function(object, new_ests){
  for (i in 1:length(new_ests)){
    est <- object@estimates@estimates[[names(new_ests)[i]]]@estimates
    stopifnot(length(est) == length(new_ests[[i]]))
    object@estimates@estimates[[names(new_ests)[i]]]@estimates <- new_ests[[i]]
  }
  object
}


# y_to_zeros-------------------------------------------------------------------

setMethod("y_to_zeros", "unmarkedFrameOccuMulti", function(object, ...){
  newy <- lapply(object@ylist, function(x){
    x[] <- 0
    x
  })
  object@ylist <- newy
  object
})


# get_fitting_function---------------------------------------------------------

setMethod("get_fitting_function", "unmarkedFrameDS", 
          function(object, model, ...){
  if(!missing(model) && identical(model, IDS)) stop("IDS not supported", call.=FALSE)
  distsamp
})

setMethod("get_fitting_function", "unmarkedFrameDSO", 
          function(object, model, ...){
  distsampOpen
})

setMethod("get_fitting_function", "unmarkedFrameGDS",
          function(object, model, ...){
  gdistsamp
})

setMethod("get_fitting_function", "unmarkedFrameGMM", 
          function(object, model, ...){
  gmultmix
})

setMethod("get_fitting_function", "unmarkedFrameGPC", 
          function(object, model, ...){
  gpcount
})

setMethod("get_fitting_function", "unmarkedFrameOccuFP",
          function(object, model, ...){
  occuFP
})

setMethod("get_fitting_function", "unmarkedFrameOccuMS",
          function(object, model, ...){
  occuMS
})

setMethod("get_fitting_function", "unmarkedFrameOccuTTD",
          function(object, model, ...){
  if(!(identical(model, occuTTD) | identical(model, nmixTTD))){
    stop("model argument must be occuTTD or nmixTTD", call.=FALSE)
  }
  model
})

setMethod("get_fitting_function", "unmarkedFrameOccuMulti",
          function(object, model, ...){
  occuMulti
})

setMethod("get_fitting_function", "unmarkedFrameMPois",
          function(object, model, ...){
  multinomPois
})

setMethod("get_fitting_function", "unmarkedFrameMMO", 
          function(object, model, ...){
  multmixOpen
})

setMethod("get_fitting_function", "unmarkedFramePCount",
          function(object, model, ...){
  pcount
})

setMethod("get_fitting_function", "unmarkedFramePCO", 
          function(object, model, ...){
  pcountOpen
})

setMethod("get_fitting_function", "unmarkedMultFrame", 
          function(object, model, ...){
  colext
})


# unmarkedFit simulate method--------------------------------------------------

# Overall simulate method for unmarkedFit objects
setMethod("simulate", "unmarkedFit", function(object, nsim = 1, seed = NULL, ...){
  simulate_internal(object, nsim = nsim)
})

# Internal methods
setMethod("simulate_internal", "unmarkedFitColExt",
  function(object, nsim){
    data <- object@data
    y <- data@y
    M <- nrow(y)	# M <- nrow(X.it)
    nY <- data@numPrimary
    J <- obsNum(data)/nY

    psiP <- predict(object, type = "psi", level = NULL, na.rm = FALSE)$Predicted
    detP <- predict(object, type = "det", level = NULL, na.rm = FALSE)$Predicted
    colP <- predict(object, type = "col", level = NULL, na.rm = FALSE)$Predicted
    extP <- predict(object, type = "ext", level = NULL, na.rm = FALSE)$Predicted

    detP <- array(detP, c(J, nY, M))
    detP <- aperm(detP, c(3, 1, 2))
    colP <- matrix(colP, M, nY, byrow = TRUE)
    extP <- matrix(extP, M, nY, byrow = TRUE)

    simList <- list()
    for(s in 1:nsim) {
        ## generate first year's data
        x <- matrix(0, M, nY)
        x[,1] <- stats::rbinom(M, 1, psiP)

        ## create transition matrices (phi^T)
        phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
        for(i in 1:M) {
            for(t in 1:(nY-1)) {
                phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t],
                    1-extP[i,t]))
                }
            }

        ## generate latent years 2:T
        for(i in 1:M) {
            for(t in 2:nY) {
                x[i,t] <- stats::rbinom(1, 1, phis[2,x[i,t-1]+1,t-1,i])
                }
            }

        ## generate observations
        y <- array(NA, c(M, J, nY))
        for(t in 1:nY) {
            y[,,t] <- stats::rbinom(M*J, 1, x[,t]*detP[,,t])
            }

        y.mat <- y[,,1]
        for(i in 2:dim(y)[3]) {
            y.mat <- cbind(y.mat,y[,,i])
            }
        simList[[s]] <- y.mat
        }
    return(simList)
})


setMethod("simulate_internal", "unmarkedFitDS",
  function(object, nsim){ 
  lambda <- predict(object, type = "state", level = NULL, na.rm=FALSE)$Predicted
  if(object@output == "density"){
    lambda <- lambda * get_ds_area(object@data, object@unitsOut)
  }

  p <- getP(object, na.rm = FALSE)
  M <- nrow(p)
  J <- ncol(p)

  lapply(1:nsim, function(i){
    yvec <- stats::rpois(M*J, lambda * p)
    matrix(yvec, M, J)
  })
})


setMethod("simulate_internal", "unmarkedFitGDS",
    function(object, nsim){

    umf <- object@data
    db <- umf@dist.breaks
    w <- diff(db)
    mixture <- object@mixture
    keyfun <- object@keyfun

    M <- nrow(umf@y)
    T <- umf@numPrimary
    R <- ncol(umf@y)
    J <- R / T

    if(T>1){
      phi <- predict(object, type = "phi", level = NULL, na.rm=FALSE)$Predicted
    } else {
      phi <- rep(1, M*T)
    }
    phi <- matrix(phi, M, T, byrow=TRUE)
    
    detPars <- coef(object, type="det")

    lambda <- predict(object, type = "lambda", level = NULL, na.rm=FALSE)$Predicted
    if(object@output == "density"){
      lambda <- lambda * get_ds_area(object@data, object@unitsOut)
    }

    cp <- getP(object, na.rm = FALSE)
    ysim <- cpa <- array(cp, c(M, J, T))

    simList <- list()
    for(s in 1:nsim) {
        for(i in 1:M) {
            switch(mixture,
                P = Ns <- rpois(1, lambda[i]),
                NB = {
                    alpha <- exp(coef(object, type="alpha"))
                    Ns <- rnbinom(1, mu=lambda[i], size=alpha)
                    },
                ZIP = {
                    psi <- plogis(coef(object['psi']))
                    Ns <- rzip(1, lambda[i], psi)
                    }
            )
            for(t in 1:T) {
                N <- rbinom(1, Ns, phi[i,t])
                cp.it <- cpa[i,,t]
                cp.it[J+1] <- 1-sum(cp.it)
                y.it <- as.integer(rmultinom2(1, N, prob=cp.it))
                ysim[i,,t] <- y.it[1:J]
                }
            }
        y.mat <- matrix(ysim, nrow=M)
        simList[[s]] <- y.mat
        }
    return(simList)
})


setMethod("simulate_internal", "unmarkedFitGMM",
  function(object, nsim){
    umf <- object@data
    mixture <- object@mixture
    y <- umf@y
    n <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T

    lam <- predict(object, type = "lambda", level = NULL, na.rm=FALSE)$Predicted
    
    if(T == 1){
      phi <- rep(1, n*T)
    } else {
      phi <- predict(object, type = "phi", level = NULL, na.rm=FALSE)$Predicted
    }
    phi.mat <- matrix(phi, nrow=n, ncol=T, byrow=TRUE)

    cp.arr <- array(NA, c(n, T, J+1))
    cp.mat <- getP(object, na.rm = FALSE)
    cp.mat[is.na(y)] <- NA
    cp.temp <- array(cp.mat, c(n, J, T))
    cp.arr[,,1:J] <- aperm(cp.temp, c(1,3,2))
    cp.arr[,, 1:J][is.na(y)]<- NA   # Andy added 5/30
    cp.arr[,,J+1] <- 1 - apply(cp.arr[,,1:J,drop=FALSE], 1:2, sum, na.rm=TRUE)

    simList <- list()
    for(s in 1:nsim) {
        switch(mixture,
            P = M <- rpois(n=n, lambda=lam),
            NB = M <- rnbinom(n=n, mu=lam,
                size=exp(coef(object, type="alpha"))),
            ZIP = {
              psi <- plogis(coef(object['psi']))
              M <- rzip(n, lambda=lam, psi=psi)
            }
        )

        N <- rbinom(n*T, size=M, prob=phi.mat)
        # bug fix 3/16/2010
        N <- matrix(N, nrow=n, ncol=T, byrow=FALSE) # , byrow=TRUE)

        y.sim <- array(NA, c(n, J, T))
        for(i in 1:n) {
            for(t in 1:T) {
                if(is.na(N[i,t]))
                    next
                pi.it <- cp.arr[i,t,]
                na.it <- is.na(pi.it)
                pi.it[na.it] <- 0
                y.sim[i,,t] <- drop(rmultinom2(1, N[i,t], pi.it))[1:J]
                y.sim[i,na.it[1:J],t] <- NA
            }
        }
        simList[[s]] <- matrix(y.sim, nrow=n, ncol=J*T) # note, byrow=F
        }
    return(simList)
})


setMethod("simulate_internal", "unmarkedFitGPC",
    function(object, nsim){
    umf <- object@data
    mixture <- object@mixture
    y <- umf@y
    R <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T

    lam <- predict(object, type = "lambda", level = NULL, na.rm=FALSE)$Predicted

    if(T == 1){
      phi <- rep(1, M * T)
    } else {
      phi <- predict(object, type = "phi", level = NULL, na.rm = FALSE)$Predicted
    }
    phi.mat <- matrix(phi, nrow=R, ncol=T, byrow=TRUE)

    p <- getP(object, na.rm=FALSE)
    p <- array(p, c(R, J, T))

    simList <- list()
    for(s in 1:nsim) {
        switch(mixture,
            P = M <- rpois(n=R, lambda=lam),
            NB = M <- rnbinom(n=R, mu=lam,
               size=exp(coef(object, type="alpha"))),   
            ZIP = {
              psi <- plogis(coef(object['psi']))
              M <- rzip(R, lambda=lam, psi=psi)
            }
        )

        N <- rbinom(R*T, size=M, prob=phi.mat)
        N <- matrix(N, nrow=R, ncol=T, byrow=FALSE)

        y.sim <- array(NA, c(R, J, T))
        for(i in 1:R) {
            for(t in 1:T) {
                if(is.na(N[i,t]))
                    next
                y.sim[i,,t] <- rbinom(J, N[i,t], p[i,,t])
            }
        }
        y.sim <- matrix(y.sim, nrow=R, ncol=J*T)
        #y.sim[is.na(y)] <- NA # Not necessary if covariates exist!
        simList[[s]] <- y.sim
    }
    return(simList)
})


setMethod("simulate_internal", "unmarkedFitMPois",
  function(object, nsim){
  p <- getP(object, na.rm = FALSE)
  M <- nrow(p)
  J <- ncol(p)
  p <- as.vector(t(p))
  lam <- predict(object, type = "state", level = NULL, na.rm=FALSE)$Predicted  
  lam <- rep(lam, each = J)
  lam_p <- lam * p

  lapply(1:nsim, function(i){
    yvec <- stats::rpois(M*J, lam_p)
    matrix(yvec, M, J, byrow = TRUE)
  })
})


setMethod("simulate_internal", "unmarkedFitNmixTTD",
  function(object,  nsim){

  M <- nrow(object@data@y)
  J <- ncol(object@data@y)
  tdist <- ifelse("shape" %in% names(object@estimates), "weibull", "exp")
  mix <- ifelse("alpha" %in% names(object@estimates), "NB", "P")

  #Get predicted values
  abun <- predict(object, 'state', level = NULL, na.rm=FALSE)$Predicted
  lam <- predict(object, 'det', level = NULL, na.rm=FALSE)$Predicted
  tmax <- as.vector(t(object@data@surveyLength))
  not_na <- !is.na(lam)

  simlist <- list()
  for(s in 1:nsim){

    if(mix=="P"){
      N <- rpois(M, abun)
    } else if(mix=="NB"){
      alpha <- exp(coef(object, "alpha"))
      N <- rnbinom(M, mu=abun, size=alpha)
    }

    lamN <- lam*rep(N, each=J)
    is_zero <- lamN==0
    to_sim <- !is_zero & not_na

    ttd <- rep(NA, length(lamN))
    ttd[is_zero] <- tmax[is_zero]

    if(tdist == "weibull"){
      k <- exp(coef(object)['k(k)'])
      ttd[to_sim] <- stats::rweibull(sum(to_sim),k,1/lamN[to_sim])
    } else {
      ttd[to_sim] <- stats::rexp(sum(to_sim), lamN[to_sim])
    }
    #Truncate
    ttd[ttd>tmax] <- tmax[ttd>tmax]

    ttd <- matrix(ttd, nrow=M, byrow=T)

    #Add NAs
    nas <- which(is.na(object@data@y))
    ttd[nas] <- NA
    simlist[[s]] <- ttd
  }
  simlist
})


setMethod("simulate_internal", "unmarkedFitOccu",
  function(object, nsim){
  p <- getP(object, na.rm = FALSE)
  M <- nrow(p)
  J <- ncol(p)
  p <- as.vector(t(p))
  psi <- predict(object, type = "state", level = NULL, na.rm = FALSE)$Predicted
  
  lapply(1:nsim, function(i){
    z <- stats::rbinom(M, 1, psi)
    z[object@knownOcc] <- 1
    z <- rep(z, each = J)
    yvec <- z * stats::rbinom(M*J, 1, p)
    matrix(yvec, M, J, byrow = TRUE)
  })
})


setMethod("simulate_internal", "unmarkedFitOccuFP",
  function(object, nsim){
  psi <- predict(object, type = "state", level = NULL, na.rm = FALSE)$Predicted
  p <- getP(object, na.rm=FALSE)
  M <- nrow(p)
  J <- ncol(p)
  p <- as.vector(t(p))
  fp <- as.vector(t(getFP(object)))
  b <- as.vector(t(getB(object)))
            
  simList <- list()
  for(i in 1:nsim) {
    Z <- rbinom(M, 1, psi)
    Z[object@knownOcc] <- 1
    Z <- rep(Z, each = J)
    P <- matrix(0,M*J,3)
    P[,1] <- Z*rbinom(M * J, 1, prob = (1-p)) + (1-Z)*rbinom(M * J, 1, prob = (1-fp))
    P[,2] <- (1-P[,1])*(1-Z) + (1-P[,1])*rbinom(M * J, 1, prob = (1-b))*Z
    P[,3] <- 1 - P[,1]-P[,2]
    yvec <- sapply(1:(M*J),function(x) which(as.logical(rmultinom2(1,1,P[x,])))-1)
    simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
  }
  return(simList)
})


setMethod("simulate_internal", "unmarkedFitOccuMS",
  function(object, nsim){

  S <- object@data@numStates
  N <- numSites(object@data)
  T <- object@data@numPrimary
  J <- obsNum(object@data) / T

  prm <- object@parameterization
  # TODO: update syntax here
  psi_raw <- predict(object, "psi", level=NULL)
  psi_raw <- sapply(psi_raw, function(x) x$Predicted)
  p <- getP(object)

  guide <- matrix(NA,nrow=S,ncol=S)
  guide <- lower.tri(guide,diag=TRUE)
  guide[,1] <- FALSE
  guide <- which(guide,arr.ind=TRUE)

  out <- vector("list",nsim)

  for (i in 1:nsim){

  #State process
  if(prm == "multinomial"){
    psi <- cbind(1-apply(psi_raw,1,sum),psi_raw)
  } else if (prm == "condbinom"){
    psi <- matrix(NA, nrow=N, ncol=S)
    psi[,1] <- 1-psi_raw[,1]
    psi[,2] <- (1-psi_raw[,2])*psi_raw[,1]
    psi[,3] <- psi_raw[,1]*psi_raw[,2]
  }

  z <- matrix(NA, nrow=N, ncol=T)

  #initial occupancy
  for (n in 1:N){
    z[n,1] <- sample(0:(S-1), 1, prob=psi[n,])
  }

  #transitions if T>1----------------------------------------------------------
  get_phimat <- function(prob_vec){
    if(prm=="multinomial"){
      out <- matrix(NA, nrow=S, ncol=S)
      out[outer(1:S, 1:S, function(i,j) i!=j)] <- prob_vec
      out <- t(out)
      diag(out) <- 1 - rowSums(out,na.rm=T)
      return(out)
    } else if(prm == "condbinom"){
      out <- matrix(prob_vec, nrow=S)
      return(cbind(1-out[,1], out[,1]*(1-out[,2]), out[,1]*out[,2]))
    }
  }

  if(T>1){
    phi_raw <- predict(object, "phi", level=NULL)
    phi_raw <- sapply(phi_raw, function(x) x$Predicted)
    phi_index <- 1
    for (n in 1:N){
      for (t in 2:T){
        phimat <- get_phimat(phi_raw[phi_index,])
        phi_t <- phimat[(z[n,(t-1)]+1),]
        z[n,t] <- sample(0:(S-1), 1, prob=phi_t)
        phi_index <- phi_index+1
      }
    }
  }
  #----------------------------------------------------------------------------

  #Detection process
  y <- matrix(0, nrow=N, ncol=J*T)
  for (n in 1:N){
    yindex <- 1
    for (t in 1:T){
      for (j in 1:J){
        if(prm == "multinomial"){
          probs_raw <- sapply(p, function(x) x[n,yindex])
          # Make sure output is NA if probs have NA
          if(any(is.na(probs_raw))){
            y[n,yindex] <- NA
            yindex <- yindex + 1
            next
          }

          sdp <- matrix(0, nrow=S, ncol=S)
          sdp[guide] <- probs_raw
          sdp[,1] <- 1 - rowSums(sdp)

          probs <- sdp[z[n,t]+1,]

        } else if (prm == "condbinom"){
          p11 <- p[[1]][n,yindex]
          p12 <- p[[2]][n,yindex]
          p22 <- p[[3]][n,yindex]
          # Trap NAs in probability of detection
          if(any(is.na(c(p11, p12, p22)))){
            y[n,yindex] <- NA
            next
          }
          probs <- switch(z[n,t]+1,
                          c(1,0,0),
                          c(1-p11,p11,0),
                          c(1-p12,p12*(1-p22),p12*p22))
        }
        # this NA trap probably isn't necessary but leaving it in just in case
        if(all(!is.na(probs))){
          y[n,yindex] <- sample(0:(S-1), 1, prob=probs)
        } else {
          y[n,yindex] <- NA
        }
        yindex <- yindex + 1
      }
    }
  }

  out[[i]] <- y
  }

  out
})


setMethod("simulate_internal", "unmarkedFitOccuMulti",
    function(object, nsim){
    data <- object@data
    ynames <- names(object@data@ylist)
    maxOrder <- object@call$maxOrder
    # TODO: put maxOrder in output object?
    if(is.null(maxOrder)) maxOrder <- length(object@data@ylist)
    dm <- getDesign(object@data, object@formlist, maxOrder)
    # TODO: standardize this
    psi <- predict(object, "state", level=NULL)$Predicted
    p <- getP(object)
    stopifnot(nrow(psi) == nrow(object@data@y))

    simList <- list()
    for (s in 1:nsim){
      #True state
      ztruth <- matrix(NA,nrow=dm$N,ncol=dm$S)
      nZ <- nrow(dm$z)
      for (i in 1:dm$N){
        ztruth[i,] <- as.matrix(dm$z[sample(nZ,1,prob=psi[i,]),])
      }

      y <- list()
      for (i in 1:dm$S){
        y[[i]] <- matrix(NA,dm$N,dm$J)
        for (j in 1:dm$N){
          for (k in 1:dm$J){
            if(!is.na(p[[i]][j,k])){
              y[[i]][j,k] <- stats::rbinom(1,1,ztruth[j,i]*p[[i]][j,k])
            }
          }
        }
      }
      names(y) <- ynames
      simList[[s]] <- y
    }
    simList
})


setMethod("simulate_internal", "unmarkedFitOccuRN",
  function(object, nsim){
  lambda <- predict(object, type = "state", level = NULL, na.rm = FALSE)$Predicted
  r <- getP(object, na.rm = FALSE) # individual-level detection prob
  M <- nrow(r)
  J <- ncol(r)
  r <- as.vector(t(r))

  lapply(1:nsim, function(i){
    N <- stats::rpois(M, lambda)
    N <- rep(N, each = J)
    p <- 1 - (1 - r)^N 
    yvec <- stats::rbinom(M*J, 1, p)
    matrix(yvec, M, J, byrow = TRUE)
  })
})


setMethod("simulate_internal", "unmarkedFitOccuTTD",
  function(object,  nsim){

  N <- nrow(object@data@y)
  T <- object@data@numPrimary
  J <- ncol(object@data@y)/T
  tdist <- ifelse("shape" %in% names(object@estimates), "weibull", "exp")

  #Get predicted values
  psi <- predict(object, type='psi', level=NULL, na.rm=FALSE)$Predicted
  lam <- predict(object, type='det', level=NULL, na.rm=FALSE)$Predicted
  if(T>1){
    p_col <- predict(object, type='col', level = NULL, na.rm=FALSE)$Predicted
    p_col <- matrix(p_col, N, T, byrow=TRUE)
    p_ext <- predict(object, type='ext', level = NULL, na.rm=FALSE)$Predicted
    p_ext <- matrix(p_ext, N, T, byrow=TRUE)
  }

  tmax <- object@data@surveyLength
  not_na <- which(!is.na(lam))

  simlist <- list()
  for(s in 1:nsim){
    ttd <- rep(NA, length(lam))
    if(tdist == "weibull"){
      k <- exp(coef(object)['k(k)'])
      ttd[not_na] <- stats::rweibull(length(not_na),k,1/lam[not_na])
    } else {
      ttd[not_na] <- stats::rexp(length(not_na), lam[not_na])
    }
    #Truncate
    ttd <- matrix(ttd, nrow=N, byrow=T)
    ttd[which(ttd>tmax)] <- tmax[which(ttd>tmax)]

    #Latent state
    z <- matrix(NA, N, T)
    z[,1] <- rbinom(N, 1, psi)

    if(T>1){
      for (t in 1:(T-1)){
        z_ext <- rbinom(N, 1, 1-p_ext[,t])
        z_col <- rbinom(N, 1, p_col[,t])
        z[,t+1] <- ifelse(z[,t], z_ext, z_col)
      }
    }

    #Detection process
    yout <- matrix(NA, N, J*T)
    d_ind <- 1
    for (t in 1:T){
      for (j in 1:J){
        yout[,d_ind] <- ifelse(z[,t], ttd[,d_ind], tmax[,d_ind])
        d_ind <- d_ind + 1
      }
    }

    #Add NAs
    nas <- which(is.na(object@data@y))
    yout[nas] <- NA
    simlist[[s]] <- yout
  }
  simlist
})


setMethod("simulate_internal", "unmarkedFitPCount",
  function(object, nsim){
  # Can't use predict because it will be incorrect with ZIP
  dm <- getDesign(object@data, object@formlist, na.rm = FALSE)
  lam <- exp(dm$X_state %*% coef(object, "state") + dm$offset_state)
  p <- getP(object, na.rm=FALSE)
  M <- nrow(p)
  J <- ncol(p)
  p <- as.vector(t(p))

  lapply(1:nsim, function(i){
    N <- switch(object@mixture,
      P = rpois(M, lam),
      NB = rnbinom(M, size = exp(coef(object["alpha"])), mu = lam),
      ZIP = rzip(M, lam, plogis(coef(object["psi"])))
    )
    N <- rep(N, each = J)

    yvec <- rbinom(M*J, size = N, prob = p)
    matrix(yvec, M, J, byrow=TRUE)
  })
})

# Open-population models-------------------------------------------------------

#Simulate open-population abundance
simOpenN <- function(object){
    mix <- object@mixture
    dynamics <- object@dynamics
    umf <- object@data
    #To partially handle old saved model objects
    fix <- tryCatch(object@fix, error=function(e) "none")
    immigration <- tryCatch(object@immigration, error=function(e) FALSE)
    dm <- getDesign(umf, object@formlist, na.rm = FALSE)
    delta <- dm$delta
    y <- umf@y
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T

    # Can't use predict because it will be incorrect with ZIP
    lambda <- exp(dm$X_lambda %*% coef(object, "lambda") + dm$offset_lambda)

    if(fix == "gamma"){
        gamma <- matrix(0, M, T-1)
    } else if(dynamics != "notrend"){
        gamma <- predict(object, type = "gamma", level = NULL, na.rm=FALSE)$Predicted
        gamma <- matrix(gamma, M, T-1, byrow=TRUE)
    } else {
        gamma <- matrix(NA, M, T-1)
    }
    if (fix == "omega"){
        omega <- matrix(1, M, T-1)
    } else if(dynamics == "trend"){
        omega <- matrix(-999, M, T-1) # placeholder
    } else { # this should cover both situations with log and logit link automatically
        omega <- predict(object, type = "omega", level = NULL, na.rm = FALSE)$Predicted
        omega <- matrix(omega, M, T-1, byrow = TRUE)
    }

    if(immigration){
        iota <- predict(object, type = "iota", level = NULL, na.rm = FALSE)$Predicted
        iota <- matrix(iota, M, T-1, byrow=TRUE)
    } else {
        iota <- matrix(0, M, T-1)
    }

    if(identical(mix, "ZIP")){
        psi <- plogis(coef(object, type="psi"))
    }

    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    for(i in 1:M) {
            switch(mix,
                   P = N[i, 1] <- stats::rpois(1, lambda[i]),
                   NB = N[i, 1] <- stats::rnbinom(1, size =
                       exp(coef(object["alpha"])), mu = lambda[i]),
                   ZIP = N[i,1] <- rzip(1, lambda[i], psi))
            if(delta[i, 1] > 1) {
                for(d in 2:delta[i, 1]) {
                    if(dynamics == "trend")
                        N[i,1] <- stats::rpois(1, N[i,1]*gamma[i,1]+iota[i,1])
                    else if(dynamics == "ricker")
                        N[i,1] <- stats::rpois(1, N[i,1]*exp(gamma[i, 1] * (1 -
                          N[i, 1] / omega[i, 1])) + iota[i, 1])
                    else if(dynamics == "gompertz")
                        N[i,1] <- stats::rpois(1, N[i, 1] * exp(gamma[i, 1] * (1 -
                          log(N[i, 1] + 1)/log(omega[i, 1] + 1))) +
                          iota[i, 1])
                    else if(dynamics == "autoreg")
                        N[i,1] <- stats::rbinom(1, N[i,1], omega[i,1]) +
                            stats::rpois(1, gamma[i,1]*N[i,1] + iota[i, 1])
                    else
                        N[i,1] <- stats::rbinom(1, N[i,1], omega[i,1]) +
                            stats::rpois(1, gamma[i, 1])
                }
            }
            # Might need more NA handling here...ignore gaps?
            for(t in 2:T) {
                if(is.na(omega[i, t-1]) | is.na(gamma[i,t-1]))
                    N[i, t] <- N[i, t-1] # just a place holder
                else{
                    if(!identical(dynamics, "trend"))
                        S[i, t-1] <- stats::rbinom(1, N[i, t-1], omega[i, t-1])
                    if(identical(dynamics, "autoreg"))
                        gamma[i, t-1] <- gamma[i, t-1] * N[i,t-1] + iota[i,t-1]
                    else if(identical(dynamics, "notrend"))
                        gamma[i, t-1] <- (1-omega[i, t-1]) * lambda[i]
                    G[i, t-1] <- stats::rpois(1, gamma[i, t-1])
                    if(identical(dynamics, "trend"))
                        N[i,t] <- stats::rpois(1, N[i,t-1]*gamma[i,t-1]+iota[i,t-1])
                    else if(identical(dynamics, "ricker"))
                        N[i,t] <- stats::rpois(1, N[i,t-1]*exp(gamma[i,t-1]*(1-N[i,
                            t-1]/omega[i,t-1]))+iota[i,t-1])
                    else if(identical(dynamics, "gompertz"))
                        N[i,t] <- stats::rpois(1, N[i,t-1]*exp(gamma[i,t-1]*(1-
                            log(N[i,t-1] + 1) / log(omega[i,t-1] + 1))) +
                            iota[i,t-1])
                    else
                        N[i, t] <- S[i, t-1] + G[i, t-1]
                    if(delta[i, t] > 1) {
                        for(d in 2:delta[i, 1]) {
                            if(dynamics == "trend")
                                N[i,t] <- stats::rpois(1, N[i,t]*gamma[i,t-1]+
                                    iota[i,t-1])
                            else if(identical(dynamics, "ricker"))
                                N[i,t] <- stats::rpois(1, N[i,t]*exp(gamma[i,t-1]*(1-
                                    N[i,t]/omega[i,t-1]))+iota[i,t-1])
                            else if(identical(dynamics, "gompertz"))
                                N[i,t] <- stats::rpois(1, N[i,t]*exp(gamma[i,t-1]*(1-
                                    log(N[i,t] + 1)/
                                    log(omega[i,t-1] + 1))) + iota[i,t-1])
                            else {
                                S[i,t-1] <- stats::rbinom(1, N[i,t], omega[i,t-1])
                                G[i,t-1] <- stats::rpois(1, gamma[i, t-1])
                                N[i, t] <- S[i, t-1] + G[i, t-1]
                            }
                        }
                    }
                }
            }
        }
    N
}

setMethod("simulate_internal", "unmarkedFitPCO",
  function(object, nsim){

    umf <- object@data
    M <- numSites(umf)
    T <- umf@numPrimary
    J <- ncol(getY(umf)) / T
    p <- getP(object, na.rm = FALSE) # M x TJ
  
    simList <- list()
    for(s in 1:nsim) {
        y.sim <- matrix(NA, M, J*T)
        N <- simOpenN(object) # M x T
        N <- N[,rep(1:T, each=J)]
        yvec <- rbinom(M*T*J, as.vector(N), as.vector(p))
        ymat <- matrix(yvec, M, T*J)
        simList[[s]] <- ymat
    }
    return(simList)
})


# Covers DSO and MMO
setMethod("simulate_internal", "unmarkedFitDailMadsen",
  function(object, nsim){

  umf <- object@data
  M <- numSites(umf)
  T <- umf@numPrimary
  J <- ncol(getY(umf)) / T

  p <- getP(object, na.rm = FALSE)
  p <- array(p, c(M,J,T))
  cp <- array(NA, c(M,J+1,T))
  for (i in 1:M){
    for (t in 1:T){
      cp[i, 1:J, t] <- p[i,,t]
      cp[i, J+1, t] <- 1 - sum(p[i,,t], na.rm=TRUE)
    }
  }

  simList <- list()
  for(s in 1:nsim) {
    y.sim <- matrix(NA, M, J*T)
    N <- simOpenN(object)

    for(i in 1:M) {
      yst <- 1
        for(t in 1:T) {
          yend <- yst + J - 1
          #rmultinom2 in utils.R
          y.it <- as.integer(rmultinom2(1, N[i,t], prob=cp[i,,t]))
          y.sim[i,yst:yend] <- y.it[1:J]
          yst <- yst + J
        }
    }
    simList[[s]] <- y.sim
  }
  return(simList)
})
