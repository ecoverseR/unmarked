# Square bracket selector methods (exported)-----------------------------------
setMethod("[", c("unmarkedFrame", "numericOrLogical", "missing", "missing"),
  function(x, i){
  subset_sites(x, i)
})

setMethod("[", c("unmarkedFrame", "missing", "numericOrLogical", "missing"),
  function(x, i, j){
  subset_obs(x, j)
})

setMethod("[", c("unmarkedFrame", "numericOrLogical", "numericOrLogical", "missing"),
  function(x, i, j){
  x <- x[i,]
  x[,j]
})


# Internal site-selector methods-----------------------------------------------

setMethod("subset_sites", c("unmarkedFrame", "numericOrLogical"),
  function(umf, i){
 
  M <- numSites(umf)
  J <- obsNum(umf)

  # Process selection
  keep <- process_site_index(umf, i)

  # Subset y
  y <- getY(umf)[keep,,drop=FALSE]
  umf@y <- y

  # Subset site covs
  sc <- siteCovs(umf)
  if(!is.null(sc)){
    siteCovs(umf) <- sc[keep,,drop=FALSE]
  }

  # Subset obs covs
  oc <- obsCovs(umf)
  if(!is.null(oc)){
    keep_obs <- rep(1:M, each = J)
    oc_new <- lapply(keep, function(x){
      oc[keep_obs == x,,drop=FALSE]
    })
    obsCovs(umf) <- do.call(rbind, oc_new)
  }

  # Subset yearly site covs if they exist
  if(methods::.hasSlot(umf, "yearlySiteCovs")){
    ysc <- yearlySiteCovs(umf)
    keep_per <- rep(1:M, each = umf@numPrimary)
    if(!is.null(ysc)){
      ysc_new <- lapply(keep, function(x){
        ysc[keep_per == x,,drop=FALSE]
      })
      yearlySiteCovs(umf) <- do.call(rbind, ysc_new)
    }
  }

  umf
})

setMethod("subset_sites", c("unmarkedFrameDailMadsen", "numericOrLogical"),
  function(umf, i){
  idx <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  umf@primaryPeriod <- umf@primaryPeriod[idx,,drop=FALSE]
  umf
})

setMethod("subset_sites", c("unmarkedFrameDS", "numericOrLogical"),
  function(umf, i){
  idx <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  if(umf@survey == "line"){
    umf@tlength <- umf@tlength[idx] 
  }
  umf
})

setMethod("subset_sites", c("unmarkedFrameDSO", "numericOrLogical"),
  function(umf, i){
  idx <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  if(umf@survey == "line"){
    umf@tlength <- umf@tlength[idx] 
  }
  umf
})

setMethod("subset_sites", c("unmarkedFrameGDR", "numericOrLogical"),
  function(umf, i){
  keep <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  umf@yDistance <- umf@yDistance[keep,,drop=FALSE]
  umf@yRemoval <- umf@yRemoval[keep,,drop=FALSE]
  umf
})

setMethod("subset_sites", c("unmarkedFrameGDS", "numericOrLogical"),
  function(umf, i){
  idx <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  if(umf@survey == "line"){
    umf@tlength <- umf@tlength[idx] 
  }
  umf
})

# This one is more complicated due to species covs and repeated site/obs covs
setMethod("subset_sites", c("unmarkedFrameOccuComm", "numericOrLogical"),
    function(umf, i){
  if(length(i) == 0) return(umf)
  i <- process_site_index(umf, i)
  M <- numSites(umf)
  J <- obsNum(umf)
  S <- length(umf@ylist)

  ylist <- lapply(umf@ylist,function(x) x[i,,drop=F])
 
  siteCovs <- siteCovs(umf)
  if (!is.null(siteCovs)) {
    siteCovs <- siteCovs(umf)[i, , drop = FALSE]
  }

  obsCovs <- obsCovs(umf)
  if (!is.null(obsCovs)) {
    .site <- rep(1:M, each = J)
    oc <- lapply(i, function(ind){
      obsCovs[.site==ind,,drop=FALSE]
    })
    obsCovs <- do.call(rbind, oc)
  }

  # Species covs
  spc <- umf@speciesCovs
  if(!is.null(spc)){
    # length S covs are unchanged
    spc_sp <- sapply(spc, function(x) identical(length(x), S))
    spc_sp <- spc[spc_sp]

    # M x S covs
    spc_site <- sapply(spc, function(x) identical(dim(x), c(M, S))) 
    spc_site <- spc[spc_site]
    if(length(spc_site) > 0){
      spc_site <- lapply(spc_site, function(x){
        x[i,,drop=FALSE]
      })
    }
    # M x J x S covs
    spc_obs <- sapply(spc, function(x) identical(dim(x), c(M, J, S)))
    spc_obs <- spc[spc_obs]
    if(length(spc_obs) > 0){
      spc_obs <- lapply(spc_obs, function(x){
        x[i,,,drop=FALSE]
      })
    }
    new_spc <- c(spc_site, spc_obs, spc_sp)
  } else {
    new_spc <- NULL
  }

  umf@y <- ylist[[1]]
  umf@ylist <- ylist
  umf@siteCovs <- siteCovs
  umf@obsCovs <- obsCovs
  umf@speciesCovs <- new_spc
  umf
})

setMethod("subset_sites", c("unmarkedFrameOccuCOP", "numericOrLogical"),
  function(umf, i){
  keep <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  umf@L <- umf@L[keep,,drop=FALSE]
  umf
})

setMethod("subset_sites", c("unmarkedFrameOccuMulti", "numericOrLogical"),
  function(umf, i){
  keep <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  umf@ylist <- lapply(umf@ylist, function(s) s[keep,,drop=FALSE])
  umf
})

setMethod("subset_sites", c("unmarkedFrameOccuTTD", "numericOrLogical"),
  function(umf, i){
  idx <- process_site_index(umf, i)
  umf <- methods::callNextMethod(umf, i)
  umf@surveyLength <- umf@surveyLength[idx,,drop=FALSE]
  umf
})

# Utility function to clean up indices
process_site_index <- function(umf, i){
  if(is.logical(i)){
    stopifnot(length(i) == numSites(umf))
    keep <- 1:numSites(umf)
    keep <- keep[i]
  } else if(is.numeric(i)){
    stopifnot(all(abs(i) %in% 1:numSites(umf)))
    stopifnot(all(i>0) | all(i<0))
    if(all(i<0)){
      keep <- 1:numSites(umf)
      keep <- keep[!keep %in% abs(i)]
    } else {
      keep <- i
    }
  }
  keep
}


# Internal obs/period-selector methods-----------------------------------------

setMethod("subset_obs", c("unmarkedFrame", "numericOrLogical"),
  function(umf, j){
  if(obsNum(umf) != numY(umf)){
    stop("Can't subset multinomial observations", call.=FALSE)
  }
  M <- numSites(umf)
  J <- obsNum(umf)

  # Process selection
  keep <- process_obs_index(umf, j)

  # Subset y
  y <- getY(umf)[,keep,drop=FALSE]
  umf@y <- y

  # Subset obs covs
  oc <- obsCovs(umf)
  if(!is.null(oc)){
    keep_site <- rep(1:M, each = J)
    oc_new <- lapply(1:M, function(m){
      oc_sub <- oc[keep_site == m,,drop=FALSE]
      oc_sub <- oc_sub[keep,,drop=FALSE]
    })
    obsCovs(umf) <- do.call(rbind, oc_new)
  }

  # Subset obsToY
  # This only works with non-multinomial
  obsToY(umf) <- diag(length(keep))

  umf
})

setMethod("subset_obs", c("unmarkedFrameOccuCOP", "numericOrLogical"),
  function(umf, j){
  keep <- process_obs_index(umf, j)
  umf <- methods::callNextMethod(umf, j)
  umf@L <- umf@L[,keep,drop=FALSE]
  umf
})

setMethod("subset_obs", c("unmarkedFrameOccuMulti", "numericOrLogical"),
  function(umf, j){
  keep <- process_obs_index(umf, j)
  umf <- methods::callNextMethod(umf, j)
  umf@ylist <- lapply(umf@ylist, function(s) s[,keep,drop=FALSE])
  umf
})

process_obs_index <- function(umf, j){
  if(methods::.hasSlot(umf, "numPrimary")){
    if(umf@numPrimary > 1){
      stop("Can't subset observations in multi-period umf", call.=FALSE)
    }
  }
  if(is.logical(j)){
    stopifnot(length(j) == obsNum(umf))
    keep <- 1:obsNum(umf)
    keep <- keep[j]
  } else if(is.numeric(j)){
    stopifnot(all(abs(j) %in% 1:obsNum(umf)))
    stopifnot(all(j>0) | all(j<0))
    if(all(j<0)){
      keep <- 1:obsNum(umf)
      keep <- keep[!keep %in% abs(j)]
    } else {
      keep <- j
    }
  }
  keep
}


# The following select by period instead of individual occasion----------------

setMethod("subset_obs", c("unmarkedMultFrame", "numericOrLogical"),
  function(umf, j){
  if(obsNum(umf) != numY(umf)){
    stop("Can't subset multinomial observations", call.=FALSE)
  }
  M <- numSites(umf)
  T_old <- umf@numPrimary
  J <- numY(umf) / T_old
  R <- obsNum(umf) / T_old

  # Process selection
  keep <- process_period_index(umf, j)

  umf@numPrimary <- length(keep)

  # Subset y
  keep_per <- rep(1:T_old, each = J)
  y <- lapply(keep, function(x){
    getY(umf)[,keep_per == x,drop=FALSE]
  })
  umf@y <- do.call(cbind, y)

  # Subset yearly site covs
  ysc <- yearlySiteCovs(umf)
  if(!is.null(ysc)){
    keep_site <- rep(1:M, each = T_old)
    ysc_new <- lapply(1:M, function(m){
      ysc_sub <- ysc[keep_site == m,,drop=FALSE]
      ysc_sub <- ysc_sub[keep,,drop=FALSE]
    })
    yearlySiteCovs(umf) <- do.call(rbind, ysc_new)
  }

  # Subset obs covs
  oc <- obsCovs(umf)
  if(!is.null(oc)){
    keep_site <- rep(1:M, each = obsNum(umf))
    keep_per <- rep(1:T_old, each = R)

    oc_new <- lapply(1:M, function(m){
      oc_sub <- oc[keep_site == m,,drop=FALSE]
      oc_sub <- lapply(keep, function(x){
        oc_sub[keep_per == x,,drop=FALSE] 
      })
      do.call(rbind, oc_sub)
    })
    obsCovs(umf) <- do.call(rbind, oc_new)
  }

  # Subset obsToY
  # Only works with non-multinomial
  obsToY(umf) <- diag(numY(umf)) 

  umf
})

setMethod("subset_obs", c("unmarkedFrameGDR", "numericOrLogical"),
  function(umf, j){
  keep <- process_period_index(umf, j)
  T <- umf@numPrimary
  Jdist <- ncol(umf@yDistance) / T
  Jrem <- ncol(umf@yRemoval) / T
  umf <- methods::callNextMethod(umf, j)

  # Subset y matrices
  keep_per <- rep(1:T, each = Jdist)
  yDistance <- lapply(keep, function(k){
    umf@yDistance[,keep_per == k,drop=FALSE]
  })
  umf@yDistance <- do.call(cbind, yDistance)

  keep_per <- rep(1:T, each = Jrem)
  yRemoval <- lapply(keep, function(k){
    umf@yRemoval[,keep_per == k,drop=FALSE]
  })
  umf@yRemoval <- do.call(cbind, yRemoval)
  umf
})

process_period_index <- function(umf, j){
  stopifnot(methods::.hasSlot(umf, "numPrimary"))
  T <- umf@numPrimary
  stopifnot(T > 1)
  if(is.logical(j)){
    stopifnot(length(j) == T)
    keep <- 1:T
    keep <- keep[j]
  } else if(is.numeric(j)){
    stopifnot(all(abs(j) %in% 1:T))
    stopifnot(all(j>0) | all(j<0))
    if(all(j < 0)){
      keep <- 1:T
      keep <- keep[!keep %in% abs(j)]
    } else {
      keep <- j
    }
  }
  keep
}

# head method------------------------------------------------------------------

setMethod("head", "unmarkedFrame", function(x, n) {
  if(missing(n)) n <- 6
  x[1:n,]
})
