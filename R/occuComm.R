unmarkedFrameOccuComm <- function(y, siteCovs=NULL, obsCovs=NULL, speciesCovs=NULL){

  if(is.array(y)){
    # If array, convert to list
    if(length(dim(y)) != 3){
      stop("Incorrect y array dimensions", call.=FALSE)
    }
    S <- dim(y)[3]
    ylist <- lapply(1:S, function(x) y[,,x])
    names(ylist) <- dimnames(y)[[3]]
    y <- ylist
  }

  M <- nrow(y[[1]])
  J <- ncol(y[[1]])
  S <- length(y)
  # Check species covs dimensions
  if(!is.null(speciesCovs)){
    stopifnot(is.list(speciesCovs))
    correct_dims <- sapply(speciesCovs, function(x){
      identical(dim(x), c(M, S)) | identical(dim(x), c(M, J, S)) | identical(length(x), S) 
    })
    bad_dims <- names(speciesCovs)[!correct_dims]
    if(length(bad_dims > 0)){
      stop(paste0("Species covariate(s) ", paste(bad_dims, collapse=", "), " have incorrect dimensions"),
           call.=FALSE)
    }
  }

  if(is.null(names(y))){
    names(y) <- paste0("sp", sprintf("%02d", 1:S))
  }

  obsCovs <- covsToDF(obsCovs, "obsCovs", ncol(y[[1]]), nrow(y[[1]]))
  new("unmarkedFrameOccuComm", y=y[[1]], ylist = y, siteCovs=siteCovs, 
      obsCovs=obsCovs, speciesCovs=speciesCovs, obsToY = diag(J))
}

# Unmarked Frame methods
setMethod("show", "unmarkedFrameOccuComm", function(object)
{
    df <- as(object, "data.frame")
    cat("Data frame representation of unmarkedFrame object.\n")
    cat("Only showing observation matrix for species 1.\n")
    print(df)
})

setMethod("summary", "unmarkedFrameOccuComm", function(object,...) {
    cat("unmarkedFrameOccuComm Object\n\n")
    cat(nrow(object@y), "sites\n")
    cat(length(object@ylist), "species\n")
    cat("Maximum number of observations per site:",obsNum(object),"\n")
    mean.obs <- mean(rowSums(!is.na(object@ylist[[1]])))
    cat("Mean number of observations per site:",round(mean.obs,2),"\n")
    s.one <- sapply(object@ylist,function(x) sum(rowSums(x,na.rm=TRUE)>0))
    cat("Sites with at least one detection, quantiles by species:\n")
    print(quantile(s.one))
    if(!is.null(object@siteCovs)) {
        cat("\nSite-level covariates:\n")
        print(summary(object@siteCovs))
    }
    if(!is.null(object@obsCovs)) {
        cat("\nObservation-level covariates:\n")
        print(summary(object@obsCovs))
    }
    if(!is.null(object@speciesCovs)){
        cat("\nSpecies-level covariates:\n")
        utils::str(object@speciesCovs)
    }
})

setMethod("plot", c(x="unmarkedFrameOccuComm", y="missing"),
	function (x, y, colorkey, ylab="Site", xlab="Observation", ...)
{
    y <- getY(x)
    ym <- max(y, na.rm=TRUE)
    M <- nrow(y)
    J <- ncol(y)
    S <- length(x@ylist)
    y <- as.data.frame(do.call(rbind,x@ylist))
    colnames(y) <- paste("obs",1:J)
    y$site <- rep(1:M,S)
    y$species <- as.factor(rep(names(x@ylist),each=M))
    y2 <- reshape(y, idvar=c("site", "species"), varying=list(1:obsNum(x)),
              v.names="value", direction="long")
    y2$variable <- factor(paste("obs", y2$time), levels=colnames(y))
    if(missing(colorkey))
        colorkey <- list(at=0:(ym+1), labels=list(labels=as.character(0:ym),
            at=(0:ym)+0.5))
    levelplot(value ~ variable*site | species, y2,
        scales=list(relation="free", x=list(labels=1:J)),
        colorkey=colorkey, strip=T, xlab=xlab, ylab=ylab,
        labels=names(x@ylist), ...)
})

process_multispecies_umf <- function(umf, interact_covs){
  ylist <- umf@ylist
  M <- nrow(ylist[[1]])
  J <- ncol(ylist[[1]])
  S <- length(ylist)
  sp <- factor(names(ylist), levels=names(ylist))

  y <- do.call(rbind, ylist)

  sc <- siteCovs(umf)
  sp_rep <- rep(sp, each = M)
  if(is.null(sc)){
    sc <- data.frame(species=sp_rep)
  } else {
    sc <- sc[rep(1:M, S),,drop=FALSE]
    sc[["species"]] <- sp_rep

    # Handle random effects
    interact_covs <- names(sc)[names(sc) %in% interact_covs]
    for (i in interact_covs){
      stopifnot(is.factor(sc[[i]]) | is.character(sc[[i]]))
      sc[[paste0(i, "_species")]] <- factor(paste0(sc[[i]], "_", sc[["species"]]))  
    }
  }

  oc <- obsCovs(umf)
  #sp_rep <- rep(sp, each=M*J)
  if(!is.null(oc)){
    oc <- oc[rep(1:(M*J), S),,drop=FALSE]

    # Handle random effects
    interact_covs <- names(oc)[names(oc) %in% interact_covs]
    for (i in interact_covs){
      stopifnot(is.factor(oc[[i]]) | is.character(oc[[i]]))
      sp_rep <- rep(sc$species, each=M*J)
      oc[[paste0(i, "_species")]] <- factor(paste0(oc[[i]], "_", oc[["species"]]))    
    }
  }

  # Add species level covs
  spc <- umf@speciesCovs
  if(!is.null(spc)){
    # Length S covs
    spc_sp <- sapply(spc, function(x) identical(length(x), S))
    spc_sp <- spc[spc_sp]
    if(length(spc_sp) > 0){
      spc_sp <- lapply(spc_sp, function(x) rep(x, each=M))
      spc_sp <- as.data.frame(spc_sp)
      if(is.null(sc)){
        sc <- spc_sp
      } else {
        sc <- cbind(sc, spc_sp)
      }
    }

    # Site level species covs
    spc_site <- sapply(spc, function(x) identical(dim(x), c(M, S))) 
    spc_site <- spc[spc_site]
    if(length(spc_site) > 0){
      spc_site <- as.data.frame(lapply(spc_site, as.vector))
      if(is.null(sc)){
        sc <- spc_site
      } else {
        sc <- cbind(sc, spc_site)
      }
    }

    spc_obs <- sapply(spc, function(x) identical(dim(x), c(M, J, S)))
    spc_obs <- spc[spc_obs]
    if(length(spc_obs) > 0){
      spc_obs <- lapply(spc_obs, function(x){
        x <- aperm(x, c(2, 1, 3))
        as.vector(x)
      })
      spc_obs <- as.data.frame(spc_obs)
      if(is.null(oc)){
        oc <- spc_obs
      } else {
        oc <- cbind(oc, spc_obs)
      }
    }
  }

  unmarkedFrameOccu(y = y, siteCovs=sc, obsCovs=oc)
}

multispeciesFormula <- function(form, species_covs){
  sf <- split_formula(form)

  covs_add_species <- character(0)

  # Det formula
  det_nobar <- reformulas::nobars(sf[[1]])

  # Remove length S species covs from random part of formula
  species_S <- names(species_covs)[sapply(species_covs, function(x) is.vector(x))]

  trms_obj <- stats::terms(det_nobar)
  trms <- attr(trms_obj, "term.labels")
  if(attr(trms_obj, "intercept")){
    trms <- c("1", trms)
  } else {
    trms <- c("0", trms)
  }
 
  # This is necessary to identify S covariates inside functions like scale()
  trms_S <- sapply(trms, function(x){
    trm_form <- as.formula(paste0("~", x))
    vars <- all.vars(trm_form)
    if(length(vars) == 0) return(FALSE)
    all(vars %in% species_S)
  })
  
  rand_form <- stats::reformulate(trms[!trms_S])

  nterm <- length(trms[!trms_S])
  bars <- ifelse(nterm == 1, "|", "||")
  rand <- paste0("+ (", safeDeparse(rand_form[[2]]), " ", bars, " species)")
  
  barexp <- reformulas::findbars(sf[[1]])
  if(!is.null(barexp)){
    lhs <- lapply(barexp, function(x) x[[2]])
    check_lhs <- sapply(lhs, function(x) x == 1)
    if(any(!check_lhs)) stop("Only random intercepts allowed", call.=FALSE)
    new_rhs <- lapply(barexp, function(x){
      rhs <- safeDeparse(x[[3]])
      covs_add_species <<- c(covs_add_species, rhs)
      new_rhs <- paste0(rhs, "_species")
      x[[3]] <- str2lang(new_rhs)
      x
    })
    rand2 <- lapply(new_rhs, function(x){
      paste0("(", safeDeparse(x), ")")
    })
    rand <- paste(rand, "+", paste(rand2, collapse="+"))
  }

  sf[[1]] <- paste(safeDeparse(det_nobar), rand)

  # State formula
  state_nobar <- reformulas::nobars(sf[[2]])

  # Remove species covs from random part of formula
  trms_obj <- stats::terms(state_nobar)
  trms <- attr(trms_obj, "term.labels")
  if(attr(trms_obj, "intercept")){
    trms <- c("1", trms)
  } else {
    trms <- c("0", trms)
  }

  trms_S <- sapply(trms, function(x){
    trm_form <- as.formula(paste0("~", x))
    vars <- all.vars(trm_form)
    if(length(vars) == 0) return(FALSE)
    all(vars %in% species_S)
  })

  rand_form <- stats::reformulate(trms[!trms_S])

  nterm <- length(trms[!trms_S])
  bars <- ifelse(nterm == 1, "|", "||")
  rand <- paste0("+ (", safeDeparse(rand_form[[2]]), " ", bars, " species)")

  barexp <- reformulas::findbars(sf[[2]])
  if(!is.null(barexp)){
    lhs <- lapply(barexp, function(x) x[[2]])
    check_lhs <- sapply(lhs, function(x) x == 1)
    if(any(!check_lhs)) stop("Only random intercepts allowed", call.=FALSE)
    new_rhs <- lapply(barexp, function(x){
      rhs <- safeDeparse(x[[3]])
      covs_add_species <<- c(covs_add_species, rhs)
      new_rhs <- paste0(rhs, "_species")
      x[[3]] <- str2lang(new_rhs)
      x
    })
    rand2 <- lapply(new_rhs, function(x){
      paste0("(", safeDeparse(x), ")")
    })
    rand <- paste(rand, "+", paste(rand2, collapse="+"))
  }

  sf[[2]] <- paste(safeDeparse(state_nobar), rand)
  
  list(formula = stats::as.formula(str2lang(paste0(sf[[1]], sf[[2]]))),
       covs = unique(covs_add_species))
}

safeDeparse <- function(inp) {
  out <- deparse(inp)
  paste(sapply(out, trimws), collapse=" ")
}

occuComm <- function(formula, data, starts, method="BFGS", se=TRUE, ...){
  newform <- multispeciesFormula(formula, data@speciesCovs)
  newumf <- process_multispecies_umf(data, newform$covs)
  if(missing(starts)) starts <- NULL
  out <- occu(newform$formula, data=newumf, starts = starts, method = method, se = se, ...)
  out@call <- match.call()
  out@formula <- formula
  out@data <- data
  out <- as(out, "unmarkedFitOccuComm")
  out
}

setMethod("ranef_internal", "unmarkedFitOccuComm", function(object){
  S <- length(object@data@ylist)
  M <- numSites(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula, object@data@speciesCovs)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")
  r <- ranef(new_object)
  inds <- split(1:nrow(r@post), rep(1:S, each=M))
  names(inds) <- names(object@data@ylist)
  lapply(inds, function(x){ 
         out <- r
         out@post <- out@post[x,,,drop=FALSE]
         out
  })
})

setMethod("richness", "unmarkedFitOccuComm", 
          function(object, nsims=100, posterior=FALSE){
  S <- length(object@data@ylist)
  M <- numSites(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula, object@data@speciesCovs)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")
  r <- ranef(new_object)
  post <- posteriorSamples(r, nsims=nsims)@samples
  post <- array(post, c(M, S, 1, nsims))
  rich <- apply(post, c(1,3,4), sum, na.rm=TRUE)
  if(!posterior){
    rich <- drop(rich)
    return(apply(rich, 1, mean, na.rm=TRUE))
  }
  new("unmarkedPostSamples", numSites=M, numPrimary=1, nsims=nsims,
      samples=rich)
})

setMethod("getP_internal", "unmarkedFitOccuComm", function(object){
  M <- numSites(object@data)
  J <- obsNum(object@data)
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)
  p <- lapply(p, function(x){
    matrix(x$Predicted, M, J, byrow = TRUE)
  })
  p
})

setMethod("getY_internal", "unmarkedFitOccuComm", function(object) {
            object@data@ylist
})

setMethod("SSE", "unmarkedFitOccuComm", function(fit, ...){
    r <- do.call(rbind, residuals(fit))
    return(c(SSE = sum(r^2, na.rm=T)))
})

setMethod("simulate_internal", "unmarkedFitOccuComm", function(object, nsim){
  S <- length(object@data@ylist)
  M <- numSites(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula, object@data@speciesCovs)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")
  s <- simulate(new_object, nsim = nsim)
  inds <- split(1:nrow(s[[1]]), rep(1:S, each=M))
  names(inds) <- names(object@data@ylist)

  lapply(1:nsim, function(i){
    this_sim <- s[[i]]
    lapply(inds, function(x){
      this_sim[x,,drop=FALSE] # select only one species
    })
  })
})

setMethod("replaceY", "unmarkedFrameOccuComm",
          function(object, newY, replNA=TRUE, ...){
      if(replNA){
        newY <- mapply(function(x, y){ is.na(x) <- is.na(y); x},
                       newY , object@ylist, SIMPLIFY=FALSE)
      }
      object@ylist <- newY
      object
})

setMethod("residual_plot", "unmarkedFitOccuComm", function(x, ...)
{
  r <- do.call(rbind,residuals(x))
  e <- do.call(rbind,fitted(x))
  plot(e, r, ylab = "Residuals", xlab = "Predicted values")
  abline(h = 0, lty = 3, col = "gray")
})

setMethod("summary_internal", "unmarkedFitOccuComm", function(object)
{
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  summaryOut <- summary(object@estimates)
  cat("AIC:", object@AIC,"\n")

  M <- numSites(object@data)
  S <- length(object@data@ylist)
  sr <- object@sitesRemoved
  srmat <- matrix(1:(M*S), M, S)
  srmat[srmat %in% sr] <- NA
  sr <- apply(srmat, 1, function(x) any(is.na(x)))
  sr <- which(sr)

  cat("Number of species:", S)
  cat("\nNumber of sites:", numSites(object@data) - length(sr))
  if(length(sr) > 0){
    cat("\nID of sites removed due to NA:", sr)
  }
  if(!identical(object@opt$convergence, 0L)){
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.", call.=FALSE)
  }

  # Check for potentially bad estimates
  if(!is.null(object@opt$hessian)){
    se <- SE(object)
    has_na <- any(is.na(se)) | any(is.nan(se))
    big_se <- any(abs(se) >= 5)
    if(has_na | big_se){
      warning("Large or missing SE values. Be very cautious using these results.", call.=FALSE)
    }
  }

  nboot <- length(object@bootstrapSamples)
  if(nboot > 0){
    cat("\nBootstrap iterations:", nboot)
  }
  cat("\n\n")
  invisible(summaryOut)
})
