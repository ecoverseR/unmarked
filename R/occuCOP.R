# Fit the occupancy model COP
# (Counting Occurrences Process)

# Occupancy
# z_i ~ Bernoulli(psi_i)
# 
# with:
#   z_i = Occupancy state of site i
#       = 1 if the site i is occupied
#       = 0 else
#   psi_i = Occupancy probability of site i

# Detection
# N_ij | z_i = 1 ~ Poisson(lambda_ij*L_ij)
# N_ij | z_i = 0 ~ 0
# 
# with:
#   N_ij = Number of detections of site i during observation j
#   z_i = Occupancy state of site i
#   lambda_ij = Detection rate of the observation j in site i
#   L_ij = Length/Duration of the observation j in site i

# Methods ----------------------------------------------------------------------

## getDesign method ----
# All this method does is translate the submodel names
setMethod("getDesign", "unmarkedFrameOccuCOP", function(umf, formulas, na.rm = TRUE){
  names(formulas)[names(formulas) == "psi"] <- "state"
  names(formulas)[names(formulas) == "lambda"] <- "det"
  out <- methods::callNextMethod(umf, formulas = formulas, na.rm = na.rm)
  nms <- gsub("state", "psi", names(out))
  nms <- gsub("det", "lambda", nms)
  names(out) <- nms
  out
})


## getL method ----
setMethod("getL", "unmarkedFrameOccuCOP", function(object) {
  return(object@L)
})


## show method ----
setMethod("show", "unmarkedFrameOccuCOP", function(object) {
  J <- ncol(object@L)
  df_unmarkedFrame <- as(object, "data.frame")
  df_L <- data.frame(object@L)
  colnames(df_L) <- paste0("L.", 1:J)
  if (ncol(df_unmarkedFrame) > J) {
    df <- cbind(df_unmarkedFrame[, 1:J, drop = FALSE], 
                df_L, 
                df_unmarkedFrame[, (J + 1):ncol(df_unmarkedFrame), drop = FALSE])
  } else {
    df <- cbind(df_unmarkedFrame[, 1:J], 
                df_L)
  }
  cat("Data frame representation of unmarkedFrame object.\n")
  print(df)
})
# LP note: as is defined in unmarkedFrame.R part "COERCION"


## summary method ----
setMethod("summary", "unmarkedFrameOccuCOP", function(object,...) {
  cat("unmarkedFrameOccuCOP Object\n\n")
  
  cat(nrow(object@y), "sites\n")
  cat("Maximum number of sampling occasions per site:",obsNum(object),"\n")
  mean.obs <- mean(rowSums(!is.na(getY(object))))
  cat("Mean number of sampling occasions per site:",round(mean.obs,2),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  
  cat("Tabulation of y observations:")
  print(table(object@y, exclude=NULL))
  
  if(!is.null(object@L)) {
    cat("\nTabulation of sampling occasions length:")
    print(table(object@L))
  }

  if(!is.null(object@siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object@siteCovs))
  }
  
  if(!is.null(object@obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object@obsCovs))
  }
})


## getP ----
setMethod("getP_internal", "unmarkedFitOccuCOP", function(object) {
  data <- object@data
  M = nrow(getY(data))
  J = ncol(getY(data))
  lam <- predict(object, type="lambda", level=NULL, na.rm=FALSE)$Predicted
  lam <- matrix(lam, M, J, byrow=TRUE)
  lam
})


## plot ----
setMethod("residual_plot", "unmarkedFitOccuCOP", function(x, ...) {
  y <- getY(x)
  r <- residuals(x)
  e <- fitted(x)
  
  old_graph <- graphics::par("mfrow", "mar")
  on.exit(graphics::par(mfrow = old_graph$mfrow, mar = old_graph$mar))
  
  {
    graphics::par(mfrow = c(2, 1))
    graphics::par(mar = c(0, 5, 2, 2))
    plot(e, y,
         ylab = "Observed data",
         xlab = "Predicted data",
         xaxt = 'n')
    abline(a = 0, b = 1, lty = 3, col = "red")
    title(main = "COP model - detection events count", outer = F)
    
    graphics::par(mar = c(4, 5, 0.5, 2))
    plot(e, r,
         ylab = "Residuals",
         xlab = "Predicted data")
    abline(h = 0, lty = 3, col = "red")
  }
})


## simulate ----
setMethod("simulate_internal", "unmarkedFitOccuCOP",
  function(object, nsim){
  y <- object@data@y
  M <- nrow(y)
  J <- ncol(y)
  
  # Occupancy probability psi depending on the site covariates
  psi <- predict(object, type = "psi", level = NULL, na.rm=FALSE)$Predicted
  
  # Detection rate lambda depending on the observation covariates
  lambda = getP(object = object)
  
  # Simulations
  simList <- vector("list", nsim)
  for(i in 1:nsim) {
    Z <- rbinom(M, 1, psi)
    # Z[object@knownOcc] <- 1
    # TODO: should Z be replicated J times here?
    y = matrix(rpois(n = M * J, lambda = as.numeric(t(lambda))),
               nrow = M, ncol = J, byrow = T) * Z
    simList[[i]] <- y
  }
  return(simList)
})

setMethod("get_fitting_function", "unmarkedFrameOccuCOP",
          function(object, model, ...){
  occuCOP
})

## nonparboot ----
setMethod("nonparboot_internal", "unmarkedFitOccuCOP",
  function(object, B, keepOldSamples) {
  stop("Not currently supported for unmarkedFitOccuCOP", call.=FALSE)
})


## ranef ----
setMethod("ranef_internal", "unmarkedFitOccuCOP", function(object, ...) {
  # Sites removed (srm) and sites kept (sk)
  srm <- object@sitesRemoved
  sk <- 1:numSites(getData(object))
  if (length(srm) > 0) {
    sk = sk[-srm]
  }

  # unmarkedFrame informations
  M <- length(sk)
  J <- obsNum(getData(object))
  y <- getY(getData(object))[sk,]
  
  # Estimated parameters
  psi <- predict(object, type = "psi")[sk, 1]
  lambda <- getP(object)[sk,]
  
  # Estimate posterior distributions
  z = c(0, 1)
  post <- array(0, c(M, 2, 1), dimnames = list(NULL, z))
  for (i in 1:M) {
    # psi density
    f <- dbinom(x = z,
                size = 1,
                prob = psi[i])
    
    # lambda density
    g <- c(1, 1)
    for (j in 1:J) {
      if (is.na(y[i, j]) | is.na(lambda[i, j])) {
        next
      }
      g = g * dpois(x = y[i, j], lambda = lambda[i, j] * z)
    }
    
    # psi*lambda density
    fudge <- f * g
    post[i, , 1] <- fudge / sum(fudge)
  }
  
  new("unmarkedRanef", post = post)
})


# Used by update() method
setMethod("rebuild_call", "unmarkedFitOccuCOP", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["psiformula"]] <- object@formlist$psi
  cl[["lambdaformula"]] <- object@formlist$lambda
  cl
})

# Useful functions -------------------------------------------------------------

check.integer <- function(x, na.ignore = F) {
  if (is.na(x) & na.ignore) {
    return(T)
  }
  x %% 1 == 0
}

# unmarkedFrame ----------------------------------------------------------------

unmarkedFrameOccuCOP <- function(y, L, siteCovs = NULL, obsCovs = NULL) {
  
  # Verification that they are non-NA data in y
  if (all(is.na(y))) {
    stop("y only contains NA. y needs to contain non-NA integers.")
  }
  
  # Verification that these are count data (and not detection/non-detection)
  if (max(y, na.rm = T) == 1) {
    warning("unmarkedFrameOccuCOP is for count data. ",
            "The data furnished appear to be detection/non-detection.")
  }
  
  # Number of sampling occasions
  J <- ncol(y)
  
  # If missing L: replace by matrix of 1
  # and the lambda will be the detection rate per observation length
  if (missing(L)) {
    L <- y * 0 + 1
    warning("L is missing, replacing it by a matrix of 1.")
  } else if (is.null(L)) {
    L <- y * 0 + 1
    warning("L is missing, replacing it by a matrix of 1.")
  }
  
  # Transformation observation covariates
  obsCovs <- covsToDF(
    covs = obsCovs,
    name = "obsCovs",
    obsNum = J,
    numSites = nrow(y)
  )
  
  # Create S4 object of class unmarkedFrameOccuCOP
  umf <- new(
    Class = "unmarkedFrameOccuCOP",
    y = y,
    L = L,
    siteCovs = siteCovs,
    obsCovs = obsCovs,
    obsToY = diag(J)
  )
  
  return(umf)
}


# occuCOP ----------------------------------------------------------------------

occuCOP <- function(data,
                    psiformula = ~1,
                    lambdaformula = ~1,
                    psistarts,
                    lambdastarts,
                    starts,
                    method = "BFGS",
                    se = TRUE,
                    engine = c("C", "R"),
                    na.rm = TRUE,
                    return.negloglik = NULL,
                    L1 = FALSE,
                    ...) {
  #TODO: random effects
  
  # Neg loglikelihood COP ------------------------------------------------------
  R_nll_occuCOP <- function(params) {
    
    # Reading and transforming params
    # Taking into account the covariates
    
    # psi as a function of covariates
    #   psi in params are transformed with a logit transformation (qlogis)
    #       so they're back-transformed to a proba with inverse logit (plogis)
    #   Xpsi is the matrix with occupancy covariates
    #   params is the vector with all the parameters
    #   psiIdx is the index of Occupancy Parameters in params
    psi <- plogis(Xpsi %*% params[psiIdx])
    
    # lambda as a function of covariates
    #   lambda in params are transformed with a log-transformation
    #          so they're back-transformed to a rate with exp here
    #   Xlambda is the matrix with detection covariates
    #   params is the vector with all the parameters
    #   lambdaIdx is the index of Occupancy Parameters in params
    lambda <- exp(Xlambda %*% params[lambdaIdx])
    
    # Listing sites analysed = sites not removed (due to NAs)
    if (length(sitesRemoved) > 0) {
      siteAnalysed = (1:M)[-sitesRemoved]
    } else {
      siteAnalysed = (1:M)
    }
    
    # Probability for each site (i)
    iProb <- rep(NA, M)
    
    for (i in siteAnalysed) {
      # iIdx is the index to access the vectorised vectors of all obs in site i
      iIdxall <- ((i - 1) * J + 1):(i * J)
      
      # Removing NAs
      iIdx = iIdxall[!removed_obsvec[iIdxall]]
      
      if (SitesWithDetec[i]) {
        # If there is at least one detection in site i
        iProb[i] = psi[i] * (
          (sum(lambda[iIdx] * Lvec[iIdx])) ^ sum(yvec[iIdx]) / 
            factorial(sum(yvec[iIdx])) *
            exp(-sum(lambda[iIdx] * Lvec[iIdx]))
        )
        
      } else {
        # If there is zero detection in site i
        iProb[i] = psi[i] * exp(-sum(lambda[iIdx] * Lvec[iIdx])) + (1 - psi[i]) 
        
      }
    }
    
    # log-likelihood
    ll = sum(log(iProb[siteAnalysed]))
    return(-ll)
  }
  
  # Check arguments ------------------------------------------------------------
  if (!is(data, "unmarkedFrameOccuCOP")) {
    stop("Data is not an unmarkedFrameOccuCOP object. See ?unmarkedFrameOccuCOP if necessary.")
  }
  stopifnot(class(psiformula) == "formula")
  stopifnot(class(lambdaformula) == "formula")
  if(!missing(psistarts)){stopifnot(class(psistarts) %in%  c("numeric", "double", "integer"))}
  if(!missing(lambdastarts)){stopifnot(class(lambdastarts) %in%  c("numeric", "double", "integer"))}
  se = as.logical(match.arg(
    arg = as.character(se),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  na.rm = as.logical(match.arg(
    arg = as.character(na.rm),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  engine <- match.arg(engine, c("C", "R"))
  L1 = as.logical(match.arg(
    arg = as.character(L1),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  
  
  # Do not yet manage random effects!!!
  if (has_random(psiformula) | has_random(lambdaformula)) {
    stop("occuCOP does not currently handle random effects.")
  }
  
  # Format input data ----------------------------------------------------------
  
  # Retrieve formlist
  formlist <- list(lambda = lambdaformula, psi = psiformula)
  
  # Get the design matrix (calling the getDesign method for unmarkedFrame)
  # For more informations, see: getMethod("getDesign", "unmarkedFrameOccuCOP")
  designMats <- getDesign(umf = data, formulas = formlist, na.rm = TRUE)
  
  # y is the count detection data (matrix of size M sites x J observations)
  y <- designMats$y
  
  # L is the length of observations (matrix of size M sites x J observations)
  L <- getL(data)
  L <- L[!1:nrow(L) %in% designMats$removed.sites,,drop=FALSE]
  if (!L1) {
    if (!any(is.na(L))) {
      if (all(L == 1)) {
        warning(
          "All observations lengths (L) are set to 1. ",
          "If they were not user-defined, lambda corresponds to the ",
          "detection rate multiplied by the observation length, ",
          "not just the detection rate per time-unit or space-unit.\n",
          "You can remove this warning by adding 'L1=TRUE' in the function inputs."
        )
      }
    }
  }
  
  # Xpsi is the fixed effects design matrix for occupancy
  Xpsi <- designMats$X_psi
  
  # Xlambda is the fixed effects design matrix for detection rate
  Xlambda <- designMats$X_lambda
  
  # Zpsi is the random effects design matrix for occupancy
  Zpsi <- designMats$Z_psi
  
  # Zlambda is the random effects design matrix for detection rate
  Zlambda <- designMats$Z_lambda
  
  # removed_obs is a M x J matrix of the observations removed from the analysis
  removed_obs <- is.na(y)
  sitesRemoved <- designMats$removed.sites
  
  # Number of sites
  M <- nrow(y)
  
  # Number of sampling occasions
  J <- ncol(y)
  
  # Total number of detection per site
  NbDetecPerSite = rowSums(y, na.rm=T)
  
  # Sites where there was at least one detection
  SitesWithDetec = NbDetecPerSite > 0
  
  # Number of sites where there was at least one detection
  NbSitesWithDetec = sum(SitesWithDetec)
  
  # Set up parameter names and indices-----------------------------------------
  
  # ParamPsi Occupancy parameter names
  ParamPsi <- colnames(Xpsi)
  
  # ParamLambda Detection parameter names
  ParamLambda <- colnames(Xlambda)
  
  # NbParamPsi Number of occupancy parameters
  NbParamPsi <- ncol(Xpsi)
  
  # NbParamLambda Number of detection parameters
  NbParamLambda <- ncol(Xlambda)
  
  # nP Total number of parameters
  nP <- NbParamPsi + NbParamLambda
  
  # psiIdx Index of the occupancy parameters
  psiIdx <- 1:NbParamPsi
  
  # lambdaIdx Index of the detection parameters
  lambdaIdx <- (NbParamPsi+1):nP
  
  # Re-format some variables for C and R engines
  #   From Matrix of dim MxJ to vector of length MxJ:
  #   c(ySite1Obs1, ySite1Obs2, ..., ySite1ObsJ, ysite2Obs1, ...)
  yvec <- as.numeric(t(y))
  Lvec <- as.numeric(t(L))
  removed_obsvec <- as.logical(t(removed_obs))
  
  # return.negloglik -----------------------------------------------------------
  if (!is.null(return.negloglik)) {
    df_NLL = data.frame(t(as.data.frame(return.negloglik)))
    rownames(df_NLL) = NULL
    colnames(df_NLL) = c(paste0("logit(psi).", ParamPsi),
                         paste0("log(lambda).", ParamLambda))
    df_NLL$nll = NA
    for (i in 1:nrow(df_NLL)) {
      df_NLL$nll[i] = R_nll_occuCOP(params = as.numeric(as.vector(df_NLL[i, -ncol(df_NLL)])))
    }
    return(df_NLL)
  }
  
  # nll function depending on engine -------------------------------------------
  if (identical(engine, "C")) {
    nll <- function(params) {
      nll_occuCOP(
        y = yvec,
        L = Lvec,
        Xpsi = Xpsi,
        Xlambda = Xlambda,
        beta_psi = params[psiIdx],
        beta_lambda = params[lambdaIdx],
        removed = removed_obsvec
      )
    }
  } else if (identical(engine, "R")) {
    nll <- R_nll_occuCOP
  }
  
  
  # Optimisation ---------------------------------------------------------------
  
  ## Checking the starting point for optim ----
  # Check if either (psistarts AND lambdastarts) OR starts is provided
  if (!missing(psistarts) & !missing(lambdastarts)) {
    # Both psistarts and lambdastarts provided
    if (!missing(starts)){
      if (!all(c(psistarts, lambdastarts) == starts)) {
        warning(
          "You provided psistarts, lambdastarts and starts. ",
          "Please provide either (psistarts AND lambdastarts) OR starts. ",
          "Using psistarts and lambdastarts."
        )
      }
    }
    if (length(lambdastarts) != NbParamLambda) {
      stop("lambdastarts (", paste(lambdastarts, collapse = ", "), ") ",
           "should be of length ", NbParamLambda, " with lambdaformula ", lambdaformula)
    }
    if (length(psistarts) != NbParamPsi) {
      stop("psistarts (", paste(psistarts, collapse = ", "), ") ",
           "should be of length ", NbParamPsi, " with psiformula ", psiformula)
    }
    starts <- c(psistarts, lambdastarts)
  } else if (!missing(starts)) {
    # starts provided
    if (length(starts) != nP) {
      stop("starts (", paste(starts, collapse = ", "), ") ",
           "should be of length ", nP, 
           " with psiformula ", psiformula, 
           " and lambdaformula ", lambdaformula)
    }
    
    psistarts <- starts[1:NbParamPsi]
    lambdastarts <- starts[(NbParamPsi + 1):(NbParamPsi + NbParamLambda)]
    
  } else {
    # No arguments provided, apply default values
    
    if (missing(lambdastarts)) {
      # If lambda starts argument was not given:
      # 0 by default 
      # so lambda = exp(0) = 1 by default
      lambdastarts = rep(0, NbParamLambda)
    } else if (length(lambdastarts) != NbParamLambda) {
      stop("lambdastarts (", paste(lambdastarts, collapse = ", "), ") ",
           "should be of length ", NbParamLambda, " with lambdaformula ", lambdaformula)
    }
    
    if (missing(psistarts)) {
      # If psi starts argument was not given
      # 0 by default 
      # so psi = plogis(0) = 0.5 by default
      psistarts = rep(0, NbParamPsi)
    } else if (length(psistarts) != NbParamPsi) {
      stop("psistarts (", paste(psistarts, collapse = ", "), ") ",
           "should be of length ", NbParamPsi, " with psiformula ", psiformula)
    }
    
    starts <- c(psistarts, lambdastarts)
  }
  
  ## Run optim ----
  opt <- optim(
    starts,
    nll,
    method = method,
    hessian = se,
    ...
  )
  
  # Get output -----------------------------------------------------------------
  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  tmb_mod <- NULL
  fmAIC <- 2 * opt$value + 2 * nP
  
  # Organize effect estimates
  names(ests) <- c(ParamPsi, ParamLambda)
  psi_coef <- list(ests = ests[psiIdx], cov = as.matrix(covMat[psiIdx, psiIdx]))
  lambda_coef <- list(ests = ests[lambdaIdx],
                      cov = as.matrix(covMat[lambdaIdx, lambdaIdx]))
  
  # No random effects
  psi_rand_info <- lambda_rand_info <- list()
  
  # Create unmarkedEstimates ---------------------------------------------------
  psi_uE <- unmarkedEstimate(
    name = "Occupancy probability",
    short.name = "psi",
    estimates = psi_coef$ests,
    covMat = psi_coef$cov,
    fixed = 1:NbParamPsi,
    invlink = "logistic",
    invlinkGrad = "logistic.grad",
    randomVarInfo = psi_rand_info
  )
  
  lambda_uE <- unmarkedEstimate(
    name = "Detection rate",
    short.name = "lambda",
    estimates = lambda_coef$ests,
    covMat = lambda_coef$cov,
    fixed = 1:NbParamLambda,
    invlink = "exp",
    invlinkGrad = "exp",
    randomVarInfo = lambda_rand_info
  )
  
  estimateList <- unmarkedEstimateList(list(psi = psi_uE, lambda = lambda_uE))
  
  # Create unmarkedFit object--------------------------------------------------
  umfit <- new(
    "unmarkedFitOccuCOP",
    fitType = "occuCOP",
    call = match.call(),
    formlist = formlist,
    data = data,
    estimates = estimateList,
    sitesRemoved = sitesRemoved,
    removed_obs = removed_obs,
    AIC = fmAIC,
    opt = opt,
    negLogLike = opt$value,
    nllFun = nll,
    TMB = tmb_mod
  )
  
  return(umfit)
}
