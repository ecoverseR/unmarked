
#  Fit the occupancy model of MacKenzie et al (2002).

occu <- function(formula, data, knownOcc = numeric(0),
                 linkPsi = c("logit", "cloglog"), starts = NULL, method = "BFGS",
                 se = TRUE, engine = c("C", "R", "TMB"), threads=1, ...) {

  # Check arguments------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccu"))
    stop("Data is not an unmarkedFrameOccu object.")

  engine <- match.arg(engine, c("C", "R", "TMB"))
  formulas <- split_formula(formula)
  if(any(sapply(formulas, has_random))) engine <- "TMB"
  if(length(knownOcc)>0 & engine == "TMB"){
    stop("TMB engine does not support knownOcc argument", call.=FALSE)
  }

  linkPsi <- match.arg(linkPsi, c("logit","cloglog"))
  psiLinkFunc <- ifelse(linkPsi=="cloglog", cloglog, plogis)
  psiInvLink <- ifelse(linkPsi=="cloglog", "cloglog", "logistic")

  # Format input data----------------------------------------------------------
  dm <- getDesign(data, formulas)
  y <- truncateToBinary(dm$y)

  # Re-format some variables for C and R engines
  yvec <- as.numeric(t(y))
  navec <- is.na(yvec)
  nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i

  # convert knownOcc to logical so we can correctly to handle NAs.
  knownOccLog <- rep(FALSE, numSites(data))
  knownOccLog[knownOcc] <- TRUE
  if(length(dm$removed.sites)>0) knownOccLog <- knownOccLog[-dm$removed.sites]

  # Set up submodels-----------------------------------------------------------
  state <- unmarkedEstimate("Occupancy", short.name = "psi", invlink=psiInvLink)
  det <- unmarkedEstimate("Detection", short.name = "p", invlink="logistic")
  estimateList <- unmarkedEstimateList(list(state=state,det=det))

  # Set up parameter names and indices-----------------------------------------
  par_inds <- get_parameter_inds(estimateList, dm)

  # Set up negative log likelihood functions for C++ and R engines-----_-------
  if(identical(engine, "C")) {
    nll <- function(params) {
      beta.psi <- params[par_inds$state]
      beta.p <- params[par_inds$det]
      nll_occu(
        yvec, dm$X_state, dm$X_det, beta.psi, beta.p, nd, knownOccLog, navec,
        dm$offset_state, dm$offset_det, linkPsi
      )
    }
  } else if (identical(engine, "R")){

    J <- ncol(y)
    M <- nrow(y)

    nll <- function(params) {
      psi <- psiLinkFunc(dm$X_state %*% params[par_inds$state] + dm$offset_state)
      psi[knownOccLog] <- 1
      pvec <- plogis(dm$X_det %*% params[par_inds$det] + dm$offset_det)
      cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
      cp[navec] <- 1 # so that NA's don't modify likelihood
      cpmat <- matrix(cp, M, J, byrow = TRUE) #
      loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
      -sum(loglik)
    }
  }

  # Fit model with C++ and R engines-------------------------------------------
  if(engine %in% c("C", "R")){

    fit <- fit_optim(nll, starts, method, se, estimateList, par_inds, ...)

  # Fit model with TMB engine--------------------------------------------------
  } else if(identical(engine, "TMB")){

    # Set up TMB input data
    tmb_inputs <- get_TMB_inputs(formulas, dm, par_inds, data, 
                                 no_detect=nd, link=ifelse(linkPsi=="cloglog",1,0))

    # Fit model with TMB
    fit <- fit_TMB("tmb_occu", starts, method, estimateList, par_inds,
                    tmb_inputs, data, ...)
  }

  # Create unmarkedFit object--------------------------------------------------
  new("unmarkedFitOccu", fitType = "occu", call = match.call(),
      formula = formula, formlist = formulas, data = data,
      sitesRemoved = dm$removed.sites,
      estimates = fit$estimate_list, AIC = fit$AIC, opt = fit$opt,
      negLogLike = fit$opt$value,
      nllFun = fit$nll, knownOcc = knownOccLog, TMB=fit$TMB)
}
