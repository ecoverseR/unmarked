
#  Fit the occupancy model of MacKenzie et al (2002).

occu <- function(formula, data, knownOcc = numeric(0),
                 linkPsi = c("logit", "cloglog"), starts, method = "BFGS",
                 se = TRUE, engine = c("C", "R", "TMB"), threads=1, ...) {

  # Check arguments------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccu"))
    stop("Data is not an unmarkedFrameOccu object.")

  engine <- match.arg(engine, c("C", "R", "TMB"))
  formulas <- split_formula(formula)
  names(formulas) <- c("det", "state")
  if(any(sapply(formulas, has_random))) engine <- "TMB"
  if(length(knownOcc)>0 & engine == "TMB"){
    stop("TMB engine does not support knownOcc argument", call.=FALSE)
  }

  linkPsi <- match.arg(linkPsi, c("logit","cloglog"))
  psiLinkFunc <- ifelse(linkPsi=="cloglog", cloglog, plogis)
  psiInvLink <- ifelse(linkPsi=="cloglog", "cloglog", "logistic")
  psiLinkGrad <- ifelse(linkPsi=="cloglog", "cloglog.grad", "logistic.grad")

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

  # Set up parameter names and indices-----------------------------------------
  occParms <- colnames(dm$X_state)
  detParms <- colnames(dm$X_det)
  nDP <- ncol(dm$X_det)
  nOP <- ncol(dm$X_state)
  nP <- nDP + nOP
  psiIdx <- 1:nOP
  pIdx <- (nOP+1):nP

  # Set up negative log likelihood functions for C++ and R engines-----_-------
  if(identical(engine, "C")) {
    nll <- function(params) {
      beta.psi <- params[1:nOP]
      beta.p <- params[(nOP+1):nP]
      nll_occu(
        yvec, dm$X_state, dm$X_det, beta.psi, beta.p, nd, knownOccLog, navec,
        dm$offset_state, dm$offset_det, linkPsi
      )
    }
  } else if (identical(engine, "R")){

    J <- ncol(y)
    M <- nrow(y)

    nll <- function(params) {
      psi <- psiLinkFunc(dm$X_state %*% params[1 : nOP] + dm$offset_state)
      psi[knownOccLog] <- 1
      pvec <- plogis(dm$X_det %*% params[(nOP + 1) : nP] + dm$offset_det)
      cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
      cp[navec] <- 1 # so that NA's don't modify likelihood
      cpmat <- matrix(cp, M, J, byrow = TRUE) #
      loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
      -sum(loglik)
    }
  }

  # Fit model with C++ and R engines-------------------------------------------
  if(engine %in% c("C", "R")){
    if(missing(starts)) starts <- rep(0, nP)
    if(length(starts) != nP){
      stop(paste("The number of starting values should be", nP))
    }
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    covMat <- invertHessian(fm, nP, se)
    ests <- fm$par
    tmb_mod <- NULL
    fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)

    # Organize effect estimates
    names(ests) <- c(occParms, detParms)
    state_coef <- list(ests=ests[1:nOP], cov=as.matrix(covMat[1:nOP,1:nOP]))
    det_coef <- list(ests=ests[(nOP+1):nP],
                     cov=as.matrix(covMat[(nOP+1):nP, (nOP+1):nP]))

    # No random effects for C++ and R engines
    state_rand_info <- det_rand_info <- list()

  # Fit model with TMB engine--------------------------------------------------
  } else if(identical(engine, "TMB")){

    # Set up TMB input data
    obs_all <- add_covariates(obsCovs(data), siteCovs(data), length(getY(data)))
    inps <- get_ranef_inputs(formulas, list(det=obs_all, state=siteCovs(data)),
                             list(dm$X_det, dm$X_state), dm[c("Z_det","Z_state")])

    tmb_dat <- c(list(y=y, no_detect=nd, link=ifelse(linkPsi=="cloglog",1,0),
                      offset_state=dm$offset_state, offset_det=dm$offset_det), inps$data)

    # Fit model with TMB
    if(missing(starts)) starts <- NULL
    tmb_out <- fit_TMB("tmb_occu", tmb_dat, inps$pars, inps$rand_ef,
                        starts=starts, method, ...)
    tmb_mod <- tmb_out$TMB
    fm <- tmb_out$opt
    fmAIC <- tmb_out$AIC
    nll <- tmb_mod$fn

    # Organize fixed effect estimates
    state_coef <- get_coef_info(tmb_out$sdr, "state", occParms, 1:nOP)
    det_coef <- get_coef_info(tmb_out$sdr, "det", detParms, (nOP+1):nP)

    # Organize random effect estimates
    state_rand_info <- get_randvar_info(tmb_out$sdr, "state", formulas$state,
                                        siteCovs(data))
    det_rand_info <- get_randvar_info(tmb_out$sdr, "det", formulas$det,
                                      obs_all)

    }

  # Create unmarkedEstimates---------------------------------------------------
  state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                            estimates = state_coef$est,
                            covMat = state_coef$cov,
                            fixed = 1:nOP,
                            invlink = psiInvLink,
                            invlinkGrad = psiLinkGrad,
                            randomVarInfo=state_rand_info
                            )

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates =det_coef$est,
                          covMat = det_coef$cov,
                          fixed = 1:nDP,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad",
                          randomVarInfo=det_rand_info
                          )

  estimateList <- unmarkedEstimateList(list(state=state, det=det))

  # Create unmarkedFit object--------------------------------------------------
  umfit <- new("unmarkedFitOccu", fitType = "occu", call = match.call(),
                 formula = formula, formlist = formulas, data = data,
                 sitesRemoved = dm$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = fm,
                 negLogLike = fm$value,
                 nllFun = nll, knownOcc = knownOccLog, TMB=tmb_mod)

  return(umfit)
}
