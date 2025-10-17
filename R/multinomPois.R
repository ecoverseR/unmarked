# Fit the multinomial-Poisson abundance mixture model.

multinomPois <- function(formula, data, starts, method = "BFGS",
    se = TRUE, engine = c("C","R","TMB"), ...)
{

    if(!is(data, "unmarkedFrameMPois"))
		    stop("Data is not a data frame or unmarkedFrame.")
    engine <- match.arg(engine, c("C", "R", "TMB"))
    formulas <- split_formula(formula)
    names(formulas) <- c("det", "state")
    if(any(sapply(formulas, has_random))) engine <- "TMB"
    dm <- getDesign(data, formulas)
    y <- dm$y

    J <- ncol(y)
    R <- obsNum(data)
    M <- nrow(y)
    piFun <- data@piFun

    lamParms <- colnames(dm$X_state)
    detParms <- colnames(dm$X_det)
    nDP <- ncol(dm$X_det)
    nAP <- ncol(dm$X_state)
    lamIdx <- 1:nAP
    pIdx <- (nAP+1):(nAP+nDP)
    nP <- nDP + nAP
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec <- as.numeric(y)
    navec <- is.na(yvec)

    nll_R <- function(parms) {
        lambda <- exp(dm$X_state %*% parms[1 : nAP] + dm$offset_state)
        p <- plogis(dm$X_det %*% parms[(nAP + 1) : nP] + dm$offset_det)
        p.matrix <- matrix(p, M, R, byrow = TRUE)
        pi <- do.call(piFun, list(p = p.matrix))
        logLikeSite <- dpois(y, matrix(lambda, M, J) * pi, log = TRUE)
        logLikeSite[navec] <- 0
        -sum(logLikeSite)
        }

    nll_C <- function(params) {
        nll_multinomPois(
            params,piFun,
            dm$X_state, dm$offset_state, dm$X_det, dm$offset_det,
            yC, navecC, nP,nAP
        )
    }

    if(engine=="R"){
      nll <- nll_R
    } else if(engine=="C"){
      yC <- as.numeric(t(y))
      navecC <- is.na(yC)
      nll <- nll_C
      if(!piFun%in%c('doublePiFun','removalPiFun','depDoublePiFun')){
        warning("Custom pi functions are not supported by C engine. Using R engine instead.")
        nll <- nll_R
      }
    }

    if(engine %in% c("C","R")){
      if(missing(starts)) starts <- rep(0, nP)
      fm <- optim(starts, nll, method = method, hessian = se, ...)
      covMat <- invertHessian(fm, nP, se)
      ests <- fm$par
      names(ests) <- c(lamParms, detParms)
      fmAIC <- 2 * fm$value + 2 * nP
      tmb_mod <- NULL

      # Organize fixed-effect estimates
      state_coef <- list(ests=ests[lamIdx], cov=as.matrix(covMat[lamIdx,lamIdx]))
      det_coef <- list(ests=ests[pIdx], cov=as.matrix(covMat[pIdx, pIdx]))

      # No random effects in C or R engines
      state_rand_info <- det_rand_info <- list()

    } else if(engine == "TMB"){

      obs_all <- add_covariates(obsCovs(data), siteCovs(data), numSites(data)*obsNum(data))
      inps <- get_ranef_inputs(formulas, list(det=obs_all, state=siteCovs(data)),
                               list(dm$X_det, dm$X_state), dm[c("Z_det","Z_state")])

      if(!piFun%in%c('doublePiFun','removalPiFun','depDoublePiFun')){
        stop("Custom pi functions are not supported by TMB engine.")
      }
      pifun_type <- switch(piFun, removalPiFun={0}, doublePiFun={1},
                           depDoublePiFun={2})
      tmb_dat <- c(list(y=y, pifun_type=pifun_type, offset_state=dm$offset_state,
                        offset_det=dm$offset_det), inps$data)

      # Fit model in TMB
      if(missing(starts)) starts <- NULL
      tmb_out <- fit_TMB("tmb_multinomPois", tmb_dat, inps$pars, inps$rand_ef,
                         starts=starts, method, ...)
      tmb_mod <- tmb_out$TMB
      fm <- tmb_out$opt
      fmAIC <- tmb_out$AIC
      nll <- tmb_mod$fn

      # Organize fixed-effect estimate from TMB output
      state_coef <- get_coef_info(tmb_out$sdr, "state", lamParms, lamIdx)
      det_coef <- get_coef_info(tmb_out$sdr, "det", detParms, pIdx)

      # Organize random-effect estimates from TMB output
      state_rand_info <- get_randvar_info(tmb_out$sdr, "state", formulas$state, siteCovs(data))
      det_rand_info <- get_randvar_info(tmb_out$sdr, "det", formulas$det, obs_all)

    }

    stateEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
                        estimates = state_coef$ests, covMat = state_coef$cov,
                        fixed=1:nAP, invlink = "exp", invlinkGrad = "exp",
                        randomVarInfo=state_rand_info)

    detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
                        estimates = det_coef$ests, covMat = det_coef$cov,
                        fixed = 1:nDP, invlink = "logistic", invlinkGrad = "logistic.grad",
                        randomVarInfo=det_rand_info)

    estimateList <- unmarkedEstimateList(list(state=stateEstimates,
        det=detEstimates))

    umfit <- new("unmarkedFitMPois", fitType = "multinomPois",
        call = match.call(), formula = formula, formlist = formulas, data = data,
        estimates = estimateList, sitesRemoved = dm$removed.sites,
        AIC = fmAIC, opt = fm, negLogLike = fm$value, nllFun = nll, TMB=tmb_mod)

    return(umfit)
}
