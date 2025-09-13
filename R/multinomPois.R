# Fit the multinomial-Poisson abundance mixture model.

multinomPois <- function(formula, data, starts = NULL, method = "BFGS",
    se = TRUE, engine = c("C","R","TMB"), ...){

  if(!is(data, "unmarkedFrameMPois"))
		  stop("Data is not a data frame or unmarkedFrame.")
  engine <- match.arg(engine, c("C", "R", "TMB"))
  formulas <- split_formula(formula)
  if(any(sapply(formulas, has_random))) engine <- "TMB"
  dm <- getDesign(data, formulas)
  y <- dm$y

  J <- ncol(y)
  R <- obsNum(data)
  M <- nrow(y)
  piFun <- data@piFun

  # Set up submodels
  state <- unmarkedEstimate("Abundance", short.name = "lambda", invlink="exp")
  det <- unmarkedEstimate("Detection", short.name = "p", invlink="logistic")
  estimateList <- unmarkedEstimateList(list(state=state,det=det))
    
  # Set up parameter names and indices
  par_inds <- get_parameter_inds(estimateList, dm)

  yvec <- as.numeric(y)
  navec <- is.na(yvec)

  nll_R <- function(parms) {
      lambda <- exp(dm$X_state %*% parms[par_inds$state] + dm$offset_state)
      p <- plogis(dm$X_det %*% parms[par_inds$det] + dm$offset_det)
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
      yC, navecC, length(unlist(par_inds)), length(par_inds$state)
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
    
    fit <- fit_optim(nll, starts, method, se, estimateList, par_inds, ...)

  } else if(engine == "TMB"){
    # Set up TMB input data
    if(!piFun%in%c('doublePiFun','removalPiFun','depDoublePiFun')){
      stop("Custom pi functions are not supported by TMB engine.")
    }
    pifun_type <- switch(piFun, removalPiFun={0}, doublePiFun={1},
                         depDoublePiFun={2})
    tmb_inputs <- get_TMB_inputs(formulas, dm, par_inds, data, 
                                 pifun_type=pifun_type)

    # Fit model with TMB
    fit <- fit_TMB2("tmb_multinomPois", starts, method, estimateList, par_inds,
                    tmb_inputs, data, ...)
  }
  
  new("unmarkedFitMPois", fitType = "multinomPois", call = match.call(), 
      formula = formula, formlist = formulas, data = data,
      estimates = fit$estimate_list, sitesRemoved = dm$removed.sites,
      AIC = fit$AIC, opt = fit$opt, negLogLike = fit$opt$value, nllFun = fit$nll, 
      TMB=fit$TMB)
}
