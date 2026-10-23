occuCT <- function(formula, data, interval, starts, method = "BFGS",
                   se = TRUE, ...) {

  # Check arguments------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccuCT"))
    stop("Data is not an unmarkedFrameOccuCT object.")

  formulas <- split_formula(formula)
  names(formulas) <- c("det", "state")

  gd <- getDesign(data, formulas, interval = interval)
  # Adjustments for TMB
  gd$y <- NULL
  gd$dep_state <- gd$dep_state - 1
  gd$dep_det <- gd$dep_det - 1
  gd$dep_det[,2] <- gd$dep_det[,2] - gd$dep_det[,1] + 1
  gd$dep_int <- gd$dep_int - 1
  gd$dep_int[,2] <- gd$dep_int[,2] - gd$dep_int[,1] + 1

  # Parameter list
  state_ind <- 1:ncol(gd$X_state)
  det_ind <- 1:ncol(gd$X_det) + max(state_ind)
  nP <- max(det_ind)
  params <- list(beta_state = rep(0, length(state_ind)), 
                 beta_det = rep(0, length(det_ind)))

  # Create TMB object
  tmb_mod <- TMB::MakeADFun(data = c(model = "tmb_occuCT", gd),
                            parameters = params,
                            DLL = "unmarked_TMBExports", silent = TRUE)

  # Optimize TMB object
  if(missing(starts) || is.null(starts)) starts <- tmb_mod$par
  opt <- optim(starts, fn = tmb_mod$fn, gr = tmb_mod$gr, method = method,
               hessian = se, ...)

  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  names(ests) <- c(colnames(gd$X_state), colnames(gd$X_det))
  fmAIC <- 2 * opt$value + 2 * nP

  state_est <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                                estimates = ests[state_ind],
                                covMat = covMat[state_ind, state_ind, drop=FALSE],
                                fixed = 1:ncol(gd$X_state),
                                invlink = "logistic",
                                invlinkGrad = "logistic.grad",
                                randomVarInfo=list()
                               )

  det_est <- unmarkedEstimate(name = "Detection rate", short.name = "lambda",
                              estimates = ests[det_ind],
                              covMat = covMat[det_ind, det_ind, drop=FALSE],
                              fixed = 1:ncol(gd$X_det),
                              invlink = "exp",
                              invlinkGrad = "exp",
                              randomVarInfo=list()
                             )

  estimate_list <- unmarkedEstimateList(list(state=state_est, det=det_est))

  new("unmarkedFitOccuCT", fitType = "occuCT", call = match.call(),
      formlist=formulas, data = data, sitesRemoved = numeric(0),
      estimates = estimate_list, AIC = fmAIC, opt = opt,
      negLogLike = opt$value, nllFun = tmb_mod$fn, TMB=tmb_mod)
}

# Interpolation stuff----------------------------------------------------------

linear_interpolation <- function(times, cov_times, cov_vals){
  
  stopifnot(length(cov_times) == length(cov_vals))

  # Make sure times are sorted
  times <- times[order(times)]

  ord <- order(cov_times)
  cov_times <- cov_times[ord]
  cov_vals <- cov_vals[ord] 
  
  ind <- findInterval(times, cov_times)

  time_s <- as.numeric(cov_times)
  min_time <- min(time_s)
  time_s <- time_s - min_time

  tnum <- as.numeric(times) - min_time

  out <- rep(NA, length(times))

  keepindex <- pmatch(times, cov_times)
  if(any(!is.na(keepindex))){
    insert_idx <- which(!is.na(keepindex))
    out[insert_idx] <- cov_vals[na.omit(keepindex)]
  }
  must_interp <- which(is.na(out))

  if(length(must_interp) > 0){
    interp_vals <- sapply(must_interp, function(i){
      x0 <- time_s[ind[i]]
      x1 <- time_s[ind[i]+1]
      y0 <- cov_vals[ind[i]]
      y1 <- cov_vals[ind[i]+1]
      x <- tnum[i]
      (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
    })
    out[must_interp] <- interp_vals
  }
  out
}

interpolate_covs_for_integration <- function(umf, formula, interval){
  oc <- umf@obsCovs
  if(is.null(oc)) return(NULL)
  # Function to use for interpolation
  fun <- umf@interpolationFunction
  if(is.null(fun)) fun <- linear_interpolation

  deploys <- unique(umf@deployments$deployment)

  if(is.data.frame(oc)){
    # convert to list format
    covs <- names(oc)[3:length(oc)]
    oc <- lapply(3:length(oc), function(i){
      oc[,c(1:2, i)]
    })
    names(oc) <- covs
  } else {
    covs <- names(oc)
  }

  # Only process covariates actually used in formula
  form_covs <- all.vars(formula)
  if(!is.null(umf@siteCovs)){
    # pull out covs used in det formula
    site_covs <- names(umf@siteCovs)
    form_covs <- form_covs[! form_covs %in% site_covs]
  }
  if(!is.null(umf@deploymentCovs)){
    # pull out covs used in det formula
    dep_covs <- names(umf@deploymentCovs)
    form_covs <- form_covs[! form_covs %in% dep_covs]
  }
  stopifnot(all(form_covs %in% covs))
  covs <- covs[covs %in% form_covs]
  oc <- oc[covs]

  # Iterate over deployments
  out <- lapply(deploys, function(this_deploy){
    deploy_info <- umf@deployments[umf@deployments$deployment == this_deploy,]

    # Get sequence of values for integration
    st <- deploy_info$starttime
    en <- deploy_info$endtime
    t <- seq(st, en, by = interval)

    interpolate_covs(this_deploy, oc, t, fun)
  })
  do.call(rbind, out)
}

interpolate_covs_for_y <- function(umf, formula){
  oc <- umf@obsCovs
  if(is.null(oc)) return(NULL)
  # Function to use for interpolation
  fun <- umf@interpolationFunction
  if(is.null(fun)) fun <- linear_interpolation

  deploys <- unique(umf@y$deployment)

  if(is.data.frame(oc)){
    # convert to list format
    covs <- names(oc)[3:length(oc)]
    oc <- lapply(3:length(oc), function(i){
      oc[,c(1:2, i)]
    })
    names(oc) <- covs
  } else {
    covs <- names(oc)
  }

  # Only process covariates actually used in formula
  form_covs <- all.vars(formula)
  if(!is.null(umf@siteCovs)){
    # pull out covs used in det formula
    site_covs <- names(umf@siteCovs)
    form_covs <- form_covs[! form_covs %in% site_covs]
  }
  if(!is.null(umf@deploymentCovs)){
    # pull out covs used in det formula
    dep_covs <- names(umf@deploymentCovs)
    form_covs <- form_covs[! form_covs %in% dep_covs]
  }
  stopifnot(all(form_covs %in% covs))
  covs <- covs[covs %in% form_covs]
  oc <- oc[covs]

  # Iterate over deployments
  out <- lapply(deploys, function(this_deploy){
    # Get sequence of values for integration
    deploy_obs <- umf@y[umf@y$deployment == this_deploy,]
    t <- deploy_obs$obstime

    interpolate_covs(this_deploy, oc, t, fun)
  })
  do.call(rbind, out)
}

# Common code used in both functions above
interpolate_covs <- function(this_deploy, oc, t, fun){
  
  if(length(oc) == 0){
    return(data.frame(deployment = this_deploy, time = t, dummy_cov = 0))
  }

  covs <- names(oc)

  out_cov <- lapply(covs, function(this_cov){
    covsub <- oc[[this_cov]]
    deploy_covs <- covsub[covsub$deployment == this_deploy,]
    
    cov_times <- deploy_covs$covtime
    cov_values <- deploy_covs[[this_cov]]
    stopifnot(length(cov_times) == length(cov_values))
      
    fun(t, cov_times, cov_values)
  })
  names(out_cov) <- covs
  out_cov <- as.data.frame(out_cov)
    
  data.frame(deployment = this_deploy, time = t, out_cov)
}

deployment_inds <- function(covs_df, umf){
  deploys <- umf@deployments$deployment
  out <- lapply(deploys, function(x){
    inds <- which(covs_df$deployment == x)
    if(length(inds) > 0) return(range(inds))
    c(NA, NA)
  })
  out <- as.matrix(do.call(rbind, out))
  rownames(out) <- deploys
  out
}
