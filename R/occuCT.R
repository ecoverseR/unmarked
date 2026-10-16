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
