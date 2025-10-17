# The getDesign function creates design matrices and handles missing values
# Many unmarked fitting functions have bespoke getDesign methods
# Some methods do other things like creating indices, etc.

# generic method for unmarkedFrame
# used by distsamp, multinomPois, occu, occuPEN, occuRN, pcount, pcount.spHDS, IDS
setMethod("getDesign", "unmarkedFrame", 
  function(umf, formulas, na.rm = TRUE, ...){

  M <- numSites(umf)
  J <- obsNum(umf)
  y <- getY(umf)

  # Process covariates
  covs <- clean_up_covs(umf)
  
  # Model matrices and offset for state submodel
  X_state <- get_model_matrix(formulas$state, covs$site_covs)
  offset_state <- get_offset(formulas$state, covs$site_covs)
  Z_state <- get_Z(formulas$state, covs$site_covs)

  # Model matrices and offset for detection submodel
  X_det <- get_model_matrix(formulas$det, covs$obs_covs)
  offset_det <- get_offset(formulas$det, covs$obs_covs)
  Z_det <- get_Z(formulas$det, covs$obs_covs)

  # Identify missing values in state covs
  has_na_site <- row_has_na(cbind(X_state, Z_state))
  has_na_site <- matrix(rep(has_na_site, ncol(y)), M)

  # Identify missing values in det covs
  has_na_obs <- row_has_na(cbind(X_det, Z_det))
  has_na_obs <- matrix(has_na_obs, M, J, byrow=TRUE)
  # Multiplying by obsToY handles models where the number of estimated parameters
  # is not the same as the number of observations, such as double-observer models
  has_na_obs <- has_na_obs %*% obsToY(umf) > 0

  # Combine missing value information
  has_na <- has_na_site | has_na_obs
  stopifnot(identical(dim(y), dim(has_na)))
  y[has_na] <- NA
  drop_sites <- row_all_na(y)

  # Drop sites with all NAs if requested
  if(na.rm & any(drop_sites)){
    warning("Site(s) ", paste(which(drop_sites), collapse = ","),
            " dropped due to missing values", call.=FALSE)
    y <- y[!drop_sites,,drop=FALSE]
    X_state <- X_state[!drop_sites,,drop=FALSE]
    offset_state <- offset_state[!drop_sites]
    Z_state <- Z_state[!drop_sites,,drop=FALSE]
    
    drop_sites_obs <- rep(drop_sites, each = J)
    X_det <- X_det[!drop_sites_obs,,drop=FALSE]
    offset_det <- offset_det[!drop_sites_obs]
    Z_det <- Z_det[!drop_sites_obs,,drop=FALSE]
  }

  list(y = y, 
       X_state = X_state, offset_state = offset_state, Z_state = Z_state,
       X_det = X_det, offset_det = offset_det, Z_det = Z_det,
       removed.sites = which(drop_sites))
})


# UnmarkedMultFrame
# used by colext, occuTTD, nmixTTD and base class for G3 and DailMadsen
setMethod("getDesign", "unmarkedMultFrame",
  function(umf, formulas, na.rm = TRUE){

  # Handle varying name of state variable in models that use this method
  state_type <- names(formulas)[1]
  # Process state and detection with generic umf method
  names(formulas)[1] <- "state"
  out <- methods::callNextMethod(umf, formulas = formulas, na.rm = FALSE)
 
  M <- numSites(umf)
  R <- obsNum(umf)
  T <- umf@numPrimary
  J <- R / T

  y <- out$y

  # Process covariates
  # Note drop_final = TRUE to remove unused factor levels
  # The rows are actually dropped in colext()
  covs <- clean_up_covs(umf, drop_final = TRUE)

  # Model matrix for colonization
  X_col <- get_model_matrix(formulas$col, covs$yearly_site_covs)
  offset_col <- get_offset(formulas$col, covs$yearly_site_covs)

  # Model matrix for extinction
  X_ext <- get_model_matrix(formulas$ext, covs$yearly_site_covs)
  offset_ext <- get_offset(formulas$ext, covs$yearly_site_covs)

  # Error if any offsets specified
  if(any(c(offset_col, offset_ext, out$X.offset, out$V.offset) != 0)){
    stop("Offsets not allowed in colext", call.=FALSE)
  }

  # Check missing values in transition parameters
  has_na <- row_has_na(cbind(X_col, X_ext))
  has_na <- rep(has_na, each = J) # expand to match observation dims
  has_na <- matrix(has_na, M, R, byrow=TRUE)
  # We don't care about any missing values during the last period
  # since they are not used in the model. Force has_na to be FALSE for those.
  has_na[1:M, ((T-1)*J+1):R] <- FALSE
  stopifnot(identical(dim(y), dim(has_na)))
  y[has_na] <- NA
  drop_sites <- row_all_na(y)

  # Remove missing sites if requested
  if(na.rm & any(drop_sites)){
    warning("Site(s) ", paste(which(drop_sites), collapse = ","),
            " dropped due to missing values", call.=FALSE)
    y <- y[!drop_sites,,drop=FALSE]
    out$X_state <- out$X_state[!drop_sites,,drop=FALSE]
    drop_sites_per <- rep(drop_sites, each = T)
    X_col <- X_col[!drop_sites_per,,drop=FALSE]
    X_ext <- X_ext[!drop_sites_per,,drop=FALSE]
    drop_sites_obs <- rep(drop_sites, each = R)
    out$X_det <- out$X_det[!drop_sites_obs,,drop=FALSE]
  }

  # Combine outputs
  out <- list(y = y, X_state = out$X_state, X_col = X_col, X_ext = X_ext, X_det = out$X_det,
              removed.sites = which(drop_sites))
  # Rename design matrix for state variable to correct type
  names(out)[names(out) == "X_state"] <- paste0("X_", state_type)
  out
})


# unmarkedFrameG3 (gpcount, gmultmix, gdistsamp, goccu)
setMethod("getDesign", "unmarkedFrameG3",
  function(umf, formulas, na.rm = TRUE){

  # Handle varying name of state variable in models that use this method
  state_type <- names(formulas)[1]
  
  # Process state and detection with generic umf method  
  # Have to use getMethod because this inherits from unmarkedMultFrame
  names(formulas)[1] <- "state"
  getDesign_generic <- methods::getMethod("getDesign", "unmarkedFrame")
  out <- getDesign_generic(umf, formulas = formulas, na.rm = FALSE)
 
  M <- numSites(umf)
  R <- obsNum(umf)
  T <- umf@numPrimary
  J <- numY(umf) / T

  y <- out$y

  # Process covariates
  covs <- clean_up_covs(umf)

  # Model matrix for availability
  X_phi <- get_model_matrix(formulas$phi, covs$yearly_site_covs)
  offset_phi <- get_offset(formulas$phi, covs$yearly_site_covs)

  # Check missing values in availability
  has_na <- row_has_na(X_phi)
  has_na <- rep(has_na, each = J) # expand to match observation dims
  has_na <- matrix(has_na, M, numY(umf), byrow=TRUE)
  stopifnot(identical(dim(y), dim(has_na)))
  y[has_na] <- NA
  drop_sites <- row_all_na(y)

  # Remove missing sites if requested
  if(na.rm & any(drop_sites)){
    warning("Site(s) ", paste(which(drop_sites), collapse = ","),
            " dropped due to missing values", call.=FALSE)
    y <- y[!drop_sites,,drop=FALSE]
    out$X_state <- out$X_state[!drop_sites,,drop=FALSE]
    out$offset_state <- out$offset_state[!drop_sites]
    drop_sites_per <- rep(drop_sites, each = T)
    X_phi <- X_phi[!drop_sites_per,,drop=FALSE]
    offset_phi <- offset_phi[!drop_sites_per]
    drop_sites_obs <- rep(drop_sites, each = R)
    out$X_det <- out$X_det[!drop_sites_obs,,drop=FALSE]
    out$offset_det <- out$offset_det[!drop_sites_obs]
  }

  # Combine outputs
  out <- list(y = y, X_state = out$X_state, offset_state = out$offset_state,
              X_phi = X_phi, offset_phi = offset_phi,
              X_det = out$X_det, offset_det = out$offset_det,
              removed.sites = which(drop_sites))
  # Rename design matrix and offset for state variable to correct type
  names(out)[names(out) == "X_state"] <- paste0("X_", state_type)
  names(out)[names(out) == "offset_state"] <- paste0("offset_", state_type)
  out
})


# unmarkedFrameDailMadsen (pcountOpen, multmixOpen, distsampOpen)
setMethod("getDesign", "unmarkedFrameDailMadsen",
  function(umf, formulas, na.rm = TRUE){

  M <- numSites(umf)
  T <- umf@numPrimary
  R <- obsNum(umf)
  nY <- numY(umf)
  J <- nY / T

  y <- getY(umf)
  delta <- umf@primaryPeriod

  covs <- clean_up_covs(umf, drop_final = TRUE)

  # Model matrix for abundance
  X_lambda <- get_model_matrix(formulas$lambda, covs$site_covs)
  offset_lambda <- get_offset(formulas$lambda, covs$site_covs)

  # Transition probs
  # Drop last period of transition prob design matrices
  # NOTE: this is done outside getDesign for colext
  drop_periods <- rep(1:T, M) == T
  ysc_drop <- covs$yearly_site_covs[!drop_periods,,drop=FALSE]
  X_gamma <- get_model_matrix(formulas$gamma, ysc_drop)
  offset_gamma <- get_offset(formulas$gamma, ysc_drop)

  X_omega <- get_model_matrix(formulas$omega, ysc_drop)
  offset_omega <- get_offset(formulas$omega, ysc_drop)

  X_iota <- get_model_matrix(formulas$iota, ysc_drop)
  offset_iota <- get_offset(formulas$iota, ysc_drop)

  # Detection, which differs by model type
  det_covs <- covs$obs_covs
  if(inherits(umf, "unmarkedFrameDSO")){
    # Distance sampling detection model uses yearly site covs
    # We don't want to drop the last period in this case
    covs_dso <- clean_up_covs(umf, drop_final = FALSE)
    det_covs <- covs_dso$yearly_site_covs
  }
  X_det <-  get_model_matrix(formulas$det, det_covs)
  offset_det <- get_offset(formulas$det, det_covs)

  # Identify missing values in lambda covs
  has_na_site <- row_has_na(X_lambda)
  has_na_site <- matrix(rep(has_na_site, ncol(y)), M)

  # Identify missing values in transition covs
  has_na_trans <- row_has_na(cbind(X_gamma, X_omega, X_iota))
  has_na_trans <- rep(has_na_trans, each = J) # expand to match observation dims
  has_na_trans <- matrix(has_na_trans, M, J*(T-1), byrow=TRUE)
  # Assumption is no missing values for first period, so we add a matrix of FALSE for first period
  # this seems weird but is how the old code did it
  # the effect is that a missing value for a yearly site cov in the first period actually
  # results in missing y for period 2
  has_na_trans <- cbind(matrix(FALSE, M, J), has_na_trans)

  # Identify missing values in det covs
  has_na_obs <- row_has_na(X_det)
  if(inherits(umf, "unmarkedFrameDSO")){
    # Since DSO detection model matrix is based on yearly site covs, we
    # need to expand this to match the dimensions of y
    has_na_obs <- rep(has_na_obs, each = R / T)
  }
  has_na_obs <- matrix(has_na_obs, M, R, byrow=TRUE)
  has_na_obs <- has_na_obs %*% umf@obsToY > 0

  # Combine missing value information
  has_na <- has_na_site | has_na_trans | has_na_obs
  stopifnot(identical(dim(y), dim(has_na)))
  y[has_na] <- NA

  # For DSO, if there are any NAs in a period, the whole period becomes NA
  # NOTE: It's not certain this is actually necessary
  if(inherits(umf, "unmarkedFrameDSO")){
    ymat <- array(y, c(M,J,T))
    obs_any_na <- apply(ymat, c(1,3), function(x) any(is.na(x)))
    any_na_ind <- which(obs_any_na, arr.ind=TRUE)
    if(sum(obs_any_na)>0){
      for(i in 1:nrow(any_na_ind)){
        ymat[any_na_ind[i,1], ,any_na_ind[i,2]] <- NA
      }
    }
    y <- matrix(as.vector(ymat), M, T*J)
  }

  drop_sites <- row_all_na(y)

  # Remove missing sites if requested
  if(na.rm & any(drop_sites)){
    warning("Site(s) ", paste(which(drop_sites), collapse = ","),
            " dropped due to missing values", call.=FALSE)
    y <- y[!drop_sites,,drop=FALSE]
    ya <- array(y, c(nrow(y), J, T))
    yna <- apply(is.na(ya), c(1,3), all)

    delta <- delta[!drop_sites,,drop=FALSE]
    delta <- formatDelta(delta, yna)
 
    X_lambda <- X_lambda[!drop_sites,,drop=FALSE]
    offset_lambda <- offset_lambda[!drop_sites]
    
    drop_sites_per <- rep(drop_sites, each = T-1)
    X_gamma <- X_gamma[!drop_sites_per,,drop=FALSE]
    offset_gamma <- offset_gamma[!drop_sites_per]
    X_omega <- X_omega[!drop_sites_per,,drop=FALSE]
    offset_omega <- offset_omega[!drop_sites_per]
    X_iota <- X_iota[!drop_sites_per,,drop=FALSE]
    offset_iota <- offset_iota[!drop_sites_per]
   
    if(inherits(umf, "unmarkedFrameDSO")){
      # DSO uses yearly site covs here, so dimensions are different
      drop_sites_obs <- rep(drop_sites, each = T)
    } else {
      drop_sites_obs <- rep(drop_sites, each = nY)
    }
    X_det <- X_det[!drop_sites_obs,,drop=FALSE]
    offset_det <- offset_det[!drop_sites_obs]
  } else {
    # Still need to format delta
    ya <- array(y, c(M, J, T))
    yna <- apply(is.na(ya), c(1,3), all)
    delta <- formatDelta(delta, yna)
  }

  # determine if gamma, omega, and iota are scalar, vector, or matrix valued
  # Runtime is much faster for scalars and vectors
  Xgo <- cbind(X_gamma, X_omega, X_iota)
  getGOdims <- function(x) {
    xm <- matrix(x, M, T-1, byrow=TRUE)
    col.table <- apply(xm, 2, table)
    row.table <- apply(xm, 1, table)
    if(is.vector(col.table) & !is.list(col.table)) {
      return("rowvec")
    } else if(is.vector(row.table) & !is.list(row.table)) {
      return("colvec")
    } else
      return("matrix")
  }
  if(length(all.vars(formulas$gamma)) == 0 & length(all.vars(formulas$omega)) == 0 & 
     length(all.vars(formulas$iota)) == 0){
    go.dims <- "scalar"
  } else {
    go.dims.vec <- apply(Xgo, 2, getGOdims)
    if(all(go.dims.vec == "rowvec")){
      go.dims <- "rowvec"
    } else if(all(go.dims.vec == "colvec")){
      ## NOTE: Temporary fix to the problem reported with
      ## time-only-varying covariates
      go.dims <- "matrix" ##"colvec"
    } else {
      go.dims <- "matrix"
    }
  }

	list(y = y, X_lambda = X_lambda, offset_lambda = offset_lambda,
       X_gamma = X_gamma, offset_gamma = offset_gamma,
       X_omega = X_omega, offset_omega = offset_omega,
       X_iota = X_iota, offset_iota = offset_iota,
       X_det = X_det, offset_det = offset_det, 
       delta = delta, removed.sites = which(drop_sites),
       go.dims = go.dims)
})

# Calculate time intervals acknowledging gaps due to NAs
# The first column indicates is time since first primary period + 1
formatDelta <- function(d, yna)
{
    M <- nrow(yna)
    T <- ncol(yna)
    d <- d - min(d, na.rm=TRUE) + 1
    dout <- matrix(NA, M, T)
    dout[,1] <- d[,1]
    dout[,2:T] <- t(apply(d, 1, diff))
    for(i in 1:M) {
        if(any(yna[i,]) & !all(yna[i,])) { # 2nd test for simulate
            last <- max(which(!yna[i,]))
            y.in <- yna[i, 1:last]
            d.in <- d[i, 1:last]
            if(any(y.in)) {
                for(j in last:2) { # first will always be time since 1
                    nextReal <- which(!yna[i, 1:(j-1)])
                    if(length(nextReal) > 0)
                        dout[i, j] <- d[i, j] - d[i, max(nextReal)]
                    else
                        dout[i, j] <- d[i, j] - 1
                    }
                }
            }
        }
    return(dout)
}


# Methods for specific function types------------------------------------------

# gdistremoval
setMethod("getDesign", "unmarkedFrameGDR",
  function(umf, formulas, na.rm=TRUE){

  M <- numSites(umf)
  T <- umf@numPrimary
  Rdist <- ncol(umf@yDistance)
  Jdist <- Rdist/T
  Rrem <- ncol(umf@yRemoval)
  Jrem <- Rrem/T
  yRem <- as.vector(t(umf@yRemoval))
  yDist <- as.vector(t(umf@yDistance))

  sc <- siteCovs(umf)
  oc <- obsCovs(umf)
  ysc <- yearlySiteCovs(umf)

  if(is.null(sc)) sc <- data.frame(.dummy=rep(0, M))
  if(is.null(ysc)) ysc <- data.frame(.dummy=rep(0, M*T))
  if(is.null(oc)) oc <- data.frame(.dummy=rep(0, M*Rrem))

  ysc <- cbind(ysc, sc[rep(1:M, each=T),,drop=FALSE])
  oc <- cbind(oc, ysc[rep(1:nrow(ysc), each=Jrem),,drop=FALSE])

  lam_fixed <- reformulas::nobars(formulas$lambda)
  X_lambda <- model.matrix(lam_fixed,
            model.frame(lam_fixed, sc, na.action=NULL))

  phi_fixed <- reformulas::nobars(formulas$phi)
  X_phi<- model.matrix(phi_fixed,
            model.frame(phi_fixed, ysc, na.action=NULL))

  dist_fixed <- reformulas::nobars(formulas$dist)
  X_dist <- model.matrix(dist_fixed,
            model.frame(dist_fixed, ysc, na.action=NULL))

  rem_fixed <- reformulas::nobars(formulas$rem)
  X_rem <- model.matrix(rem_fixed,
            model.frame(rem_fixed, oc, na.action=NULL))

  Z_lambda <- get_Z(formulas$lambda, sc)
  Z_phi <- get_Z(formulas$phi, ysc)
  Z_dist <- get_Z(formulas$dist, ysc)
  Z_rem <- get_Z(formulas$rem, oc)
 
  # Check if there are missing yearlySiteCovs
  ydist_mat <- apply(matrix(yDist, nrow=M*T, byrow=TRUE), 1, function(x) any(is.na(x)))
  yrem_mat <- apply(matrix(yRem, nrow=M*T, byrow=TRUE), 1, function(x) any(is.na(x)))
  ok_missing_phi_covs <- ydist_mat | yrem_mat
  missing_phi_covs <- apply(X_phi, 1, function(x) any(is.na(x)))  
  if(!all(which(missing_phi_covs) %in% which(ok_missing_phi_covs))){
    stop("Missing yearlySiteCovs values for some observations that are not missing", call.=FALSE)
  }

  # Check if there are missing dist covs
  missing_dist_covs <- apply(X_dist, 1, function(x) any(is.na(x)))
  ok_missing_dist_covs <- ydist_mat
  if(!all(which(missing_dist_covs) %in% which(ok_missing_dist_covs))){
    stop("Missing yearlySiteCovs values for some distance observations that are not missing", call.=FALSE)
  }

  # Check if there are missing rem covs
  missing_obs_covs <- apply(X_rem, 1, function(x) any(is.na(x)))
  missing_obs_covs <- apply(matrix(missing_obs_covs, nrow=M*T, byrow=TRUE), 1, function(x) any(x))
  ok_missing_obs_covs <- yrem_mat
  if(!all(which(missing_obs_covs) %in% which(ok_missing_obs_covs))){
    stop("Missing obsCovs values for some removal observations that are not missing", call.=FALSE)
  }
    
  if(any(is.na(X_lambda))){
    stop("gdistremoval does not currently handle missing values in siteCovs", call.=FALSE)
  }

  list(yDist=yDist, yRem=yRem, X_lambda=X_lambda, X_phi=X_phi, X_dist=X_dist, X_rem=X_rem,
       Z_lambda=Z_lambda, Z_phi=Z_phi, Z_dist=Z_dist, Z_rem=Z_rem)
})


# occuFP
# there are 3 observation formula which are stored in V (true positive detections),
# U (false positive detections), and W (b or probability detetion is certain)
setMethod("getDesign", "unmarkedFrameOccuFP",
  function(umf, formulas, na.rm = TRUE){

  # Process state and true detection with generic umf method
  out <- methods::callNextMethod(umf, formulas = formulas, na.rm = FALSE)
 
  M <- numSites(umf)
  J <- obsNum(umf)
  y <- out$y
  
  # Process covariates
  covs <- clean_up_covs(umf)
  
  # Model matrix and offset for false positives
  X_fp <- get_model_matrix(formulas$fp, covs$obs_covs)
  offset_fp <- get_offset(formulas$fp, covs$obs_covs)

  # Model matrices and offset for b (prob detection is certain)
  X_b <- get_model_matrix(formulas$b, covs$obs_covs)
  offset_b <- get_offset(formulas$b, covs$obs_covs)

  # Check missing values in FP and b
  has_na <- row_has_na(cbind(X_fp, X_b))
  has_na <- matrix(has_na, M, J, byrow=TRUE)
  has_na <- has_na %*% umf@obsToY > 0
  stopifnot(identical(dim(y), dim(has_na)))
  y[has_na] <- NA
  drop_sites <- row_all_na(y)

  # Remove missing sites if requested
  if(na.rm & any(drop_sites)){
    warning("Site(s) ", paste(which(drop_sites), collapse = ","),
            " dropped due to missing values", call.=FALSE)
    y <- y[!drop_sites,,drop=FALSE]

    out$X_state <- out$X_state[!drop_sites,,drop=FALSE]
    out$offset_state <- out$offset_state[!drop_sites]
    drop_sites_obs <- rep(drop_sites, each = J)
    out$X_det <- out$X_det[!drop_sites_obs,,drop=FALSE]
    out$offset_det <- out$offset_det[!drop_sites_obs]

    X_fp <- X_fp[!drop_sites_obs,,drop=FALSE]
    offset_fp <- offset_fp[!drop_sites_obs]
    X_b <- X_b[!drop_sites_obs,,drop=FALSE]
    offset_b <- offset_b[!drop_sites_obs]
  }

  # Combine outputs (U=fp, W=b)
  new_out <- list(X_fp = X_fp, offset_fp = offset_fp, X_b = X_b, offset_b = offset_b)
  out <- c(out, new_out)
  out$y <- y
  out$removed.sites <- which(drop_sites)
  out
})


# occuMS
setMethod("getDesign", "unmarkedFrameOccuMS",
    function(umf, formulas, prm, na.rm=TRUE, newdata=NULL, type="psi")
{

  N <- numSites(umf)
  S <- umf@numStates
  T <- umf@numPrimary
  R <- obsNum(umf)
  J <- R / T
  npsi <- S-1 #Number of free psi values
  nphi <- S^2 - S #Number of free phi values
  np <- S * (S-1) / 2 #Number of free p values

  psiformulas <- formulas$state
  phiformulas <- formulas$transition
  detformulas <- formulas$det
  if(length(psiformulas) != npsi){
    stop(paste(npsi,'formulas are required in psiformulas vector'))
  }

  if(is.null(phiformulas)){
    if(prm == 'condbinom') {
      phiformulas <- rep('~1',6)
    } else {
      phiformulas <- rep('~1',S^2-S)
    }
  } else if(T>1){
    if(length(phiformulas)!=nphi){
      stop(paste(nphi,'formulas are required in phiformulas vector. See data@phiOrder for help'))
    }
  }

  if(length(detformulas) != np){
    stop(paste(np,'formulas are required in detformulas vector'))
  }

  #get placeholder for empty covs if necessary
  get_covs <- function(covs, length){
    if(is.null(covs)){
      return(data.frame(placeHolder = rep(1, length)))
    }
    return(covs)
  }

  #Function to create list of design matrices from list of formulas
  get_dm <- function(formulas, covs, namevec, old_covs=NULL){

    ref_covs <- covs
    if(!is.null(old_covs)) ref_covs <- old_covs

    fac_col <- ref_covs[, sapply(ref_covs, is.factor), drop=FALSE]
    xlevs_all <- lapply(fac_col, levels)

    apply_func <- function(i){
      mf <- model.frame(as.formula(formulas[i]), ref_covs)
      xlevs <- xlevs_all[names(xlevs_all) %in% names(mf)]
      out <- model.matrix(as.formula(formulas[i]),
              model.frame(stats::terms(mf), covs,
                          na.action=stats::na.pass, xlev=xlevs))
      #improve these names
      colnames(out) <- paste(namevec[i], colnames(out))
      out
    }

    out <- lapply(seq_along(formulas), apply_func)
    names(out) <- namevec
    out
  }

  #Generate informative names for p
  get_p_names <- function(S, prm){
    if(prm=='condbinom'){
      return(c('p[1]','p[2]','delta'))
    }
    inds <- matrix(NA,nrow=S,ncol=S)
    inds <- lower.tri(inds,diag=TRUE)
    inds[,1] <- FALSE
    inds <- which(inds,arr.ind=TRUE) - 1
    paste0('p[',inds[,2],inds[,1],']')
  }

  #Informative names for phi
  get_phi_names <- function(np, prm){
    if(prm=='condbinom'){
      return(c(paste0('phi[',0:(S-1),']'),paste0('R[',0:(S-1),']')))
    }
    vals <- paste0('phi[',rep(0:(S-1),each=S),rep(0:(S-1),S),']')
    vals <- matrix(vals,nrow=S)
    diag(vals) <- NA
    c(na.omit(as.vector(vals)))
  }

  #Informative names for psi
  get_psi_names <- function(np, prm){
    if(prm=='condbinom'){
      return(c('psi','R'))
    }
    paste0('psi[',1:np,']')
  }

  #Get vector of parameter count indices from a design matrix list
  get_param_inds <- function(dm_list, offset=0){
    apply_func <- function(i){
      rep(i, ncol(dm_list[[i]]))
    }
    ind_vec <- unlist(lapply(seq_along(dm_list), apply_func))
    start_ind <- c(1,1+which(diff(ind_vec)!=0)) + offset
    stop_ind <- c(start_ind[2:length(start_ind)]-1,length(ind_vec)+offset)

    cbind(start_ind, stop_ind)
  }

  #Get param names from dm_list
  get_param_names <- function(dm_list){
    unlist(lapply(dm_list,colnames))
  }

  #Get observations with NAs across design matrices
  get_na_inds <- function(formulas, covs){
    dm_list <- get_dm(formulas, covs, rep("", length(formulas)))
    dm_mat <- Reduce(cbind, dm_list)
    which(apply(dm_mat, 1, function(x) any(is.na(x))))
  }

  site_covs <- get_covs(siteCovs(umf), N)

  y_site_covs <- get_covs(yearlySiteCovs(umf), N*T)

  #Add site covs to yearly site covs
  if(!is.null(umf@siteCovs)){
    y_site_covs <- cbind(y_site_covs,
                         site_covs[rep(1:N, each=T),,drop=FALSE])
  }

  obs_covs <- get_covs(obsCovs(umf), N*R)
  #Add yearly site covs to obs covs
  if(!is.null(umf@siteCovs) | !is.null(umf@yearlySiteCovs)){
    obs_covs <- cbind(obs_covs,
                      y_site_covs[rep(1:(N*T), each=J),,drop=FALSE])
  }

  ## in order to drop factor levels that only appear in last year,
  ## replace last year with NAs and use drop=TRUE
  y_site_covs[seq(T,N*T,by=T),] <- NA
  y_site_covs <- as.data.frame(lapply(y_site_covs, function(x) x[,drop = TRUE]))
  #Actually just remove last year
  y_site_covs <- y_site_covs[-seq(T,N*T,by=T),,drop=FALSE]

  y <- getY(umf)

  #Handle NAs
  removed.sites <- NA
  if(na.rm){

    #Det
    ylong <- as.vector(t(y))
    miss_det_cov <- get_na_inds(detformulas, obs_covs)
    miss_y <- which(is.na(ylong))
    new_na <- miss_det_cov[!miss_det_cov%in%miss_y]
    if(length(new_na)>0){
      warning('Some observations removed because covariates were missing')
      ylong[new_na] <- NA
      y <- matrix(ylong,nrow=N,ncol=R,byrow=T)
    }

    #State
    check_site_na <- function(yrow){
      if(T==1) return(all(is.na(yrow)))
      y_mat <- matrix(yrow, nrow=J)
      pp_na <- apply(y_mat,2,function(x) all(is.na(x)))
      if(any(pp_na)){
        return(TRUE)
      }
      return(FALSE)
    }

    all_y_na <- which(apply(y,1, check_site_na ))
    if(length(all_y_na)>0){
      warning("Some sites removed because all y values in a primary period were missing")
    }
    miss_covs <- get_na_inds(psiformulas, site_covs)
    if(length(miss_covs)>0){
      warning("Some sites removed because site covariates were missing")
    }
    removed.sites <- sort(unique(c(all_y_na,miss_covs)))

    if(T>1){
      ysc_na <- get_na_inds(phiformulas, y_site_covs)
      if(length(ysc_na) > 0){
        stop("Some sites are missing yearly site covs")
      }
    }

    if(length(removed.sites)>0){
      ymap <- as.vector(t(matrix(rep(1:N,each=R),ncol=R,byrow=T)))
      site_covs <- site_covs[-removed.sites,,drop=FALSE]
      obs_covs <- obs_covs[!ymap%in%removed.sites,,drop=FALSE]
      if(T>1){
        ysc_map <- as.vector(t(matrix(rep(1:N,each=(T-1)),ncol=(T-1),byrow=T)))
        y_site_covs <- y_site_covs[!ysc_map%in%removed.sites,,drop=FALSE]
      }
      y <- y[-removed.sites,,drop=FALSE]
      N <- nrow(y)
    }

  }

  site_ref <- site_covs
  ysc_ref <- y_site_covs
  obs_ref <- obs_covs

  #Assign newdata as the covariate frame if it is provided
  if(!is.null(newdata)){
    if(type == "psi"){
      site_covs <- newdata
    } else if(type == "phi"){
      y_site_covs <- newdata
    } else if(type == "det"){
      obs_covs <- newdata
    }
  }

  dm_state <- get_dm(psiformulas, site_covs,
                     get_psi_names(length(psiformulas),prm), site_ref)
  nSP <- length(get_param_names(dm_state))
  state_ind <- get_param_inds(dm_state) #generate ind matrix in function

  nPP <- 0; dm_phi <- list(); phi_ind <- c()
  if(T>1){
    dm_phi <- get_dm(phiformulas, y_site_covs,
                   get_phi_names(length(phiformulas),prm), ysc_ref)
    nPP <- length(get_param_names(dm_phi))
    phi_ind <- get_param_inds(dm_phi, offset=nSP)
  }

  dm_det <- get_dm(detformulas, obs_covs, get_p_names(S,prm), obs_ref)
  det_ind <- get_param_inds(dm_det, offset=(nSP+nPP))

  param_names <- c(get_param_names(dm_state),
                   get_param_names(dm_phi),
                   get_param_names(dm_det))

  mget(c("y","dm_state","state_ind","nSP",
         "dm_phi","phi_ind","nPP",
         "dm_det","det_ind","param_names","removed.sites"))

})


# occuMulti
setMethod("getDesign", "unmarkedFrameOccuMulti",
    function(umf, formulas, maxOrder, na.rm=TRUE, warn=FALSE,
             newdata=NULL, type="state")
{

  #Format formulas
  stateformulas <- formulas$state
  detformulas <- formulas$det
  #Workaround for parameters fixed at 0
  fixed0 <- stateformulas %in% c("~0","0")
  stateformulas[fixed0] <- "~1"
  stateformulas <- lapply(stateformulas,as.formula)
  detformulas <- lapply(detformulas,as.formula)

  #Generate some indices
  S <- length(umf@ylist) # of species
  if(missing(maxOrder)){
    maxOrder <- S
  }
  z <- expand.grid(rep(list(1:0),S))[,S:1] # z matrix
  colnames(z) <- names(umf@ylist)
  M <- nrow(z) # of possible z states

  # f design matrix
  if(maxOrder == 1){
    dmF <- as.matrix(z)
  } else {
    dmF <- model.matrix(as.formula(paste0("~.^",maxOrder,"-1")),z)
  }
  nF <- ncol(dmF) # of f parameters

  J <- ncol(umf@ylist[[1]]) # max # of samples at a site
  N <- nrow(umf@ylist[[1]]) # of sites

  #Check formulas
  if(length(stateformulas) != nF)
    stop(paste(nF,"formulas are required in stateformulas list"))
  if(length(detformulas) != S)
    stop(paste(S,"formulas are required in detformulas list"))

  if(is.null(siteCovs(umf))) {
    site_covs <- data.frame(placeHolderSite = rep(1, N))
  } else {
    site_covs <- siteCovs(umf)
  }

  if(is.null(obsCovs(umf))) {
    obs_covs <- data.frame(placeHolderObs = rep(1, J*N))
  } else {
    obs_covs <- obsCovs(umf)
  }

  #Add site covs to obs covs if we aren't predicting with newdata
  # Record future column names for obsCovs
  col_names <- c(colnames(obs_covs), colnames(site_covs))

  # add site covariates at observation-level
  obs_covs <- cbind(obs_covs, site_covs[rep(1:N, each = J),])
  colnames(obs_covs) <- col_names

  #Re-format ylist
  index <- 1
  ylong <- lapply(umf@ylist, function(x) {
                   colnames(x) <- 1:J
                   x <- cbind(x,site=1:N,species=index)
                   index <<- index+1
                   x
          })
  ylong <- as.data.frame(do.call(rbind,ylong))
  ylong <- reshape(ylong, idvar=c("site", "species"), varying=list(1:J),
                   v.names="value", direction="long")
  ylong <- reshape(ylong, idvar=c("site","time"), v.names="value",
                    timevar="species", direction="wide")
  ylong <- ylong[order(ylong$site, ylong$time), ]

  #Remove missing values
  if(na.rm){
    naSiteCovs <- which(apply(site_covs, 1, function(x) any(is.na(x))))
    if(length(naSiteCovs>0)){
      stop(paste("Missing site covariates at sites:",
                 paste(naSiteCovs,collapse=", ")))
    }

    naY <- apply(ylong, 1, function(x) any(is.na(x)))
    naCov <- apply(obs_covs, 1, function(x) any(is.na(x)))
    navec <- naY | naCov

    sites_with_missingY <- unique(ylong$site[naY])
    sites_with_missingCov <- unique(ylong$site[naCov])

    ylong <- ylong[!navec,,drop=FALSE]
    obs_covs <- obs_covs[!navec,,drop=FALSE]

    no_data_sites <- which(! 1:N %in% ylong$site)
    if(length(no_data_sites>0)){
      stop(paste("All detections and/or detection covariates are missing at sites:",
                  paste(no_data_sites,collapse=", ")))
    }

    if(sum(naY)>0&warn){
      warning(paste("Missing detections at sites:",
                    paste(sites_with_missingY,collapse=", ")))
    }
    if(sum(naCov)>0&warn){
      warning(paste("Missing detection covariate values at sites:",
                    paste(sites_with_missingCov,collapse=", ")))
    }

  }

  #Start-stop indices for sites
  yStart <- c(1,1+which(diff(ylong$site)!=0))
  yStop <- c(yStart[2:length(yStart)]-1,nrow(ylong))

  y <- as.matrix(ylong[,3:ncol(ylong)])

  #Indicator matrix for no detections at a site
  Iy0 <- do.call(cbind, lapply(umf@ylist,
                               function(x) as.numeric(rowSums(x, na.rm=T)==0)))

  #Save formatted covariate frames for use in model frames
  #For predicting with formulas etc
  site_ref <- site_covs
  obs_ref <- obs_covs

  #Assign newdata as the covariate frame if it is provided
  if(!is.null(newdata)){
    if(type == "state"){
      site_covs <- newdata
    } else if(type == "det"){
      obs_covs <- newdata
    }
  }

  #Design matrices + parameter counts
  #For f/occupancy
  fInd <- c()
  sf_no0 <- stateformulas[!fixed0]
  var_names <- colnames(dmF)[!fixed0]
  dmOcc <- lapply(seq_along(sf_no0),function(i){
                    fac_col <- site_ref[, sapply(site_ref, is.factor), drop=FALSE]
                    mf <- model.frame(sf_no0[[i]], site_ref)
                    xlevs <- lapply(fac_col, levels)
                    xlevs <- xlevs[names(xlevs) %in% names(mf)]
                    out <- model.matrix(sf_no0[[i]],
                                        model.frame(stats::terms(mf), site_covs, na.action=stats::na.pass, xlev=xlevs))
                    colnames(out) <- paste('[',var_names[i],'] ',
                                           colnames(out), sep='')
                    fInd <<- c(fInd,rep(i,ncol(out)))
                    out
          })
  fStart <- c(1,1+which(diff(fInd)!=0))
  fStop <- c(fStart[2:length(fStart)]-1,length(fInd))
  occParams <- unlist(lapply(dmOcc,colnames))
  nOP <- length(occParams)

  #For detection
  dInd <- c()
  dmDet <- lapply(seq_along(detformulas),function(i){
                    fac_col <- obs_ref[, sapply(obs_ref, is.factor), drop=FALSE]
                    mf <- model.frame(detformulas[[i]], obs_ref)
                    xlevs <- lapply(fac_col, levels)
                    xlevs <- xlevs[names(xlevs) %in% names(mf)]
                    out <- model.matrix(detformulas[[i]],
                                        model.frame(stats::terms(mf), obs_covs, na.action=stats::na.pass, xlev=xlevs))
                    colnames(out) <- paste('[',names(umf@ylist)[i],'] ',
                                           colnames(out),sep='')
                    dInd <<- c(dInd,rep(i,ncol(out)))
                    out
          })
  dStart <- c(1,1+which(diff(dInd)!=0)) + nOP
  dStop <- c(dStart[2:length(dStart)]-1,length(dInd)+nOP)
  detParams <- unlist(lapply(dmDet,colnames))
  #nD <- length(detParams)

  #Combined
  paramNames <- c(occParams,detParams)
  nP <- length(paramNames)

  mget(c("N","S","J","M","nF","fStart","fStop","fixed0","dmF","dmOcc","dmDet",
         "dStart","dStop","y","yStart","yStop","Iy0","z","nOP","nP","paramNames"))
})


# Utility functions used in several methods------------------------------------

# Convert NULL data frames to dummy data frames of proper dimension
# Add site covs to yearlysitecovs, ysc to obs covs, etc.
# Drop final year of ysc if necessary
# Add observation number (for backwards compatibility)
clean_up_covs <- function(object, drop_final=FALSE, addObsNum = TRUE){
  M <- numSites(object)
  R <- obsNum(object)
  T <- 1
  J <- R
  is_mult <- methods::.hasSlot(object, "numPrimary")
  if(is_mult){
    T <- object@numPrimary
    J <- R/T
  }

  sc <- siteCovs(object)
  if(is.null(sc)) sc <- data.frame(.dummy=rep(1,M))
  out <- list(site_covs=sc)

  if(is_mult){
    ysc <- yearlySiteCovs(object)
    if(is.null(ysc)) ysc <- data.frame(.dummy2=rep(1,M*T))
    ysc <- cbind(ysc, sc[rep(1:M, each=T),,drop=FALSE])
  }

  if(methods::.hasSlot(object, "obsCovs")){
    oc <- obsCovs(object)
    if(is.null(oc)) oc <- data.frame(.dummy3=rep(1,M*T*J))
    if(is_mult){
      oc <- cbind(oc, ysc[rep(1:(M*T), each=J),,drop=FALSE])
    } else {
      oc <- cbind(oc, sc[rep(1:M, each=J),,drop=FALSE])
    }

    # Add observation number (mainly for backwards compatibility)
    if(addObsNum & !"obsNum" %in% names(oc)){
      oc$obsNum <- as.factor(rep(1:R, M))
    }

    out$obs_covs=oc
  }

  if(is_mult){
    if(drop_final & (T > 1)){
      # Drop final year of data at each site
      # Also drop factor levels only found in last year of data
      ysc <- drop_final_year(ysc, M, T)
    }
    out$yearly_site_covs <- ysc
  }

  out
}

#Remove data in final year of yearlySiteCovs (replacing with NAs)
#then drop factor levels found only in that year
drop_final_year <- function(dat, nsites, nprimary){
  dat[seq(nprimary, nsites*nprimary, by=nprimary), ] <- NA
  dat <- lapply(dat, function(x) x[,drop = TRUE])
  as.data.frame(dat)
}

get_model_matrix <- function(formula, covs){
  miss_vars <- all.vars(formula)[!all.vars(formula) %in% names(covs)]
  if(length(miss_vars) > 0){
    stop(paste("Variable(s)", paste(miss_vars, collapse=", "), 
                "not found in covariates"), call.=FALSE)
  }

  formula <- reformulas::nobars(formula)

  X_mf <- model.frame(formula, covs, na.action = stats::na.pass)
  model.matrix(formula, X_mf)
}

get_offset <- function(formula, covs){
  formula <- reformulas::nobars(formula)
  X_mf <- model.frame(formula, covs, na.action = stats::na.pass)
  offset <- as.vector(model.offset(X_mf))
  if (!is.null(offset)) {
    offset[is.na(offset)] <- 0
  } else {
    offset <- rep(0, nrow(X_mf)) 
  }
  offset
}

row_has_na <- function(mat){
  apply(mat, 1, function(x) any(is.na(x))) 
}

row_all_na <- function(mat){
  apply(mat, 1, function(x) all(is.na(x))) 
}
