# General predict method-------------------------------------------------------

# Common predict function for all fit types
# with exception of occuMulti and occuMS (at the end of this file)
# Calls internal fit-type specific function
setMethod("predict", "unmarkedFit",
  function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
           appendData = FALSE, level=0.95, re.form=NULL, ...){
  predict_internal(object, type, newdata, backTransform, na.rm, appendData,
                   level, re.form, ...)
})

# Basic internal predict method for unmarkedFit objects
setMethod("predict_internal", "unmarkedFit",
  function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
           appendData = FALSE, level=0.95, re.form=NULL, ...){
  
  # If no newdata, get actual data
  if(missing(newdata) || is.null(newdata)) newdata <- object@data

  # Check inputs
  check_predict_arguments(object, type, newdata)

  # Get model matrix (X) and offset
  is_raster <- FALSE
  # If newdata is an unmarkedFrame, use getDesign via predict_inputs_from_umf()
  if(inherits(newdata, "unmarkedFrame")){
    # Generate model matrix and offsets
    pred_inps <- predict_inputs_from_umf(object, type, newdata, na.rm, re.form)
  } else {
    # If newdata is provided
    # 1. Get original covariates and appropriate formula for type
    orig_data <- get_covariates(object, type)
    orig_formula <- get_formula(object, type)

    # 2. If newdata is raster, get newdata from raster as data.frame
    if(inherits(newdata, c("RasterLayer","RasterStack","SpatRaster"))){
      is_raster <- TRUE
      orig_raster <- newdata
      check_vars <- all.vars(orig_formula)
      if(!is.null(re.form) && is.na(re.form)) check_vars <- all.vars(reformulas::nobars(orig_formula))
      newdata <- newdata_from_raster(newdata, check_vars)
    }

    # 3. Make model matrix and offset with newdata, informed by original data
    pred_inps <- make_mod_matrix(orig_formula, orig_data, newdata, re.form)
  }

  # Calculate predicted values in chunks (for speed) based on X and offset
  out <- predict_by_chunk(object, type, level, pred_inps$X, pred_inps$offset,
                          chunk_size = 70, backTransform, re.form)

  # Convert output to raster if newdata was raster
  if(is_raster){
    out <- raster_from_predict(out, orig_raster, appendData)
  } else if(appendData){
    # Append data if needed
    out <- data.frame(out, as(newdata, "data.frame"))
  }

  out

})

check_predict_arguments <- function(object, type, newdata){
  # Check if type is supported (i.e., is it in names(object)?)
  opts <- names(object)
  if(!type %in% opts){
    stop("Valid types are ", paste(opts, collapse=", "), call.=FALSE)
  }

  # Check newdata class
  if(!inherits(newdata, c("unmarkedFrame", "data.frame", "RasterLayer", "RasterStack", "SpatRaster"))){
    stop("newdata must be unmarkedFrame, data.frame, RasterLayer, RasterStack, or SpatRaster", call.=FALSE)
  }
  invisible(TRUE)
}

# Nested functions in formulas result in incorrect predictions.
# This seems to be a problem with model.matrix and not something we can fix.
check_nested_formula_functions <- function(formula){
  check <- 
    grepl("I(scale(", safeDeparse(formula), fixed=TRUE) |
    grepl("scale(I(", safeDeparse(formula), fixed=TRUE)
  if(check){
    warning("Predictions based on a model with nested functions in the formula like I(scale(...)^2) are probably incorrect. Scale the covariates manually instead.", call.=FALSE)
  }
  invisible()
}

# Generating inputs from an unmarked frame (e.g. when no newdata) using getDesign
predict_inputs_from_umf <- function(object, type, newdata, na.rm, re.form){

  dm <- getDesign(newdata, object@formlist, na.rm = na.rm)
  select_names <- paste0(c("X_", "Z_", "offset_"), type)
  select <- dm[select_names]
  names(select) <- c("X", "Z", "offset")

  # Scalar parameters like scale, alpha
  if(type %in% names(object) & !type %in% names(object@formlist)){
    return(list(X = matrix(1, nrow=1, ncol=1), offset = NULL))
  }

  # Add random effects matrix if needed
  if(is.null(re.form) & has_random(object@formlist[[type]])){
    stopifnot(!is.null(select$Z))
    select$X <- cbind(select$X, select$Z)
  }

  list(X = select$X, offset = select$offset)
}

# Convert a raster into a data frame to use as newdata
newdata_from_raster <- function(object, vars){
  if(inherits(object, "Raster")){
    if(!requireNamespace("raster", quietly=TRUE)) stop("raster package required", call.=FALSE)
    nd <- raster::as.data.frame(object)
    # Handle factor rasters
    is_fac <- raster::is.factor(object)
    rem_string <- paste(paste0("^",names(object),"_"), collapse="|")
    names(nd)[is_fac] <- gsub(rem_string, "", names(nd)[is_fac])
  } else if(inherits(object, "SpatRaster")){
    if(!requireNamespace("terra", quietly=TRUE)) stop("terra package required", call.=FALSE)
    nd <- terra::as.data.frame(object, na.rm=FALSE)
  }
  # Check if variables are missing
  no_match <- vars[! vars %in% names(nd)]
  if(length(no_match) > 0){
    stop(paste0("Variable(s) ",paste(no_match, collapse=", "), " not found in raster(s)"),
         call.=FALSE)
  }
  return(nd)
}

# Convert predict output into a raster
raster_from_predict <- function(pr, object, appendData){
  if(inherits(object, "Raster")){
    new_rast <- data.frame(raster::coordinates(object), pr)
    new_rast <- raster::stack(raster::rasterFromXYZ(new_rast))
    raster::crs(new_rast) <- raster::crs(object)
    if(appendData) new_rast <- raster::stack(new_rast, object)
  } else if(inherits(object, "SpatRaster")){
    new_rast <- data.frame(terra::crds(object, na.rm=FALSE), pr)
    new_rast <- terra::rast(new_rast, type="xyz")
    terra::crs(new_rast) <- terra::crs(object)
    if(appendData) new_rast <- c(new_rast, object)
  }
  new_rast
}

# Function to make model matrix and offset from formula, newdata and original data
# This function makes sure factor levels in newdata match, and that
# any functions in the formula are handled properly (e.g. scale)
make_mod_matrix <- function(formula, data, newdata, re.form=NULL){
  check_nested_formula_functions(formula)
  form_nobars <- reformulas::nobars(formula)
  mf <- model.frame(form_nobars, data, na.action=stats::na.pass)
  X.terms <- stats::terms(mf)
  fac_cols <- data[, sapply(data, is.factor), drop=FALSE]
  xlevs <- lapply(fac_cols, levels)
  xlevs <- xlevs[names(xlevs) %in% names(mf)]
  nmf <- model.frame(X.terms, newdata, na.action=stats::na.pass, xlev=xlevs)
  #X <- model.matrix(X.terms, newdata, xlev=xlevs)
  X <- model.matrix(form_nobars, nmf)
  offset <- model.offset(nmf)
  if(is.null(re.form) & !is.null(reformulas::findbars(formula))){
    Z <- get_Z(formula, data, newdata)
    X <- cbind(X, Z)
  }
  list(X=X, offset=offset)
}

# Take inputs (most importantly model matrix and offsets) and generate prediction
# done in chunks for speed, 70 was optimal after tests
setMethod("predict_by_chunk", "unmarkedFit",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){

  if(is.vector(xmat)) xmat <- matrix(xmat, nrow=1)
  nr <- nrow(xmat)
  ind <- rep(1:ceiling(nr/chunk_size), each=chunk_size, length.out=nr)

  # should find a way to keep xmat sparse, but it doesn't
  # work with linearComb
  #x_chunk <- lapply(split(as.data.frame(xmat), ind), as.matrix)
  x_chunk <- lapply(unique(ind),
                    function(i) as.matrix(xmat[ind==i,,drop=FALSE]))

  if(is.null(offsets)) offsets <- rep(0, nr)
  off_chunk <- split(offsets, ind)
  out <- mapply(function(x_i, off_i){
    has_na <- apply(x_i, 1, function(x_i) any(is.na(x_i)))
    # Work around linearComb bug where there can't be NAs in inputs
    x_i[has_na,] <- 0
    off_i[has_na] <- 0
    lc <- linearComb(object, x_i, type, offset=off_i, re.form=re.form)
    if(backTransform) lc <- backTransform(lc)
    out <- data.frame(Predicted=coef(lc))
    if(!is.null(level)){
      se <- SE(lc)
      ci <- confint(lc, level=level)
      out$SE <- se
      out$lower <- ci[,1]
      out$upper <- ci[,2]
    }
    out[has_na,] <- NA
    out
    }, x_chunk, off_chunk, SIMPLIFY=FALSE)
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
})

# Special predict approach for ZIP distribution in pcount
setMethod("predict_by_chunk", "unmarkedFitPCount",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){
  if(type == "state" & object@mixture == "ZIP"){

    out <- data.frame(matrix(NA, nrow(xmat), 4))
    names(out) <- c("Predicted", "SE", "lower", "upper")

    psi.hat <- plogis(coef(object, type="psi"))
    lamEst <- object["state"]
    psiEst <- object["psi"]
    fixedOnly <- !is.null(re.form)
    lam.mle <- coef(lamEst, fixedOnly=fixedOnly)
    lam_vcov <- vcov(lamEst, fixedOnly=fixedOnly)
    if(is.null(offsets)) offsets <- rep(0, nrow(xmat))

    for(i in 1:nrow(xmat)) {
      if(nrow(xmat) > 5000) {
          if(i %% 1000 == 0)
              cat("  doing row", i, "of", nrow(xmat), "\n")
      }
      if(any(is.na(xmat[i,]))) next
      ## for the ZIP model the predicted values on the log scale have us
      ## add log(1-psi.hat) to the normal linear prediction
      out$Predicted[i] <-   xmat[i,] %*% lam.mle + offsets[i] + log(1 - psi.hat)
      ## to compute the approximate SE, I compute the variance of the usual
      ## linear part -- that is easy, and to that I add the variance of
      ## log(1-psi.hat) obtained by the delta approximation
      logit.psi<-coef(object,type="psi")
      #  To do that I took derivative of log(1-psi.hat) using application
      #  of chain rule.... hopefully correctly.
      delta.approx.2ndpart<-   ( ((1/(1-psi.hat))*(exp(logit.psi)/((1+exp(logit.psi))^2)))^2 ) * (SE(psiEst)^2)
      ## now the SE is the sqrt of the whole thing
      out$SE[i]<- sqrt( t(xmat[i,])%*% lam_vcov %*%xmat[i,] + delta.approx.2ndpart   )

      #From Mike Meredith
      alf <- (1 - level) / 2
      crit<-qnorm(c(alf, 1 - alf))
      ci <- out$Predicted[i] + crit * out$SE[i]
      out$lower[i]<- ci[1]
      out$upper[i]<- ci[2]
      if(backTransform){
        out$Predicted[i] <- exp(out$Predicted[i])
        ### If back-transform, delta approx says var = (exp(linear.predictor)^2)*Var(linear.predictor)
        ### also I exponentiate the confidence interval.....
        out$SE[i]<- out$Predicted[i]*out$SE[i]
        ci<-exp(ci)
      }
      out$lower[i] <- ci[1]
      out$upper[i] <- ci[2]
    }
    return(out)
  }
  methods::callNextMethod(object, type, level, xmat, offsets, chunk_size,
                          backTransform, re.form, ...)
})

# Special handling for ZIP distribution
setMethod("predict_by_chunk", "unmarkedFitDailMadsen",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){
  if(type == "lambda" & object@mixture == "ZIP"){
    warning("Method to compute SE for ZIP model has not been written", call.=FALSE)
    out <- data.frame(matrix(NA, nrow(xmat), 4))
    names(out) <- c("Predicted", "SE", "lower", "upper")
    lam.mle <- coef(object, type="lambda")
    psi.hat <- plogis(coef(object, type="psi"))
    if(is.null(offsets)) offsets <- rep(0, nrow(xmat))
    out$Predicted <- as.numeric(xmat %*% lam.mle + offsets + log(1 - psi.hat))
    if(backTransform) out$Predicted <- exp(out$Predicted)
    return(out)
  }
  methods::callNextMethod(object, type, level, xmat, offsets, chunk_size,
                          backTransform, re.form, ...)
})

# IDS--------------------------------------------------------------------------

# Uses IDS_convert_class to allow pass-through to prediction using an unmarkedFitDS object
setMethod("predict_internal", "unmarkedFitIDS", function(object, type, newdata,
          backTransform=TRUE, na.rm=FALSE, appendData=FALSE, level=0.95, re.form=NULL, ...){
  stopifnot(type %in% names(object))

  # Special case of phi and  no newdata
  # We need a separate prediction for each detection dataset
  if(type == "phi" & missing(newdata)){

    dists <- names(object)[names(object) %in% c("ds", "pc", "oc")]
    out <- lapply(dists, function(x){
      conv <- IDS_convert_class(object, "phi", ds_type=x)
      predict(conv, "det", backTransform=backTransform, appendData=appendData,
              level=level, ...)
    })
    names(out) <- dists

  } else { # Regular situation
    conv <- IDS_convert_class(object, type)
    type <- switch(type, lam="state", ds="det", pc="det", oc="det", phi="det")
    out <- predict(conv, type=type, newdata=newdata, backTransform=backTransform, appendData=appendData,
                   level=level, ...)
  }
  out
})


# occuComm---------------------------------------------------------------------

setMethod("predict_internal", "unmarkedFitOccuComm",
  function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
           appendData = FALSE, level=0.95, re.form=NULL, ...){

  na.rm <- FALSE
  S <- length(object@data@ylist)
  M <- numSites(object@data)
  J <- obsNum(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula, object@data@speciesCovs)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")

  if(missing(newdata)) newdata <- NULL
  pr <- predict(new_object, type=type, newdata=newdata, backTransform=backTransform,
                na.rm=na.rm, appendData=appendData, level=level, re.form=re.form, ...)

  if(!is.null(newdata)){ # if using newdata, return now
    return(pr) 
  }

  if(type == "state"){
    inds <- split(1:nrow(pr), rep(1:length(object@data@ylist), each=M))
  } else if(type == "det"){
    inds <- split(1:nrow(pr), rep(1:length(object@data@ylist), each=M*J))
  }
  names(inds) <- names(object@data@ylist)
  lapply(inds, function(x){ 
         out <- pr[x,,drop=FALSE]
         rownames(out) <- NULL
         out
  })
})


# occuMulti--------------------------------------------------------------------

# bespoke predict method since it has numerious unusual options
# and requires bootstrapping
setMethod("predict", "unmarkedFitOccuMulti",
     function(object, type, newdata,
              level=0.95, species=NULL, cond=NULL, nsims=100,
              ...)
  {

  type <- match.arg(type, c("state", "det"))

  if(is.null(hessian(object))){
    level <- NULL
  }
  
  #For backwards compatability
  se.fit <- list(...)$se.fit
  if(is.null(se.fit)) se.fit <- TRUE

  species <- name_to_ind(species, names(object@data@ylist))
  cond <- name_to_ind(cond, names(object@data@ylist))

  if(missing(newdata)){
    newdata <- NULL
  } else {
    if(! class(newdata) %in% c('data.frame')){
      stop("newdata must be a data frame")
    }
  }

  maxOrder <- object@call$maxOrder
  if(is.null(maxOrder)) maxOrder <- length(object@data@ylist)
  dm <- getDesign(object@data, object@formlist,
                  maxOrder, na.rm=F, newdata=newdata, type=type)

  params <- coef(object)
  low_bound <- (1-level)/2
  up_bound <- level + (1-level)/2


  if(type=="state"){
    N <- nrow(dm$dmOcc[[1]]); nF <- dm$nF; dmOcc <- dm$dmOcc;
    fStart <- dm$fStart; fStop <- dm$fStop; fixed0 <- dm$fixed0
    t_dmF <- t(dm$dmF)

    calc_psi <- function(params){

      f <- matrix(NA,nrow=N,ncol=nF)
      index <- 1
      for (i in 1:nF){
        if(fixed0[i]){
          f[,i] <- 0
        } else {
          f[,i] <- dmOcc[[index]] %*% params[fStart[index]:fStop[index]]
          index <- index + 1
        }
      }
      psi <- exp(f %*% t_dmF)
      as.matrix(psi/rowSums(psi))
    }

    psi_est <- calc_psi(params)

    if(!is.null(level)){
      message('Bootstrapping confidence intervals with ',nsims,' samples')
      Sigma <- vcov(object)
      samp <- array(NA,c(dim(psi_est),nsims))
      for (i in 1:nsims){
        samp[,,i] <- calc_psi(mvrnorm(1, coef(object), Sigma))
      }
    }

    if(!is.null(species)){

      sel_col <- species

      if(!is.null(cond)){
        if(any(sel_col %in% abs(cond))){
          stop("Species can't be conditional on itself")
        }
        ftemp <- object@data@fDesign
        swap <- -1*cond[which(cond<0)]
        ftemp[,swap] <- 1 - ftemp[,swap]
        num_inds <- apply(ftemp[,c(sel_col,abs(cond))] == 1,1,all)
        denom_inds <- apply(ftemp[,abs(cond),drop=F] == 1,1,all)
        est <- rowSums(psi_est[,num_inds,drop=F]) /
          rowSums(psi_est[,denom_inds, drop=F])
        if(!is.null(level)){
          samp_num <- apply(samp[,num_inds,,drop=F],3,rowSums)
          samp_denom <- apply(samp[,denom_inds,,drop=F],3,rowSums)
          samp <- samp_num / samp_denom
        }

      } else {
        num_inds <- apply(object@data@fDesign[,sel_col,drop=FALSE] == 1,1,all)
        est <- rowSums(psi_est[,num_inds,drop=F])
        if(!is.null(level)){
          samp <- samp[,num_inds,,drop=F]
          samp <- apply(samp, 3, rowSums)
        }
      }

      if(!is.null(level)){
        if(!is.matrix(samp)) samp <- matrix(samp, nrow=1)
        boot_se <- apply(samp,1,sd, na.rm=T)
        boot_low <- apply(samp,1,quantile,low_bound, na.rm=T)
        boot_up <- apply(samp,1,quantile,up_bound, na.rm=T)
      } else{
        boot_se <- boot_low <- boot_up <- NA
      }
      return(data.frame(Predicted=est,
                        SE=boot_se,
                        lower=boot_low,
                        upper=boot_up))

    } else {
      codes <- apply(dm$z,1,function(x) paste(x,collapse=""))
      colnames(psi_est)  <- paste('psi[',codes,']',sep='')
      if(!is.null(level)){
        boot_se <- apply(samp,c(1,2),sd, na.rm=T)
        boot_low <- apply(samp,c(1,2),quantile,low_bound, na.rm=T)
        boot_up <- apply(samp,c(1,2),quantile,up_bound, na.rm=T)
        colnames(boot_se) <- colnames(boot_low) <- colnames(boot_up) <-
          colnames(psi_est)
      } else {
        boot_se <- boot_low <- boot_up <- NA
      }
      return(list(Predicted=psi_est,
                  SE=boot_se,
                  lower=boot_low,
                  upper=boot_up))
    }
  }

  if(type=="det"){
    S <- dm$S; dmDet <- dm$dmDet
    dStart <- dm$dStart; dStop <- dm$dStop

    out <- list()
    for (i in 1:S){
      #Subset estimate to species i
      inds <- dStart[i]:dStop[i]
      new_est <- object@estimates@estimates$det
      new_est@estimates <- coef(object)[inds]
      new_est@fixed <- 1:length(inds)
      if(!is.null(level)){
        new_est@covMat <- vcov(object)[inds,inds,drop=FALSE]
        new_est@covMatBS <- object@covMatBS[inds,inds,drop=FALSE]
      } else{
        new_est@covMat <- matrix(NA, nrow=length(inds), ncol=length(inds))
        new_est@covMatBS <- matrix(NA, nrow=length(inds), ncol=length(inds))
      }

      chunk_size <- 70
      xmat <- dmDet[[i]]
      if(is.vector(xmat)) xmat <- matrix(xmat, nrow=1)
      nr <- nrow(xmat)
      ind <- rep(1:ceiling(nr/chunk_size), each=chunk_size, length.out=nr)

      x_chunk <- lapply(unique(ind),
                    function(i) as.matrix(xmat[ind==i,,drop=FALSE]))

      prmat <- lapply(x_chunk, function(x_i){
        has_na <- apply(x_i, 1, function(x_i) any(is.na(x_i)))
        # Work around linearComb bug where there can't be NAs in inputs
        x_i[has_na,] <- 0
        lc <- linearComb(new_est, x_i)
        lc <- backTransform(lc)
        out <- data.frame(Predicted=coef(lc), SE=NA, lower=NA, upper=NA)
        if(!is.null(level)){
          se <- SE(lc)
          ci <- confint(lc, level=level)
          out$SE <- se
          out$lower <- ci[,1]
          out$upper <- ci[,2]
        }
        out[has_na,] <- NA
        out
      })
      prmat <- do.call(rbind, prmat)
      rownames(prmat) <- NULL
      out[[i]] <- as.data.frame(prmat)
    }
    names(out) <- names(object@data@ylist)
    if(!is.null(species)){
      return(out[[species]])
    }
    return(out)
  }
  stop("type must be 'det' or 'state'")
})


# occuMS-----------------------------------------------------------------------

# bespoke predict method since it requires bootstrapping
setMethod("predict", "unmarkedFitOccuMS",
     function(object, type, newdata, level=0.95, nsims=100, ...)
{

  #Process input---------------------------------------------------------------
  if(! type %in% c("psi","phi", "det")){
    stop("type must be 'psi', 'phi', or 'det'")
  }

  if(is.null(hessian(object))){
    level = NULL
  }

  #For backwards compatability
  se.fit <- list(...)$se.fit
  if(is.null(se.fit)) se.fit <- TRUE

  if(missing(newdata)){
    newdata <- NULL
  } else {
    if(! class(newdata) %in% c('data.frame')){
      stop("newdata must be a data frame")
    }
  }

  S <- object@data@numStates
  gd <- getDesign(object@data, object@formlist, object@parameterization, na.rm=F,
                  newdata=newdata, type=type)

  #Index guide used to organize p values
  guide <- matrix(NA,nrow=S,ncol=S)
  guide <- lower.tri(guide,diag=TRUE)
  guide[,1] <- FALSE
  guide <- which(guide,arr.ind=TRUE)
  #----------------------------------------------------------------------------

  #Utility functions-----------------------------------------------------------
  #Get matrix of linear predictor values
  get_lp <- function(params, dm_list, ind){
    L <- length(dm_list)
    out <- matrix(NA,nrow(dm_list[[1]]),L)
    for (i in 1:L){
      out[,i] <- dm_list[[i]] %*% params[ind[i,1]:ind[i,2]]
    }
    out
  }

  #Get SE/CIs for conditional binomial using delta method
  split_estimate <- function(object, estimate, inds, level){
    out <- estimate
    out@estimates <- coef(object)[inds]
    if(!is.null(level)){
      out@covMat <- vcov(object)[inds,inds,drop=FALSE]
    } else{
      out@covMat <- matrix(NA, nrow=length(inds), ncol=length(inds))
    }
    out
  }

  lc_to_predict <- function(object, estimate, inds, dm, level){

    new_est <- split_estimate(object, estimate, inds[1]:inds[2], level)

    out <- t(apply(dm, 1, function(x){
      bt <- backTransform(linearComb(new_est, x))
      if(is.null(level)) return(c(Predicted=bt@estimate, SE=NA, lower=NA, upper=NA))
      ci <- confint(bt, level=level)
      names(ci) <- c("lower", "upper")
      c(Predicted=bt@estimate, SE=SE(bt), ci)
    }))
    rownames(out) <- NULL
    as.data.frame(out)
  }


  #Calculate row-wise multinomial logit prob
  #implemented in C++ below as it is quite slow
  get_mlogit_R <- function(lp_mat){
    if(type == 'psi'){
      out <- cbind(1,exp(lp_mat))
      out <- out/rowSums(out)
      out <- out[,-1]
    } else if(type == 'phi'){ #doesn't work
      np <- nrow(lp_mat)
      out <- matrix(NA,np,ncol(lp_mat))
      ins <- outer(1:S, 1:S, function(i,j) i!=j)
      for (i in 1:np){
        phimat <- diag(S)
        phimat[ins] <- exp(lp_mat[i,])
        phimat <- t(phimat)
        phimat <- phimat/rowSums(phimat)
        out[i,] <- phimat[ins]
      }
    } else {
      R <- nrow(lp_mat)
      out <- matrix(NA,R,ncol(lp_mat))
      for (i in 1:R){
        sdp <- matrix(0,nrow=S,ncol=S)
        sdp[guide] <- exp(lp_mat[i,])
        sdp[,1] <- 1
        sdp <- sdp/rowSums(sdp)
        out[i,] <- sdp[guide]
      }
    }
    out
  }

  get_mlogit_C <- function(lp_mat){
    get_mlogit(lp_mat, type, S, guide-1) # via Rcpp
  }

  #----------------------------------------------------------------------------

  if(type=="psi"){
    dm_list <- gd$dm_state
    ind <- gd$state_ind
    est <- object@estimates@estimates$state
  } else if(type=="phi"){
    dm_list <- gd$dm_phi
    ind <- gd$phi_ind
    est <- object@estimates@estimates$transition
  } else {
    dm_list <- gd$dm_det
    ind <- gd$det_ind
    est <- object@estimates@estimates$det
  }

  P <- length(dm_list)

  low_bound <- (1-level)/2
  z <- qnorm(low_bound,lower.tail=F)

  out <- vector("list", P)
  names(out) <- names(dm_list)

  if(object@parameterization == 'condbinom'){
    out <- lapply(1:length(dm_list), function(i){
      lc_to_predict(object, est, ind[i,], dm_list[[i]], level)
    })
    names(out) <- names(dm_list)
    return(out)

  } else if (object@parameterization == "multinomial"){
    lp <- get_lp(coef(object), dm_list, ind)
    pred <- get_mlogit_C(lp)

    M <- nrow(pred)
    upr <- lwr <- se <- matrix(NA,M,P)

    if(!is.null(level)){
      message('Bootstrapping confidence intervals with',nsims,'samples')

      sig <- vcov(object)
      param_mean <- coef(object)
      rparam <- mvrnorm(nsims, param_mean, sig)

      get_pr <- function(i){
        lp <- get_lp(rparam[i,], dm_list, ind)
        get_mlogit_C(lp)
      }
      samp <- sapply(1:nsims, get_pr, simplify='array')

      for (i in 1:M){
        for (j in 1:P){
          dat <- samp[i,j,]
          se[i,j] <- sd(dat, na.rm=TRUE)
          quants <- quantile(dat, c(low_bound, (1-low_bound)),na.rm=TRUE)
          lwr[i,j] <- quants[1]
          upr[i,j] <- quants[2]
        }
      }

    }
  }

  for (i in 1:P){
    out[[i]] <- data.frame(Predicted=pred[,i], SE=se[,i],
                           lower=lwr[,i], upper=upr[,i])
  }

  out
})
