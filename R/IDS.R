get_ds_info <- function(db){
  J <- length(db) - 1
  a <- u <- rep(NA, J)
  a[1] <- pi*db[2]^2
  if(J > 1){
    for (j in 2:J){
      a[j] <- pi*db[j+1]^2 - sum(a[1:(j-1)])
    }
  }

  total_area <- sum(a)
  u <- a/total_area
  w <- diff(db)
  list(total_area=total_area, a=a, u=u, w=w)
}

IDS <- function(lambdaformula = ~1,
                detformulaDS = ~1,
                detformulaPC = NULL,
                detformulaOC = NULL,
                dataDS,
                dataPC = NULL,
                dataOC = NULL,
                availformula = NULL,
                durationDS = NULL,
                durationPC = NULL,
                durationOC = NULL,
                keyfun = "halfnorm",
                maxDistPC,
                maxDistOC,
                K = 100,
                unitsOut = "ha",
                starts = NULL,
                method = "BFGS",
                ...
                ){

  # Process inputs-------------------------------------------------------------

  form_hds <- as.formula(paste(c(as.character(detformulaDS),
                                 as.character(lambdaformula)), collapse=""))

  if(is.null(detformulaPC)){
    form_pc <- as.formula(paste(c(as.character(detformulaDS),
                                 as.character(lambdaformula)), collapse=""))
  } else{
    form_pc <- as.formula(paste(c(as.character(detformulaPC),
                                 as.character(lambdaformula)), collapse=""))
  }

  if(is.null(detformulaOC)){
    form_oc <- as.formula(paste(c(as.character(detformulaDS),
                                 as.character(lambdaformula)), collapse=""))
  } else{
    form_oc <- as.formula(paste(c(as.character(detformulaOC),
                                 as.character(lambdaformula)), collapse=""))
  }

  formlist <- list(lam=lambdaformula, ds=form_hds, pc=form_pc, oc=form_oc,
                   phi=availformula)

  stopifnot(inherits(dataDS, "unmarkedFrameDS"))
  stopifnot(inherits(dataPC, c("unmarkedFramePCount", "NULL")))
  stopifnot(inherits(dataOC, c("unmarkedFrameOccu", "NULL")))
  #stopifnot(!is.null(dataPC) | !is.null(dataOC))

  has_avail <- FALSE
  if(!is.null(durationDS) | !is.null(durationPC) | !is.null(durationOC)){
    has_avail <- TRUE
    if(is.null(availformula)) availformula <- ~1
    form_avail <- as.formula(paste("~1", paste(as.character(availformula), collapse="")))
    stopifnot(!is.null(durationDS))
    if(!is.null(dataPC)) stopifnot(!is.null(durationPC))
    if(!is.null(dataOC)) stopifnot(!is.null(durationOC))
  }

  if(has_avail & !is.null(dataOC)){
    stop("Availability estimation doesn't work with detection-nondetection data", call.=FALSE)
  }

  stopifnot(is.null(durationDS) || (length(durationDS) == numSites(dataDS)))
  stopifnot(is.null(durationPC) || (length(durationPC) == numSites(dataPC)))
  stopifnot(is.null(durationOC) || (length(durationOC) == numSites(dataOC)))
  surveyDurations <- list(ds=durationDS, pc=durationPC, oc=durationOC)

  stopifnot(keyfun %in% c("halfnorm", "exp", "hazard"))
  if(keyfun == "hazard"){
    warning("Support for hazard key function is experimental and may yield poor estimates.", call.=FALSE)
  }
  keyidx <- switch(keyfun, "halfnorm"={1}, "exp"={2}, "hazard"={3})

  if(missing(maxDistPC)) maxDistPC <- max(dataDS@dist.breaks)
  if(missing(maxDistOC)) maxDistOC <- max(dataDS@dist.breaks)


  # Design matrices------------------------------------------------------------

  # Need to add offset support here eventually
  form_hds_split <- split_formula(form_hds)
  names(form_hds_split) <- c("det", "state")
  gd_hds <- getDesign(dataDS, form_hds_split)
  if(any(is.na(gd_hds$y))){
    stop("Missing values in only some distance bins is not supported", call.=FALSE)
  } 
  ds_hds <- get_ds_info(dataDS@dist.breaks)
  Xavail_ds <- matrix(0,0,0)
  if(has_avail){
    form_avail_split <- split_formula(form_avail)
    names(form_avail_split) <- c("det", "state")
    Xavail_ds <- getDesign(dataDS, form_avail_split)$X_state
  }
  if(is.null(durationDS)) durationDS <- rep(0,0)

  gd_pc <- list(y=matrix(0,0,0), X_state=matrix(0,0,0), X_det=matrix(0,0,0))
  ds_pc <- list(total_area=0, db=c(0,0), a=0, w=0, u=0)
  Xavail_pc <- matrix(0,0,0)
  if(is.null(durationPC)) durationPC <- rep(0,0)
  if(!is.null(dataPC)){
    form_pc_split <- split_formula(form_pc)
    names(form_pc_split) <- c("det", "state")
    gd_pc <- getDesign(dataPC, form_pc_split)
    ds_pc <- get_ds_info(c(0, maxDistPC))
    if(has_avail) Xavail_pc <- getDesign(dataPC, form_avail_split)$X_state
  }

  gd_oc <- list(y=matrix(0,0,0), X_state=matrix(0,0,0), X_det=matrix(0,0,0))
  ds_oc <- list(total_area=0, db=c(0,0), a=0, w=0, u=0)
  Kmin_oc <- rep(0,0)
  Xavail_oc <- matrix(0,0,0)
  if(is.null(durationOC)) durationOC <- rep(0,0)
  if(!is.null(dataOC)){
    form_oc_split <- split_formula(form_pc)
    names(form_oc_split) <- c("det", "state")
    gd_oc <- getDesign(dataOC, form_oc_split)
    ds_oc <- get_ds_info(c(0, maxDistOC))
    Kmin_oc <- apply(gd_oc$y, 1, max, na.rm=T)
    if(has_avail) Xavail_oc <- getDesign(dataOC, form_avail_split)$X_state
  }

  # Density conversion and unequal area correction
  lam_adjust <- c(ds_hds$total_area, ds_pc$total_area, ds_oc$total_area)
  names(lam_adjust) <- c("hds", "pc", "oc")

  switch(dataDS@unitsIn,
    m = lam_adjust <- lam_adjust / 1e6,
    km = lam_adjust <- lam_adjust)

  stopifnot(unitsOut %in% c("m","ha","kmsq"))
  switch(unitsOut,
    m = lam_adjust <- lam_adjust * 1e6,
    ha = lam_adjust <- lam_adjust * 100,
    kmsq = lam_adjust <- lam_adjust)

  # Parameter stuff------------------------------------------------------------
  pind_mat <- matrix(0, nrow=8, ncol=2)
  pind_mat[1,] <- c(1, ncol(gd_hds$X_state))
  pind_mat[2,] <- max(pind_mat) + c(1, ncol(gd_hds$X_det))
  if(!is.null(detformulaPC) & !is.null(dataPC)){
    pind_mat[3,] <- max(pind_mat) + c(1, ncol(gd_pc$X_det))
  }
  if(!is.null(detformulaOC) & !is.null(dataOC)){
    pind_mat[4,] <- max(pind_mat) + c(1, ncol(gd_oc$X_det))
  }
  if(has_avail){
    pind_mat[5,] <- max(pind_mat) + c(1, ncol(Xavail_ds))
  }
  # Hazard-rate scale parameter(s)
  if(keyfun == "hazard"){
    pind_mat[6,] <- max(pind_mat) + 1
    if(!is.null(detformulaPC) & !is.null(dataPC)){
      pind_mat[7,] <- max(pind_mat) + 1
    }
    if(!is.null(detformulaOC) & !is.null(dataOC)){
      pind_mat[8,] <- max(pind_mat) + 1
    }
  }

  if(is.null(starts)){
    lam_init <- log(mean(apply(dataDS@y, 1, sum, na.rm=TRUE)) / lam_adjust[1])
    params_tmb <- list(beta_lam = c(lam_init, rep(0, ncol(gd_hds$X_state)-1)),
                     beta_hds = c(log(median(dataDS@dist.breaks)),rep(0, ncol(gd_hds$X_det)-1)),
                     beta_pc = rep(0,0),
                     beta_oc = rep(0,0),
                     beta_avail = rep(0,0),
                     beta_schds = rep(0,0),  # hazard rate scale inits
                     beta_scpc = rep(0,0), 
                     beta_scoc = rep(0,0))

    if(keyfun == "hazard") params_tmb$beta_schds <- 0
    
    if(!is.null(detformulaPC) & !is.null(dataPC)){
      params_tmb$beta_pc <- c(log(maxDistPC/2), rep(0, ncol(gd_pc$X_det)-1))
      if(keyfun == "hazard") params_tmb$beta_scpc <- 0
    }
    if(!is.null(detformulaOC) & !is.null(dataOC)){
      params_tmb$beta_oc <- c(log(maxDistOC/2), rep(0, ncol(gd_oc$X_det)-1))
      if(keyfun == "hazard") params_tmb$beta_scoc <- 0
    }
    if(has_avail){
      params_tmb$beta_avail <- rep(0, ncol(Xavail_ds))
    }
    starts <- unlist(params_tmb)
  } else {
    if(length(starts) != max(pind_mat)){
      stop("Length of starts should be ", max(pind_mat), call.=FALSE)
    }
    params_tmb <- list(beta_lam = starts[pind_mat[1,1]:pind_mat[1,2]],
                     beta_hds = starts[pind_mat[2,1]:pind_mat[2,2]],
                     beta_pc = rep(0,0),
                     beta_oc = rep(0,0),
                     beta_avail = rep(0,0),
                     beta_schds = rep(0,0),  # hazard rate scale inits
                     beta_scpc = rep(0,0), 
                     beta_scoc = rep(0,0))
    
    if(keyfun == "hazard") params_tmb$beta_schds <- pind_mat[6,1]

    if(!is.null(detformulaPC) & !is.null(dataPC)){
      params_tmb$beta_pc <- starts[pind_mat[3,1]:pind_mat[3,2]]
      if(keyfun == "hazard") params_tmb$beta_scpc <- starts[pind_mat[7,1]]
    }
    if(!is.null(detformulaOC) & !is.null(dataOC)){
      params_tmb$beta_oc <- starts[pind_mat[4,1]:pind_mat[4,2]]
      if(keyfun == "hazard") params_tmb$beta_scoc <- starts[pind_mat[8,1]]
    }
    if(has_avail){
      params_tmb$beta_avail <- starts[pind_mat[5,1]:pind_mat[5,2]]
    }
  }

  # Construct TMB data list----------------------------------------------------
  tmb_dat <- list(
    # Common
    pind = pind_mat, lam_adjust = lam_adjust,

    # HDS data
    y_hds = gd_hds$y, X_hds = gd_hds$X_state, V_hds = gd_hds$X_det, key_hds = keyidx,
    db_hds = dataDS@dist.breaks, a_hds = ds_hds$a, w_hds = ds_hds$w,
    u_hds = ds_hds$u,

    # PC data
    y_pc = gd_pc$y, X_pc = gd_pc$X_state, V_pc = gd_pc$X_det, key_pc = keyidx,
    db_pc = c(0, maxDistPC), a_pc = ds_pc$a, w_pc = ds_pc$w, u_pc = ds_pc$u,

    # occ data
    y_oc = gd_oc$y, X_oc = gd_oc$X_state, V_oc = gd_oc$X_det, key_oc = keyidx,
    db_oc = c(0, maxDistOC), a_oc = ds_oc$a, w_oc = ds_oc$w, u_oc = ds_oc$u,
    K_oc = K, Kmin_oc = Kmin_oc,

    # avail data
    durationDS = durationDS, durationPC = durationPC, durationOC = durationOC,
    Xa_hds = Xavail_ds, Xa_pc = Xavail_pc, Xa_oc = Xavail_oc
  )

  tmb_obj <- TMB::MakeADFun(data = c(model = "tmb_IDS", tmb_dat), parameters = params_tmb,
                            DLL = "unmarked_TMBExports", silent=TRUE)

  opt <- optim(unlist(params_tmb), fn=tmb_obj$fn, gr=tmb_obj$gr, method=method, ...)

  fmAIC <- 2 * opt$value + 2 * length(unlist(params_tmb))

  sdr <- TMB::sdreport(tmb_obj)

  lam_coef <- get_coef_info(sdr, "lam", colnames(gd_hds$X_state),
                                       pind_mat[1,1]:pind_mat[1,2])

  lam_est <- unmarkedEstimate(name="Density", short.name="lam",
    estimates = lam_coef$ests, covMat = lam_coef$cov, fixed=1:ncol(gd_hds$X_state),
    invlink = "exp", invlinkGrad = "exp")

  dist_coef <- get_coef_info(sdr, "hds", colnames(gd_hds$X_det),
                                       pind_mat[2,1]:pind_mat[2,2])

  dist_est <- unmarkedEstimate(name="Distance sampling detection", short.name="ds",
    estimates = dist_coef$ests, covMat = dist_coef$cov, fixed=1:ncol(gd_hds$X_det),
    invlink = "exp", invlinkGrad = "exp")

  est_list <- list(lam=lam_est, ds=dist_est)

  if(!is.null(detformulaPC) & !is.null(dataPC)){
    pc_coef <- get_coef_info(sdr, "pc", colnames(gd_pc$X_det),
                                        pind_mat[3,1]:pind_mat[3,2])

    pc_est <- unmarkedEstimate(name="Point count detection", short.name="pc",
      estimates = pc_coef$ests, covMat = pc_coef$cov, fixed=1:ncol(gd_pc$X_det),
      invlink = "exp", invlinkGrad = "exp")
    est_list <- c(est_list, list(pc=pc_est))
  }

  if(!is.null(detformulaOC) & !is.null(dataOC)){
    oc_coef <- get_coef_info(sdr, "oc", colnames(gd_oc$X_det),
                                        pind_mat[4,1]:pind_mat[4,2])
    oc_est <- unmarkedEstimate(name="Presence/absence detection", short.name="oc",
      estimates = oc_coef$ests, covMat = oc_coef$cov, fixed=1:ncol(gd_oc$X_det),
      invlink = "exp", invlinkGrad = "exp")
    est_list <- c(est_list, list(oc=oc_est))
  }

  if(has_avail){
    avail_coef <- get_coef_info(sdr, "avail", colnames(Xavail_ds),
                                pind_mat[5,1]:pind_mat[5,2])
    avail_est <- unmarkedEstimate(name="Availability", short.name="phi",
      estimates=avail_coef$ests, covMat=avail_coef$cov, fixed=1:ncol(Xavail_ds),
      invlink="exp", invlinkGrad="exp")
    est_list <- c(est_list, list(phi=avail_est))
  }

  if(keyfun == "hazard"){
    ds_haz_coef <- get_coef_info(sdr, "schds", "(Intercept)", pind_mat[6,1]:pind_mat[6,2])
    ds_est_haz <- unmarkedEstimate("Distance sampling scale", short.name="ds_scale",
      estimates = ds_haz_coef$ests[1], covMat = ds_haz_coef$cov, fixed=1,
      invlink="exp", invlinkGrad="exp")
    est_list <- c(est_list, list(schds=ds_est_haz))

    if(!is.null(detformulaPC) & !is.null(dataPC)){
      pc_haz_coef <- get_coef_info(sdr, "scpc", "(Intercept)", pind_mat[7,1]:pind_mat[7,2])
      pc_est_haz <- unmarkedEstimate("Point count scale", short.name="pc_scale",
        estimates = pc_haz_coef$ests[1], covMat = pc_haz_coef$cov, fixed=1,
        invlink="exp", invlinkGrad="exp")
      est_list <- c(est_list, list(scpc=pc_est_haz))
    }

    if(!is.null(detformulaOC) & !is.null(dataOC)){
      oc_haz_coef <- get_coef_info(sdr, "scoc", "(Intercept)", pind_mat[8,1]:pind_mat[8,2])
      oc_est_haz <- unmarkedEstimate("Presence/absence scale", short.name="oc_scale",
        estimates = oc_haz_coef$ests[1], covMat = oc_haz_coef$cov, fixed=1,
        invlink="exp", invlinkGrad="exp")
      est_list <- c(est_list, list(scoc=oc_est_haz))
    }
  }

  est_list <- unmarkedEstimateList(est_list)

  new("unmarkedFitIDS", fitType = "IDS", call = match.call(),
    opt = opt, formlist=formlist,
    data = dataDS, dataPC=dataPC, dataOC=dataOC, K=K,
    surveyDurations=surveyDurations,
    maxDist = list(pc=maxDistPC, oc=maxDistOC),
    keyfun=keyfun,
    # this needs to be fixed
    sitesRemoved = gd_hds$removed.sites,
    unitsOut=unitsOut,
    estimates = est_list, AIC = fmAIC, negLogLike = opt$value,
    nllFun = tmb_obj$fn, TMB=tmb_obj)

}

setMethod("summary_internal", "unmarkedFitIDS", function(object)
{
    cat("\nCall:\n")
    print(object@call)
    cat("\n")
    summaryOut <- summary(object@estimates)
    cat("AIC:", object@AIC,"\n")
    cat("Number of distance sampling sites:", numSites(object@data))
    if(!is.null(object@dataPC)){
      cat("\nNumber of point count sites:", numSites(object@dataPC))
    }
    if(!is.null(object@dataOC)){
      cat("\nNumber of presence/absence sites:", numSites(object@dataOC))
    }
    cat("\noptim convergence code:", object@opt$convergence)
    cat("\noptim iterations:", object@opt$counts[1], "\n")
    if(!identical(object@opt$convergence, 0L))
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
    cat("Bootstrap iterations:", length(object@bootstrapSamples), "\n\n")
    invisible(summaryOut)
})

# Need special method since an included dataset may not have a corresponding
# unique submodel in the unmarked estimates
setMethod("names", "unmarkedFitIDS", function(x){
  out <- c("lam","ds")
  if(!is.null(x@dataPC)) out <- c(out, "pc")
  if(!is.null(x@dataOC)) out <- c(out, "oc")
  if("phi" %in% names(x@estimates)) out <- c(out, "phi")
  out
})

# This function converts IDS objects into simpler distsamp objects
# so we can re-use distsamp methods
# the abundance model (lam) is combined with one of the detection models
# to yield a distsamp model
IDS_convert_class <- function(inp, type, ds_type=NULL){
  stopifnot(type %in% names(inp))
  if(is.null(ds_type)) ds_type <- type
  if(type == "lam") type <- "ds"
  data <- inp@data
  if(ds_type == "pc"){
    tempdat <- inp@dataPC
    maxDist <- inp@maxDist$pc
  }
  if(ds_type == "oc"){
    tempdat <- inp@dataOC
    maxDist <- inp@maxDist$oc
  }
  if(ds_type %in% c("pc", "oc")){
    if(is.null(maxDist)) maxDist <- max(data@dist.breaks)
    data <- unmarkedFrameDS(y=tempdat@y, siteCovs=siteCovs(tempdat),
                            dist.breaks=c(0, maxDist), survey="point",
                            unitsIn=data@unitsIn)
  }

  # Select appropriate model estimates; if sigma was not estimated
  # separately for this model, pick the DS model for detection
  det_mod <- type
  if(! det_mod %in% names(inp@estimates)) det_mod <- "ds"
  est <- inp@estimates@estimates[c("lam", det_mod)]
  names(est) <- c("state", "det")

  form <- inp@formlist[[type]]
  if(type=="phi") form <- as.formula(paste(c(as.character(form), "~1"), collapse=""))
  formlist <- split_formula(form)
  names(formlist) <- c("det", "state")

  new("unmarkedFitDS", fitType="IDS", opt=inp@opt, formula=form, formlist=formlist,
      data=data, keyfun=inp@keyfun, unitsOut=inp@unitsOut,
      estimates=unmarkedEstimateList(est),
      AIC=inp@AIC, output="density", TMB=inp@TMB)
}

# Get availability for each data type and site as a probability
setMethod("getAvail", "unmarkedFitIDS", function(object, ...){
  stopifnot("phi" %in% names(object))
  phi <- predict(object, "phi", level = NULL, na.rm=FALSE)
  dur <- object@surveyDurations
  out <- lapply(names(phi), function(x){
    1 - exp(-1 * dur[[x]] * phi[[x]]$Predicted)
  })
  names(out) <- names(phi)
  out
})

# getP returns detection probability WITHOUT availability
setMethod("getP_internal", "unmarkedFitIDS", function(object){

  dets <- names(object)[! names(object) %in% c("lam","phi")]

  out <- lapply(dets, function(x){
    conv <- IDS_convert_class(object, type=x)
    getP(conv)
  })
  names(out) <- dets
  out
})

setMethod("hist", "unmarkedFitIDS", function(x, lwd=1, lty=1, ...){

  conv <- IDS_convert_class(x, type='ds')
  hist(conv, lwd=lwd, lty=lty, ...)

})

setMethod("residual_plot", "unmarkedFitIDS", function(x, ...){

  r <- residuals(x)
  f <- fitted(x)
  nr <- length(r)
  long_names <- sapply(x@estimates@estimates, function(x) x@name)
  long_names <- long_names[long_names != "Density"]

  old_par <- graphics::par()$mfrow
  graphics::par(mfrow=c(nr,1))

  for (i in 1:nr){
    plot(f[[i]], r[[i]], ylab = "Residuals", xlab = "Predicted values",
         main=long_names[i])
    abline(h = 0, lty = 3, col = "gray")
  }
  graphics::par(mfrow=old_par)

})


setMethod("simulate_internal", "unmarkedFitIDS",
          function(object,  nsim){

  dets <- c("ds","pc","oc")

  if("phi" %in% names(object)) avail <- getAvail(object)

  temp <- lapply(dets, function(x){
    if(! x %in% names(object)) return(NULL)
    conv <- IDS_convert_class(object, type=x)
    sims <- simulate(conv, nsim=nsim, na.rm=FALSE)
    # availability process
    if("phi" %in% names(object)){
      sims <- lapply(sims, function(z){
        avail_expand <- matrix(rep(avail[[x]], ncol(z)), ncol=ncol(z))
        out <- rbinom(length(z), z, avail_expand)
        matrix(out, nrow=nrow(z), ncol=ncol(z))
      })
    }
    if(x=="oc"){
      sims <- lapply(sims, function(z){
                    z[z>1] <- 1
                    z
                  })
    }
    sims
  })

  #"permute"
  lapply(1:nsim, function(i){
    sim <- lapply(temp, function(x) x[[i]])
    names(sim) <- c("ds","pc","oc")
    sim
  })

})


setMethod("rebuild_call", "unmarkedFitIDS", function(object){           
  cl <- object@call
  cl[["dataDS"]] <- quote(object@data)
  cl[["dataPC"]] <- quote(object@dataPC)
  cl[["dataOC"]] <- quote(object@dataOC)
  cl[["lambdaformula"]] <- object@formlist$lam
  cl[["detformulaDS"]] <- split_formula(object@formlist$ds)[[1]]
  if(!is.null(cl[["detformulaPC"]])){
    cl[["detformulaPC"]] <- split_formula(object@formlist$pc)[[1]]
  }
  if(!is.null(cl[["detformulaOC"]])){
    cl[["detformulaOC"]] <- split_formula(object@formlist$oc)[[1]]
  }
  if(!is.null(cl[["availformula"]])){
    cl[["availformula"]] <- object@formlist$phi
  }
  cl[["K"]] <- object@K
  cl[["keyfun"]] <- object@keyfun
  cl[["unitsOut"]] <- object@unitsOut
  cl[["durationDS"]] <- quote(object@surveyDurations$ds)
  cl[["durationPC"]] <- quote(object@surveyDurations$pc)
  cl[["durationOC"]] <- quote(object@surveyDurations$oc)
  cl[["maxDistPC"]] <- object@maxDist$pc
  cl[["maxDistOC"]] <- object@maxDist$oc
  cl
})


setMethod("SSE", "unmarkedFitIDS", function(fit, ...){
  stop("Not currently supported for unmarkedFitIDS", call.=FALSE)
})

setMethod("replaceY", "unmarkedFitIDS", function(object, newY, replNA, ...){
  stop("Not currently supported for unmarkedFitIDS", call.=FALSE)
})

setMethod("nonparboot_internal", "unmarkedFitIDS",
    function(object, B, keepOldSamples)
{
   stop("Not currently supported for unmarkedFitIDS", call.=FALSE)
})

setMethod("ranef_internal", "unmarkedFitIDS",
    function(object, ...)
{
   stop("Not currently supported for unmarkedFitIDS", call.=FALSE)
})
