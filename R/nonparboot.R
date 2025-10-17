# Overall nonparboot method; this is the one exported
setMethod("nonparboot", "unmarkedFit", function(object, B = 1, keepOldSamples = TRUE, ...){
  bsType <- list(...)$bsType
  if(!is.null(bsType) && bsType == "both"){
    warning("Only bsType = 'site' is supported")
  }
  stopifnot(B > 0)
  # Calls fit-specific internal method
  nonparboot_internal(object, B=B, keepOldSamples=keepOldSamples)
})

# Fit-specific internal methods
setMethod("nonparboot_internal", "unmarkedFit", 
          function(object, B, keepOldSamples){

  M <- numSites(object@data)
  
  # Note: no sites are removed. I think this is more accurate as it
  # results in comparable sample sizes to original data.
  boot_iter <- lapply2(1:B, function(i){
    finish <- FALSE
    while(!finish){
      sites <- sort(sample(1:M, M, replace=TRUE))
      new_data <- object@data[sites,]
      ran <- TRUE
      tryCatch(fit <- nonparboot_update(object, data = new_data),
               error = function(e) ran <<- FALSE)
      conv <- object@opt$convergence == 0
      if(ran & conv) finish <- TRUE
    }
    fit
  })

  if(!keepOldSamples) object@bootstrapSamples <- NULL
  object@bootstrapSamples <- c(object@bootstrapSamples, boot_iter)
   
  coefs <- t(sapply(object@bootstrapSamples, coef))
  v <- stats::cov(coefs)
  object@covMatBS <- v

  short_names <- regmatches(colnames(v), regexpr("^[^\\(]+", colnames(v)))
  
  for (i in 1:length(object@estimates@estimates)){
    match_name <- object@estimates@estimates[[i]]@short.name
    inds <- short_names == match_name
    new_v <- v[inds, inds, drop=FALSE]
    object@estimates@estimates[[i]]@covMatBS <- new_v
  }

  object
})

# Fit-specific update methods (to handle occuPEN)
setMethod("nonparboot_update", "unmarkedFit", function(object, data){
  update(object, data = data, se = FALSE)
})

setMethod("nonparboot_update", "unmarkedFitOccuPEN", function(object, data){
  update(object, data = data)
})

setMethod("nonparboot_update", "unmarkedFitOccuPEN_CV", function(object, data){
  if(object@pen.type == "MPLE"){
	  MPLElambda <- computeMPLElambda(object@formula, data)
    out <- update(object, data = data, lambda = MPLElambda)
  } else {
    out <- update(object, data = data)
  }
  out
})

# Colext has some additional project stuff
setMethod("nonparboot_internal", "unmarkedFitColExt", 
          function(object, B, keepOldSamples){

  # Call base method
  object <- methods::callNextMethod(object, B=B, keepOldSamples=keepOldSamples)
  
  # Add smoothed/projected stuff
  smoothed.occ <- t(sapply(object@bootstrapSamples,
                    function(x) x@smoothed.mean[1,]))
  smoothed.unocc <- t(sapply(object@bootstrapSamples,
                      function(x) x@smoothed.mean[2,]))
  object@smoothed.mean.bsse <- rbind(sqrt(diag(cov(smoothed.occ))),
                                     sqrt(diag(cov(smoothed.unocc))))
  projected.occ <- t(sapply(object@bootstrapSamples,
                     function(x) x@projected.mean[1,]))
  projected.unocc <- t(sapply(object@bootstrapSamples,
                    function(x) x@projected.mean[2,]))
  object@projected.mean.bsse <- rbind(sqrt(diag(cov(projected.occ))),
                                      sqrt(diag(cov(projected.unocc))))
  object
})
