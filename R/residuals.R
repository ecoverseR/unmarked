# Method for all fit types, exported
setMethod("residuals", "unmarkedFit", function(object, ...){
  residuals_internal(object)
})

# Internal method for specific fit types, not exported
setMethod("residuals_internal", "unmarkedFit", function(object) {
    y <- getY(object)
    e <- fitted(object)
    r <- y - e
    return(r)
})

setMethod("residuals_internal", "unmarkedFitGDR", function(object){
  ft <- fitted(object)
  list(dist=object@data@yDistance - ft$dist, rem=object@data@yRemoval-ft$rem)
})

setMethod("residuals_internal", "unmarkedFitIDS", function(object){

  dists <- names(object)[names(object) %in% c("ds", "pc")]

  # distance and N-mix data
  out <- lapply(dists, function(x){
    conv <- IDS_convert_class(object, type=x)
    residuals(conv)
  })
  names(out) <- dists

  # occupancy data
  if("oc" %in% names(object)){
    y <- object@dataOC@y
    ft <- fitted(object)$oc
    out$oc <- y - ft
  }

  out
})

setMethod("residuals_internal", "unmarkedFitOccuComm", function(object) {
  ylist <- getY(object)
  fitlist <- fitted(object)

  mapply(function(x, y){
    x - y
  }, x = ylist, y = fitlist, SIMPLIFY = FALSE)
})

setMethod("residuals_internal", "unmarkedFitOccuMulti", function(object) {
  res_list <- list()
  ylist <- getY(object)
  fitlist <- fitted(object)

  for (i in seq_along(ylist)){
    res_list[[i]] <- ylist[[i]] - fitlist[[i]]
  }
  names(res_list) <- names(ylist)
  res_list
})

setMethod("residuals_internal", "unmarkedFitOccuTTD", function(object) {
  tmax <- object@data@surveyLength
  yraw <- object@data@y
  y <- ifelse(yraw<tmax,1,0)
  e <- fitted(object)
  y - e
})
