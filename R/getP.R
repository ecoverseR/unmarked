# getP creates a matrix of estimated detection probabilities corresponding 
# to the observed data. For example if the data are MxJ detection/non-detection,
# the result will be an MxJ matrix of detection probabilities.

# Exported method
# This method applies to all unmarkedFit types, then internally
# calls a specific method for each unmarkedFit type (minimizes documentation)
setMethod("getP", "unmarkedFit", function(object, ...){
  getP_internal(object)
})

setMethod("getP_internal", "unmarkedFit", function(object){
  M <- numSites(object@data)
  J <- obsNum(object@data)
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  p <- matrix(p, M, J, byrow = TRUE)
  p
})

setMethod("getP_internal", "unmarkedFitColExt", function(object){
  data <- object@data
  M <- numSites(data)
  nY <- data@numPrimary
  J <- obsNum(data)/nY
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  p <- array(p, c(J, nY, M))
  p <- aperm(p, c(3, 1, 2))
  p <- matrix(p, nrow=M)
  p
})

setMethod("getP_internal", "unmarkedFitDS", function(object){
  cp <- get_dist_prob(object)
  cp
})

setMethod("getP_internal", "unmarkedFitDSO", function(object){
  cp <- get_dist_prob(object)
  cp
})

# Should this return p or pi. Right now it's pi without phi.
setMethod("getP_internal", "unmarkedFitGDS", function(object){
  cp <- get_dist_prob(object)
  cp
})

setMethod("getP_internal", "unmarkedFitGMM", function(object){
  piFun <- object@data@piFun
  M <- numSites(object@data)
  T <- object@data@numPrimary
  R <- numY(object@data) / T
  J <- obsNum(object@data) / T

  p <- predict(object, type = "det", level=NULL, na.rm=FALSE)$Predicted
  p <- matrix(p, nrow=M, byrow=TRUE)
  p <- array(p, c(M, J, T))
  p <- aperm(p, c(1,3,2))

  cp <- array(as.numeric(NA), c(M, T, R))
  for(t in 1:T) cp[,t,] <- do.call(piFun, list(matrix(p[,t,], M, J)))
  cp <- aperm(cp, c(1,3,2))
  cp <- matrix(cp, nrow=M, ncol=numY(object@data))
  cp
})

setMethod("getP_internal", "unmarkedFitGPC", function(object){
  M <- numSites(object@data)
  R <- ncol(object@data@y)
  
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  p <- matrix(p, M, R, byrow=TRUE)
  p
})

setMethod("getP_internal", "unmarkedFitMMO", function(object){
  umf <- object@data
  M <- numSites(umf)
  T <- umf@numPrimary
  J <- ncol(getY(umf)) / T

  plong <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  pmat <- aperm(array(plong, c(J,T,M)), c(3,1,2))

  pout <- array(NA, c(M,J,T))
  for (t in 1:T){
    pout[,,t] <- do.call(umf@piFun, list(p=pmat[,,t]))
  }
  matrix(aperm(pout,c(2,3,1)), M, J*T, byrow=TRUE)
})

setMethod("getP_internal", "unmarkedFitMPois", function(object){
  p <- methods::callNextMethod(object)
  piFun <- object@data@piFun
  pi <- do.call(piFun, list(p = p))
  pi
})

setMethod("getP_internal", "unmarkedFitNmixTTD", function(object){
  stop("getP is not implemented for nmixTTD at this time", call.=FALSE)
})

setMethod("getP_internal", "unmarkedFitOccuMS", function(object){
  J <- ncol(object@data@y)
  N <- nrow(object@data@y)
  pred <- predict(object, type='det', level=NULL)
  lapply(pred, function(x) matrix(x$Predicted, nrow=N, ncol=J, byrow=T))
})

setMethod("getP_internal", "unmarkedFitOccuMulti", function(object){
  ylist <- object@data@ylist
  S <- length(ylist)
  N <- nrow(ylist[[1]])
  J <- ncol(ylist[[1]])
  pred <- predict(object,'det', level=NULL)
  stopifnot(nrow(pred) == N*J)
  out <- list()
  for (i in 1:S){
    pmat <- matrix(pred[[i]]$Predicted, nrow=N, ncol=J, byrow=TRUE)
    out[[i]] <- pmat
  }
  names(out) <- names(ylist)
  out
})

setMethod("getP_internal", "unmarkedFitOccuTTD", function(object){
  N <- nrow(object@data@y)
  lam <- predict(object, 'det', na.rm=FALSE)$Predicted
  tmax <- as.numeric(t(object@data@surveyLength))
  tdist <- ifelse("shape" %in% names(object@estimates), "weibull", "exp")

  not_na <- !is.na(lam)
  est_p <- rep(NA, length(lam))
  if(tdist == "weibull"){
    k <- exp(coef(object)['k(k)'])
    est_p[not_na] <- stats::pweibull(tmax[not_na], k, 1/lam[not_na])
  } else {
    est_p[not_na] <- stats::pexp(tmax[not_na], lam[not_na])
  }

  matrix(est_p, nrow=N, byrow=TRUE)
})

setMethod("getP_internal", "unmarkedFitPCO", function(object){
  umf <- object@data
  M <- numSites(umf)
  T <- umf@numPrimary
  J <- ncol(umf@y) / T
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  p <- matrix(p, M, J*T, byrow = TRUE)
  p
})
