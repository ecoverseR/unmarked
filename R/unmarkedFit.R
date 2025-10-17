# -------------------------- Show and Summary ----------------------------


setMethod("show", "unmarkedFit", function(object)
{
  summary(object)
})


setMethod("summary", "unmarkedFit", function(object)
{
  summary_internal(object)
})


setMethod("summary_internal", "unmarkedFit", function(object)
{
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  summaryOut <- summary(object@estimates)
  cat("AIC:", object@AIC,"\n")
  cat("Number of sites:", sampleSize(object))
  if(length(object@sitesRemoved) > 0){
    cat("\nID of sites removed due to NA:", object@sitesRemoved)
  }
  if(!identical(object@opt$convergence, 0L)){
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.", call.=FALSE)
  }

  # Check for potentially bad estimates
  if(!is.null(object@opt$hessian)){
    se <- SE(object)
    has_na <- any(is.na(se)) | any(is.nan(se))
    big_se <- any(abs(se) >= 5)
    if(has_na | big_se){
      warning("Large or missing SE values. Be very cautious using these results.", call.=FALSE)
    }
  }

  nboot <- length(object@bootstrapSamples)
  if(nboot > 0){
    cat("\nBootstrap iterations:", nboot)
  }
  cat("\n\n")
  invisible(summaryOut)
})



setMethod("summary_internal", "unmarkedFitDS", function(object)
{
  out <- callNextMethod()
  cat("Survey design: ", object@data@survey, "-transect", sep="")
  cat("\nDetection function:", object@keyfun)
  cat("\nUnitsIn:", object@data@unitsIn)
  cat("\nUnitsOut:", object@unitsOut, "\n\n")
  invisible(out)
})




# Compute linear combinations of estimates in unmarkedFit objects.
setMethod("linearComb",
    signature(obj = "unmarkedFit", coefficients = "matrixOrVector"),
    function(obj, coefficients, type, offset = NULL, re.form=NULL)
{
    stopifnot(!missing(type))
    stopifnot(type %in% names(obj))
    estimate <- obj@estimates[type]
    linearComb(estimate, coefficients, offset, re.form)
})


setMethod("backTransform", "unmarkedFit", function(obj, type)
{
    est <- obj[type]
    if(length(est@estimates) == 1) {
        lc <- linearComb(est, 1)
        return(backTransform(lc))
    } else {
        stop('Cannot directly backTransform an unmarkedEstimate with length > 1.')
        }
})


setMethod("[", "unmarkedFit",
          function(x, i, j, drop) {
              x@estimates[i]
          })


setMethod("names", "unmarkedFit",
          function(x) {
              names(x@estimates)
          })



# ---------------------- coef, vcov, and SE ------------------------------


setMethod("coef", "unmarkedFit",
    function(object, type, altNames = TRUE, fixedOnly=TRUE)
{
    if(missing(type)) {
        co <- lapply(object@estimates@estimates,
            function(x) coef(x, altNames=altNames, fixedOnly=fixedOnly))
        names(co) <- NULL
        co <- unlist(co)
    } else {
        co <- coef(object[type], altNames=altNames, fixedOnly=fixedOnly)
        }
    co
})


setMethod("vcov", "unmarkedFit",
    function (object, type, altNames = TRUE, method = "hessian", fixedOnly=TRUE, ...)
{
    method <- match.arg(method, c("hessian", "nonparboot"))
    if(.hasSlot(object, "TMB") && !is.null(object@TMB)) method <- "TMB"
    switch(method,
           hessian = {
            if (is.null(object@opt$hessian)) {
                stop("Hessian was not computed for this model.")
            }
            v <- solve(hessian(object))
        },
        nonparboot = {
            if (is.null(object@bootstrapSamples)) {
                stop("No bootstrap samples have been drawn. Use nonparboot first.")
                }
            v <- object@covMatBS
        },
        TMB = {
          return(vcov_TMB(object, type, fixedOnly))
        })
    rownames(v) <- colnames(v) <- names(coef(object, altNames=altNames))
    if (missing(type)) {
        return (v)
    } else {
      object[type]
        inds <- .estimateInds(object)[[type]]
        return (v[inds, inds, drop = FALSE])
        }
})

## A helper function to return a list of indices for each estimate type
##
.estimateInds <- function(umf) {
  ## get length of each estimate
  estimateLengths <- sapply(umf@estimates@estimates, function(est) {
    length(coef(est))
  })
  ## recurse function to generate list of indices
  estimateInds <- function(type) {
    if(type==1) {
      return(list(seq(length=estimateLengths[1])))
    } else {
      prev.list <- estimateInds(type-1)
      prev.max <- max(prev.list[[type-1]])
      return(c(prev.list, list(seq(prev.max+1, prev.max +
                                   estimateLengths[type]))))
    }
  }
  retlist <- estimateInds(length(estimateLengths))
  names(retlist) <- names(umf)
  retlist
}

setMethod("vcov", "unmarkedFitOccuMulti",
    function (object, type, altNames = TRUE, method = "hessian", ...)
{

  pen <- object@call[["penalty"]]
  if(is.null(pen)) pen <- 0
  if(pen>0) method <- "nonparboot"
  callNextMethod(object, type, altNames, method=method, ...)

})


setMethod("SE", "unmarkedFit", function(obj,...)
{
    v <- vcov(obj,...)
    sqrt(diag(v))
})



setMethod("logLik", "unmarkedFit", function(object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    ll <- -object@negLogLike
    #attr(ll, "df") <- length(coef(object))
    #class(ll) <- "logLik"
    return(ll)
})



setMethod("LRT", c(m1="unmarkedFit", m2="unmarkedFit"), function(m1, m2)
{
    ll1 <- logLik(m1)
    ll2 <- logLik(m2)
    chisq <- 2 * abs(ll1 - ll2)
    DF <- abs(length(coef(m1)) - length(coef(m2)))
    pval <- pchisq(chisq, DF, lower.tail=FALSE)
    return(data.frame(Chisq=chisq, DF = DF, 'Pr(>Chisq)' = pval,
        check.names=F))
})





setMethod("confint", "unmarkedFit", function(object, parm, level = 0.95,
    type, method = c("normal", "profile"))
{
    method <- match.arg(method)
    if(missing(type))
        stop(paste("Must specify type as one of (", paste(names(object), collapse=", "),").",sep=""))
    if(missing(parm))
        parm <- object[type]@fixed
    if(method == "normal") {
        callGeneric(object[type],parm = parm, level = level)
    } else {
        nllFun <- nllFun(object)
        ests <- mle(object)
        nP <- length(parm)
        ci <- matrix(NA, nP, 2)

        ## create table to match parm numbers with est/parm numbers
        types <- names(object)
        numbertable <- data.frame(type = character(0), num = numeric(0))
        for(i in seq(length=length(types))) {
            length.est <- length(object[i]@estimates)
            numbertable <- rbind(numbertable, data.frame(type =
                rep(types[i], length.est), num = seq(length=length.est)))
            }
        parm.fullnums <- which(numbertable$type == type &
            numbertable$num %in% parm)

        for(i in seq(length=nP)) {
            cat("Profiling parameter",i,"of",nP,"...")
            se <- SE(object[type])
            whichPar <- parm.fullnums[i]
            ci[i,] <- profileCI(nllFun, whichPar=whichPar, MLE=ests,
                interval=ests[whichPar] + 10*se[i]*c(-1,1), level=level)
            cat(" done.\n")
            }
        rownames(ci) <- names(coef(object[type]))[parm]
        colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
        return(ci)
        }
})


setMethod("profile", "unmarkedFit",
    function(fitted, type, parm, seq)
{
    stopifnot(length(parm) == 1)
    MLE <- mle(fitted)
    SE(fitted[type])
    nPar <- length(mle(fitted))
    nll <- nllFun(fitted)

    ## create table to match parm numbers with est/parm numbers
    types <- names(fitted)
    numbertable <- data.frame(type = character(0), num = numeric(0))
    for(i in seq(length=length(types))) {
        length.est <- length(fitted[i]@estimates)
        numbertable <- rbind(numbertable,
        data.frame(type = rep(types[i], length.est),
            num = seq(length=length.est)))
        }
    parm.fullnums <- which(numbertable$type == type &
        numbertable$num == parm)

    f <- function(value) {
        fixedNLL <- genFixedNLL(nll, parm.fullnums, value)
        mleRestricted <- optim(rep(0,nPar), fixedNLL)$value
        mleRestricted
        }
        prof.out <- sapply(seq, f)
        prof.out <- cbind(seq, prof.out)
        new("profile", prof = prof.out)
})



setMethod("hessian", "unmarkedFit",
    function(object)
{
    object@opt$hessian
})


setMethod("sampleSize", "unmarkedFit", function(object) {
    data <- getData(object)
    M <- numSites(data)
    M <- M - length(object@sitesRemoved)
    M
})


setMethod("getData", "unmarkedFit",function(object) {
    object@data
})


setMethod("nllFun", "unmarkedFit", function(object) object@nllFun)


setMethod("mle", "unmarkedFit", function(object) object@opt$par)

setMethod("smoothed","unmarkedFitColExt",
    function(object, mean) {
        if(mean) object@smoothed.mean
        else object@smoothed
    })


setMethod("projected","unmarkedFitColExt", function(object, mean) {
    if(mean) object@projected.mean
    else object@projected
})

setMethod("plot", c("profile", "missing"), function(x) {
    plot(x@prof[,1], x@prof[,2], type = "l")
})


setMethod("plot", c(x = "unmarkedFit", y = "missing"), function(x, y, ...)
{
  residual_plot(x, ...)
})

setMethod("residual_plot", "unmarkedFit", function(x, ...)
{
    r <- residuals(x)
    e <- fitted(x, na.rm = FALSE)
    plot(e, r, ylab = "Residuals", xlab = "Predicted values")
    abline(h = 0, lty = 3, col = "gray")
})

setMethod("residual_plot", "unmarkedFitOccuMulti", function(x, ...)
{
  r <- do.call(rbind,residuals(x))
  e <- do.call(rbind,fitted(x))
  plot(e, r, ylab = "Residuals", xlab = "Predicted values")
  abline(h = 0, lty = 3, col = "gray")
})

setMethod("hist", "unmarkedFitDS", function(x, lwd=1, lty=1, ...) {
    ymat <- getY(getData(x))
    dbreaks <- getData(x)@dist.breaks
    nb <- length(dbreaks)
    mids <- (dbreaks[-1] - dbreaks[-nb]) / 2 + dbreaks[-nb]
    distances <- rep(mids, times=colSums(ymat))
    h <- hist(distances, plot=F, breaks=dbreaks)
    key <- x@keyfun
    survey <- x@data@survey
    switch(key,
           halfnorm = {
               sigma <- exp(coef(x, type="det"))
               if(length(sigma) > 1)
               stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- 2 * integrate(dnorm, dbreaks[1],
                                               dbreaks[nb], sd=sigma)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) 2 * dnorm(x, mean=0, sd=sigma),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drhn, dbreaks[1], dbreaks[nb],
                                           sigma=sigma)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drhn(r, sigma=sigma),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      })
           },
           exp = {		# This doesn't work on example fm4
               rate <- exp(coef(x, type="det"))
               if(length(rate) > 1)
                   stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- integrate(dxexp, dbreaks[1], dbreaks[nb],
                                           rate=rate)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) dxexp(x, rate=rate),
                               min(dbreaks),
                               max(dbreaks), add=T, lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drexp, dbreaks[1], dbreaks[nb],
                                           rate=rate)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drexp(r, rate=rate),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      })
           },
           hazard = {
               shape <- exp(coef(x, type="det"))
               scale <- exp(coef(x, type="scale"))
               if(length(shape) > 1)
                   stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- integrate(dxhaz, dbreaks[1], dbreaks[nb],
                                           shape=shape, scale=scale)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) dxhaz(x, shape=shape,
                                                 scale=scale),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drhaz, dbreaks[1], dbreaks[nb],
                                           shape=shape, scale=scale)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drhaz(r, shape=shape,
                                                 scale=scale),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      })
           },
           uniform = {
               switch(survey,
                      line = {
                          plot(h, freq=F, ...)
                          abline(h=1/max(dbreaks), lwd=lwd, lty=lty)
                      },
                      point = {
                          plot(h, freq=F, ...)
                          plot(function(r) (pi*r^2) / (pi*max(dbreaks)^2),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      }
                      )}
           )
})




# ----------------------- CHILD CLASS METHODS ---------------------------

# Extract detection probs
setMethod("getFP", "unmarkedFitOccuFP", function(object, na.rm = TRUE)
{
  designMats <- getDesign(getData(object), object@formlist, na.rm = na.rm)
  type = object@type
  y <- designMats$y
  M <- nrow(y)
  J <- ncol(y)
  fpars <- coef(object, type = "fp")
  f <- plogis(designMats$X_fp %*% fpars + designMats$offset_fp)
  f <- matrix(f, M, J, byrow = TRUE)
  if (type[1]!=0){
    f[,1:type[1]] = 0
  }
  return(f)
})

setMethod("getB", "unmarkedFitOccuFP", function(object, na.rm = TRUE)
{
  designMats <- getDesign(getData(object), object@formlist, na.rm = na.rm)
  y <- designMats$y
  M <- nrow(y)
  J <- ncol(y)
  type = object@type
  if (type[3]!=0){
    bpars <- coef(object, type = "b")
  b <- plogis(designMats$X_b %*% bpars + designMats$offset_b)
  b <- matrix(b, M, J, byrow = TRUE)
  }
  if (type[3]==0){
    b <- matrix(0, M, J)
  }
  return(b)
})

#Y extractors for unmarkedFit objects
setMethod("getY", "unmarkedFit", function(object) getY_internal(object))

setMethod("getY_internal", "unmarkedFit", function(object) object@data@y)
setMethod("getY_internal", "unmarkedFitOccu", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY_internal", "unmarkedFitOccuRN", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY_internal", "unmarkedFitColExt", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY_internal", "unmarkedFitOccuMulti", function(object) {
            object@data@ylist
})

# Extract formula(s) for a particular submodel
setMethod("get_formula", "unmarkedFit", function(object, type, ...){
  if(methods::.hasSlot(object, "formlist")){
    formlist <- object@formlist
  } else {
    # For partial backwards compatibility
    if(methods::.hasSlot(object, "formula") & type %in% c("state", "det")){
      formlist <- split_formula(object@formula)
      names(formlist) <- c("det", "state")
    } else {
      stop("Fitted model is incompatible with this version of unmarked. Try re-fitting the model.")  
    }
  }

  formlist[[type]]
})

# Extract the covariates used for a particular submodel
setMethod("get_covariates", "unmarkedFit", function(object, type, ...){
  get_covariates(getData(object), type, ...)
})
