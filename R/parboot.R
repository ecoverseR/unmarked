

# ----------------------- PARAMETRIC BOOTSTRAP --------------------------
setMethod("replaceY", "unmarkedFrame", function(object, newY, replNA=TRUE, ...){
    if(replNA) is.na(newY) <- is.na(object@y)
    object@y <- newY
    object
})
setMethod("replaceY", "unmarkedFrameOccuMulti",
          function(object, newY, replNA=TRUE, ...){
      if(replNA){
        newY <- mapply(function(x, y){ is.na(x) <- is.na(y); x},
                       newY , object@ylist, SIMPLIFY=FALSE)
      }
      object@ylist <- newY
      object
})


setMethod("parboot", "unmarkedFit",  function(object, statistic=SSE, nsim=10,
          report, parallel = FALSE, ncores, ...){

  if(!missing(report)){
    warning("report argument is non-functional and will be deprecated in the next version", call.=FALSE)
  }

  dots <- list(...)
  call <- match.call(call = sys.call(-1))
  stopifnot(is.function(statistic))
  starts <- as.numeric(coef(object))
  # Get rid of starting values if model was fit with TMB
  if(methods::.hasSlot(object, "TMB") && !is.null(object@TMB)) starts <- NULL

  t0 <- statistic(object, ...)

  simList <- simulate(object, nsim = nsim, na.rm = FALSE)

  availcores <- parallel::detectCores() - 1
  if(missing(ncores) || ncores > availcores) ncores <- availcores

  cl <- NULL
  if(parallel){
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(unmarked))
    env_vars <- c("dots", "replaceY")
    fm.nms <- all.names(object@call)
    if (!any(grepl("~", fm.nms))) env_vars <- c(env_vars, fm.nms[2])
    if(.hasSlot(object@data, "piFun")) env_vars <- c(env_vars, object@data@piFun)
    parallel::clusterExport(cl, env_vars, envir = environment())
    parallel::clusterEvalQ(cl, list2env(dots))
  }

  run_sim <- function(x, object, statistic, starts, t0, ...){
    simdata <- replaceY(object@data, x)
    tryCatch({
      #if(runif(1,0,1) < 0.5) stop("fail") # for testing error trapping
      fit <- update(object, data = simdata, starts = starts, se = FALSE)
      statistic(fit, ...)
    }, error=function(e){
      t0[] <- NA
      t0
    })
  }

  # Uses pbapply if available, or parSapply if not (see utils.R)
  t.star <- t(sapply2(simList, run_sim, object=object,
                              statistic=statistic, starts=starts, t0=t0,
                              cl=cl, ...))
  if(length(t0) == 1) t.star <- matrix(t.star, ncol=1)

  failed <- apply(t.star, 1, function(x) any(is.na(x)))
  if(all(failed)){
    stop("Model fitting failed in all sims.", call.=FALSE)
  }
  if(sum(failed) > 0){
    warning(paste0("Model fitting failed in ",sum(failed), " sims."), call.=FALSE)
    t.star <- t.star[!failed,,drop=FALSE]
  }

  new("parboot", call = call, t0 = t0, t.star = t.star)

})


setMethod("show", "parboot", function(object)
{
    t.star <- object@t.star
    t0 <- object@t0
    nsim <- nrow(t.star)
    biasMat <- pMat <- matrix(NA, nsim, length(t0))
    for(i in 1:nsim) {
        biasMat[i,] <- t0 - t.star[i,]
        pMat[i,] <- abs(t.star[i,] - 1) > abs(t0 - 1)
        }
    bias <- colMeans(biasMat)
    bias.se <- apply(biasMat, 2, sd)
    p.val <- colSums(pMat) / (1 + nsim)
    stats <- data.frame("t0" = t0, "mean(t0 - t_B)" = bias,
        "StdDev(t0 - t_B)" = bias.se, "Pr(t_B > t0)" = p.val,
        check.names = FALSE)
    cat("\nCall:", deparse(object@call, width.cutoff=500), fill=T)
    cat("\nParametric Bootstrap Statistics:\n")
    print(stats, digits=3)
    cat("\nt_B quantiles:\n")
    print(t(apply(t.star, 2, quantile,
        probs=c(0, 2.5, 25, 50, 75, 97.5, 100) / 100)), digits=2)
    cat("\nt0 = Original statistic computed from data\n")
    cat("t_B = Vector of bootstrap samples\n\n")
})




setMethod("plot", signature(x="parboot", y="missing"),
    function(x, y, xlab, main = "Parametric Bootstrapped Samples", xlim,
        ...)
{
    t.star <- x@t.star
    t0 <- x@t0
    if(length(t0) > 1) {
      oldask <- devAskNewPage(ask = dev.interactive(orNone = TRUE))
      on.exit(devAskNewPage(oldask))
    }
    for(i in 1:length(t0)) {
        if(missing(xlab))
            xlab <- colnames(t.star)[i]
        h <- hist(t.star[,i], plot = FALSE)
        if(missing(xlim))
            xl <- c(min(h$breaks[1], t0[i]), max(max(h$breaks), t0[i]))
        else
            xl <- xlim
        hist(t.star[,i], xlab=xlab, xlim = xl, main = main, ...)
        abline(v=t0[i], lty=2)
    }
})
