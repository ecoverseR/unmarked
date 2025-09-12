setMethod("show", "unmarkedEstimateList",
    function(object) {
      for(est in object@estimates) {
        show(est)
        cat("\n")
      }
    })

setMethod("summary", "unmarkedEstimateList",
    function(object)
{
    sumList <- list()
    for(i in 1:length(object@estimates)) {
        sumList[[i]] <- summary(object@estimates[[i]])
        cat("\n")
        }
    names(sumList) <- names(object@estimates)
    invisible(sumList)
})


setMethod("estimates", "unmarkedEstimate",
    function(object) {
      object@estimates
    })

setMethod("estimates", "unmarkedEstimateList",
		function(object) {
			object@estimates
		})

unmarkedEstimateList <- function(l) {
  new("unmarkedEstimateList", estimates = l)
}



unmarkedEstimate <- function(name, short.name, estimates = numeric(0), 
                             covMat = matrix(0, 0, 0), fixed=NULL,
                             invlink, invlinkGrad, randomVarInfo=list())
{

    if(is.null(fixed)) fixed <- 1:length(estimates)
    new("unmarkedEstimate",
        name = name,
        short.name = short.name,
        estimates = estimates,
        covMat = covMat,
        fixed = fixed,
        invlink = invlink,
        invlinkGrad = get_invlinkGrad(invlink),
        randomVarInfo = randomVarInfo)
}

get_invlinkGrad <- function(invlink){
  switch(invlink, logistic = "logistic.grad", cloglog = "cloglog.grad",
         exp = "exp", multinomial = "multinomial")
}

setMethod("show", signature(object = "unmarkedEstimate"),
    function(object)
{

    has_random <- methods::.hasSlot(object, "randomVarInfo") &&
      length(object@randomVarInfo) > 0

    fixed <- 1:length(object@estimates)
    if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
    ests <- object@estimates[fixed]
    SEs <- SE(object)[fixed]
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <- !(length(ests) == 1 |
                       identical(names(ests), "(Intercept)") |
                       identical(names(ests), "1"))

    cat(object@name,":\n", sep="")

    if(has_random){
      print_randvar_info(object@randomVarInfo)
      cat("\nFixed effects:\n")
    }

    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)

})



setMethod("summary", signature(object = "unmarkedEstimate"),
    function(object)
{
    fixed <- 1:length(object@estimates)
    if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
    ests <- object@estimates[fixed]
    SEs <- SE(object)[fixed]
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <-
        !(length(ests) == 1 | identical(names(ests), "(Intercept)") | identical(names(ests), "1"))
    invlink <- object@invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit",
                   cloglog = "cloglog")
    cat(object@name, " (", link, "-scale)", ":\n", sep="")

    has_random <- methods::.hasSlot(object, "randomVarInfo") && 
      length(object@randomVarInfo) > 0
    if(has_random){
      print_randvar_info(object@randomVarInfo)
      cat("\nFixed effects:\n")
    }

    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)
    invisible(outDF)
})



# Compute linear combinations of estimates in unmarkedEstimate objects.

setMethod("linearComb",
    signature(obj = "unmarkedEstimate", coefficients = "matrixOrVector"),
    function(obj, coefficients, offset = NULL, re.form = NULL)
{
    if(!is(coefficients, "matrix"))
        coefficients <- t(as.matrix(coefficients))
    est <- obj@estimates
    covMat <- obj@covMat
    if(!is.null(re.form) & .hasSlot(obj, "fixed")){
      est <- est[obj@fixed]
      covMat <- covMat[obj@fixed, obj@fixed, drop=FALSE]
    }
    stopifnot(ncol(coefficients) == length(est))
    if (is.null(offset))
        offset <- rep(0, nrow(coefficients))
    e <- as.vector(coefficients %*% est) + offset
    v <- coefficients %*% covMat %*% t(coefficients)
    if (!is.null(obj@covMatBS)) {
        v.bs <- coefficients %*% obj@covMatBS %*% t(coefficients)
    } else {
        v.bs <- NULL
    }
    umelc <- new("unmarkedLinComb", parentEstimate = obj,
                 estimate = e, covMat = v, covMatBS = v.bs,
                 coefficients = coefficients)
    umelc
})


# Transform an unmarkedEstimate object to it's natural scale.
#setMethod("backTransform",
#    signature(obj = "unmarkedEstimate"),
#    function(obj) {
#      stopifnot(length(obj@estimates) == 1)
#      e <- eval(call(obj@invlink,obj@estimates))
#      v <- (eval(call(obj@invlinkGrad,obj@estimates)))^2 * obj@covMat
#
#      if(is(obj, "unmarkedEstimateLinearComb")) {
#        coef <- obj@coefficients
#        orig <- obj@originalEstimate
#      } else {
#        coef <- 1
#        orig <- obj
#      }
#
#      umebt <- new("unmarkedEstimateBackTransformed",
#          name = paste(obj@name,"transformed to native scale"),
#          estimates = e, covMat = v,
#          invlink = "identity", invlinkGrad = "identity",
#          originalEstimate = orig, coefficients = coef,
#          transformation = obj@invlink)
#      umebt
#    })



# backTransform is only valid for an unmarkedEstimate of length = 1.
# can backtranform a fit directly if it has length 1
# o.w. give error
setMethod("backTransform", "unmarkedEstimate", function(obj)
{
    if(length(obj@estimates) == 1) {
        lc <- linearComb(obj, 1)
        return(backTransform(lc))
    } else {
        stop("Cannot directly back-transform an unmarkedEstimate with length > 1.\nUse linearComb() and then backTransform() the resulting scalar linear combination.")
    }
})


# Compute standard error of an unmarkedEstimate object.
setMethod("SE", signature(obj = "unmarkedEstimate"), function(obj, fixedOnly=TRUE)
{
    sqrt(diag(vcov(obj, fixedOnly=fixedOnly)))
})


setMethod("[", signature("unmarkedEstimateList"),
    function(x, i, j, drop)
{
    x@estimates[[i]]
})


setMethod("names", "unmarkedEstimateList",
    function(x)
{
    names(x@estimates)
})


setMethod("coef", "unmarkedEstimate",
    function(object, altNames = TRUE, fixedOnly=TRUE, ...)
{
    coefs <- object@estimates
    if(fixedOnly){
      fixed <- 1:length(coefs)
      if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
      coefs <- coefs[fixed]
    }
    names(coefs)[names(coefs) == "(Intercept)"] <- "Int"
    if(altNames) {
        names(coefs) <- paste(object@short.name, "(", names(coefs), ")",
                              sep="")
    }
    coefs
})


setMethod("vcov", "unmarkedEstimate",
    function(object, fixedOnly=TRUE, ...)
{
        v <- object@covMat
        if(fixedOnly){
          fixed <- 1:nrow(v)
          if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
          v <- as.matrix(v[fixed,fixed])
        }

        rownames(v) <- colnames(v) <- names(coef(object, fixedOnly=fixedOnly))
        v
})


setMethod("confint", "unmarkedEstimate",
    function(object, parm, level = 0.95)
{
    if(missing(parm)) parm <- object@fixed
    ests <- object@estimates[parm]
    ses <- SE(object)[parm]
    z <- qnorm((1-level)/2, lower.tail = FALSE)
    lower.lim <- ests - z*ses
    upper.lim <- ests + z*ses
    ci <- as.matrix(cbind(lower.lim, upper.lim))
    rownames(ci) <- names(coef(object))[parm]
    colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
    ci
})


# New estimate stuff (temporary location?)
fit_optim <- function(nll, starts, method, se, estimate_list, par_inds, ...){
  nP <- max(unlist(par_inds))
  if(is.null(starts)) starts <- rep(0, nP)
  if(length(starts) != nP){
    stop(paste("The number of starting values should be", nP))
  }
  names(starts) <- unlist(lapply(par_inds, names))
  fm <- optim(starts, nll, method = method, hessian = se, ...)

  estimate_list <- insert_estimates(estimate_list, fm, par_inds)
  
  list(opt = fm, TMB = NULL, nll = nll,
       AIC = 2 * fm$value + 2 * nP,
       estimate_list = estimate_list)
}

get_parameter_inds <- function(estimate_list, design_mats){
  types <- names(estimate_list)
  inds <- list()
  idx <- 0
  for (i in types){
    mat_names <- names(design_mats)[grepl("X_", names(design_mats))]
    mat_names <- gsub("X_", "", mat_names)
    if(!i %in% mat_names){
      # Scalar parameter
      npar <- 1
      idx_names <- "(Intercept)"
    } else {
      npar <- ncol(design_mats[[paste0("X_", i)]])
      idx_names <- colnames(design_mats[[paste0("X_", i)]])
    }
    inds[[i]] <- 1:npar + idx
    names(inds[[i]]) <- idx_names
    idx <- max(inds[[i]])
  }
  inds
}

insert_estimates <- function(estimate_list, opt, par_inds){
  covMat <- invertHessian(opt, max(unlist(par_inds)), !is.null(opt$hessian))
  for (i in names(estimate_list)){
    idx <- par_inds[[i]]
    est <- opt$par[idx]
    estimate_list@estimates[[i]]@estimates <- est
    estimate_list@estimates[[i]]@covMat <- covMat[idx, idx, drop = FALSE]
    estimate_list@estimates[[i]]@fixed <- 1:length(est)
  }
  estimate_list
}

insert_TMB_estimates <- function(estimate_list, tmb_out, par_inds, formulas, umf){ 
  for (i in names(estimate_list)){
    coefs <- get_coef_info(tmb_out$sdr, i, names(par_inds[[i]]), par_inds[[i]])

    if(i %in% names(formulas)){   
      rand <- get_randvar_info(tmb_out$sdr, i, formulas[[i]],
                               get_covariates(umf, i))
      estimate_list@estimates[[i]]@randomVarInfo <- rand
    }

    estimate_list@estimates[[i]]@estimates <- coefs$est
    estimate_list@estimates[[i]]@covMat <- coefs$cov
    estimate_list@estimates[[i]]@fixed <- 1:length(par_inds[[i]])

  }

  estimate_list
}

# Temporary wrapper until we can replace fit_TMB entirely
fit_TMB2 <- function(model, starts, method, estimate_list, par_inds, 
                     tmb_inputs, umf, ...){

  tmb_out <- fit_TMB(model, tmb_inputs$data, tmb_inputs$pars, tmb_inputs$rand_ef,
                     starts=starts, method, ...)
  
  estimate_list <- insert_TMB_estimates(estimate_list, tmb_out, par_inds,
                                        tmb_inputs$formulas, umf)
  
  list(opt = tmb_out$opt, TMB = tmb_out$TMB, nll = tmb_out$TMB$fn,
       AIC = tmb_out$AIC,
       estimate_list = estimate_list)
}

get_TMB_inputs <- function(formulas, dm, par_inds, umf, ...){

  datalist <- lapply(names(formulas), function(x) get_covariates(umf, x))
  names(datalist) <- names(formulas)
  dms <- dm[paste0("X_", names(formulas))] 
  Zs <- dm[paste0("Z_", names(formulas))] 

  mods <- names(formulas)
  ngv <- mapply(get_group_vars, formulas, datalist)
  names(ngv) <- paste0("n_group_vars_",mods)
  ngroup <- mapply(get_nrandom, formulas, datalist, SIMPLIFY=FALSE)
  names(ngroup) <- paste0("n_grouplevels_",mods)
  names(dms) <- paste0("X_", mods)
  names(Zs) <- paste0("Z_", mods)

  dat <- c(list(y = dm$y), ngv, ngroup, dms, Zs, list(...))

  if(any(grepl("offset", names(dm)))){
    dat <- c(dat, dm[paste0("offset_", mods)])
  }

  beta <- lapply(par_inds, function(x) rep(0, length(x)))
  names(beta) <- paste0("beta_", names(par_inds))
  b <- lapply(ngroup, function(x) rep(0, sum(x)))
  names(b) <- paste0("b_", mods)
  lsigma <- lapply(ngv, function(x) rep(0, x))
  names(lsigma) <- paste0("lsigma_", mods)

  pars <- c(beta, b, lsigma)

  rand_ef <- paste0(names(b))[sapply(formulas, has_random)]
  if(length(rand_ef) == 0) rand_ef <- NULL

  list(data=dat, pars=pars, rand_ef=rand_ef, formulas = formulas)
}
