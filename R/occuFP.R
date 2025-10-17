#  Fit occupancy data where false positives oc cur
#  three data types are allowed
#  1) standard occupancy data as in MacKenzie et al (2002).
#  2) Royle-Link (2006) type
#  3) certain-uncertain data like Miller (2011)



occuFP <- function(detformula = ~ 1,FPformula = ~ 1,Bformula = ~ 1,stateformula = ~ 1, data, starts,
                 method = "BFGS", se = TRUE, engine = "R", ...) {


    if(!is(data, "unmarkedFrameOccuFP"))   stop("Data is not an unmarkedFrameOccuFP object.")

    type <- data@type
    if(sum(type[2:3])==0)   stop("Only type 1 data. No data types with false positives. Use occu instead.")

    formulas <- list(det = detformula, fp = FPformula, b = Bformula, state = stateformula)
    check_no_support(formulas)
    dm <- getDesign(data, formulas)
    y <- dm$y
    if(any(type[1:2]>0)) if(any(y[,1:sum(type[1:2])]>1,na.rm = TRUE))   stop("Values of y for type 1 and type 2 data must be 0 or 1.")
    if(type[3]>0) if(any(y[,1:sum(type[3])]>2,na.rm = TRUE))   stop("Values of y for type 3 data must be 0, 1, or 2.")


    J <- ncol(y)
    M <- nrow(y)

    occParms <- colnames(dm$X_state)
    detParms <- colnames(dm$X_det)
    FPParms <- colnames(dm$X_fp)
    if(type[3]!=0) BParms <- colnames(dm$X_b) else BParms = NULL
    nDP <- ncol(dm$X_det)
    nFP <- ncol(dm$X_fp)
    nBP <- ifelse(type[3]!=0,ncol(dm$X_b),0)
    nOP <- ncol(dm$X_state)

    nP <- nDP + nOP + nFP + nBP
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec0 <- as.numeric(t(y==0))
    yvec1 <- as.numeric(t(y==1))
    yvec2 <- as.numeric(t(y==2))
    navec <- c(t(is.na(y)))


        nll <- function(params) {
            psi <- plogis(dm$X_state %*% params[1 : nOP] + dm$offset_state)
            pvec <- plogis(dm$X_det %*% params[(nOP + 1) : (nOP + nDP)] + dm$offset_det)
            fvec <- plogis(dm$X_fp %*% params[(nOP + nDP + 1) : (nOP + nDP + nFP)] + dm$offset_fp)
            if (type[1]!=0) fvec[rep(c(rep(TRUE,type[1]),rep(FALSE,sum(type[2:3]))),M)] = 0
            if (type[3]!=0){
              bvec <- plogis(dm$X_b %*% params[(nOP + nDP + nFP + 1) : nP] + dm$offset_b)
              if (type[1]!=0|type[2]!=0) bvec[rep(c(rep(TRUE,sum(type[1:2])),rep(FALSE,type[3])),M)] = 0}
            if (type[3]==0){
              bvec <- matrix(0,M*J,1)
            }
            cp0 <- (1-fvec)^(yvec0) * (fvec)^(yvec1) * (1-yvec2)
            cp1 <- (1 - pvec)^(yvec0)*(pvec*(1-bvec))^(yvec1)*(pvec*bvec)^(yvec2)
            cp0[navec] <- 1 # so that NA's don't modify likelihood
            cp1[navec] <- 1 # so that NA's don't modify likelihood
            cp0mat <- matrix(cp0, M, J, byrow = TRUE) #
            cp1mat <- matrix(cp1, M, J, byrow = TRUE) #
            loglik <- log(rowProds(cp0mat) *(1-psi) + rowProds(cp1mat) * psi)
            -sum(loglik)
        }


    if(missing(starts)) starts <- rep(0, nP)
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    covMat <- invertHessian(fm, nP, se)
    ests <- fm$par
    fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)
    names(ests) <- c(occParms, detParms,FPParms,BParms)

    state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                              estimates = ests[1:nOP],
                              covMat = as.matrix(covMat[1:nOP,1:nOP]),
                              invlink = "logistic",
                              invlinkGrad = "logistic.grad")

    det <- unmarkedEstimate(name = "Detection", short.name = "p",
                            estimates = ests[(nOP + 1) : (nOP+nDP)],
                            covMat = as.matrix(covMat[(nOP + 1) : (nOP+nDP),
                                                      (nOP + 1) : (nOP+nDP)]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

    fp <- unmarkedEstimate(name = "false positive", short.name = "fp",
                            estimates = ests[(nOP+nDP+1) : (nOP+nDP+nFP)],
                            covMat = as.matrix(covMat[(nOP+nDP+1) : (nOP+nDP+nFP),
                                                      (nOP+nDP+1) : (nOP+nDP+nFP)]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")
    if (type[3]!=0){
    b <- unmarkedEstimate(name = "Pcertain", short.name = "b",
                            estimates = ests[(nOP+nDP+nFP+1): nP],
                            covMat = as.matrix(covMat[(nOP+nDP+nFP+1): nP,
                                                      (nOP+nDP+nFP+1): nP]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")
    }



    if (type[3]!=0) {
      estimateList <- unmarkedEstimateList(list(state=state, det=det,fp=fp,b=b))}
    if (type[3]==0) {
      estimateList <- unmarkedEstimateList(list(state=state, det=det,fp=fp))
    }

    umfit <- new("unmarkedFitOccuFP", fitType = "occuFP", call = match.call(),
                 formlist = formulas, type = type, data = data,
                 sitesRemoved = dm$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = fm,
                 negLogLike = fm$value,
                 nllFun = nll)

    return(umfit)
}
