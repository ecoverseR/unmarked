
# data will need to be an unmarkedMultFrame
gpcount <- function(lambdaformula, phiformula, pformula, data,
    mixture=c('P', 'NB', 'ZIP'), K, starts, method = "BFGS", se = TRUE,
    engine=c('C', 'R'), threads=1, ...)
{
if(!is(data, "unmarkedFrameGPC"))
    stop("Data is not of class unmarkedFrameGPC.")
mixture <- match.arg(mixture)
engine <- match.arg(engine)

formulas <- list(lambda = lambdaformula, phi = phiformula, det = pformula)
check_no_support(formulas)

D <- getDesign(data, formulas)
ym <- D$y  # MxJT

if(missing(K) || is.null(K)) {
    K <- max(ym, na.rm=TRUE) + 100
    warning("K was not specified, so was set to max(y)+100 =", K)
}
M <- N <- 0:K
lM <- length(M)
I <- nrow(ym)
T <- data@numPrimary
if(T==1)
    stop("use pcount instead")
J <- numY(data) / T

y <- array(ym, c(I, J, T))

lamPars <- colnames(D$X_lambda)
detPars <- colnames(D$X_det)
nLP <- ncol(D$X_lambda)
nPP <- ncol(D$X_phi)
phiPars <- colnames(D$X_phi)
nDP <- ncol(D$X_det)
nP <- nLP + nPP + nDP + (mixture%in%c('NB','ZIP'))
if(!missing(starts) && length(starts) != nP)
    stop("There should be", nP, "starting values, not", length(starts))

if(identical(engine, "R")) {
# Minus negative log-likelihood
nll <- function(pars) {
    lam <- exp(D$X_lambda %*% pars[1:nLP] + D$offset_lambda)
    phi <- plogis(D$X_phi %*% pars[(nLP+1):(nLP+nPP)] + D$offset_phi)
    phi <- matrix(phi, I, T, byrow=TRUE)
    p <- plogis(D$X_det %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + D$offset_det)
    p <- matrix(p, I, byrow=TRUE)
    p <- array(p, c(I, J, T))  # byrow?
    L <- rep(NA, I)
    for(i in 1:I) {
        f <- switch(mixture,
            P = dpois(M, lam[i], log=TRUE),
            NB = dnbinom(M, mu=lam[i], size=exp(pars[nP]), log=TRUE),
            ZIP = log(dzip(M, lambda=lam[i], psi=plogis(pars[nP])))
        )
        ghi <- rep(0, lM)
        for(t in 1:T) {
            gh <- matrix(-Inf, lM, lM)
            for(m in M) {
                if(m < max(y[i,,], na.rm=TRUE)) {
                    gh[,m+1] <- -Inf
                    next
                }
                if(is.na(phi[i,t])) {
                    g <- rep(0, lM)
                    g[N>m] <- -Inf
                }
                else
                    g <- dbinom(N, m, phi[i,t], log=TRUE)
                h <- rep(0, lM)
                for(j in 1:J) {
                    if(is.na(y[i,j,t]))
                        next
                    h <- h + dbinom(y[i,j,t], N, p[i,j,t], log=TRUE)
                }
                gh[,m+1] <- g + h
            }
            ghi <- ghi + log(colSums(exp(gh))) # sum over N(t)
        }
        fgh <- f + ghi
        L[i] <- sum(exp(fgh)) # sum over M
    }
    return(-sum(log(L)))
}
} else
if(identical(engine, "C")) {
    nll <- function(pars) {
        beta_lambda <- pars[1:nLP]
        beta_phi <- pars[(nLP+1):(nLP+nPP)]
        beta_det <- pars[(nLP+nPP+1):(nLP+nPP+nDP)]
        log.alpha <- 1
        if(mixture %in% c("NB", "ZIP"))
            log.alpha <- pars[nP]
        nll_gpcount(ym, D$X_lambda, D$X_phi, D$X_det, beta_lambda, beta_phi, beta_det,
                    log.alpha, D$offset_lambda, D$offset_phi, D$offset_det,
                    as.integer(K), mixture, T, threads)
    }
}

if(missing(starts)) starts <- rep(0, nP)
fm <- optim(starts, nll, method = method, hessian = se, ...)
covMat <- invertHessian(fm, nP, se)
ests <- fm$par
fmAIC <- 2 * fm$value + 2 * nP

nbParm <- switch(mixture, P={character(0)}, NB={"alpha"}, ZIP={"psi"})

names(ests) <- c(lamPars, phiPars, detPars, nbParm)

lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = ests[1:nLP],
    covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
    invlinkGrad = "exp")

phiEstimates <- unmarkedEstimate(name = "Availability",
                                 short.name = "phi",
                                 estimates = ests[(nLP+1):(nLP+nPP)],
                                 covMat = as.matrix(covMat[(nLP+1) :
                                 (nLP+nPP), (nLP+1):(nLP+nPP)]),
                                 invlink = "logistic",
                                 invlinkGrad = "logistic.grad")

detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
    estimates = ests[(nLP+nPP+1):(nLP+nPP+nDP)],
    covMat = as.matrix(
        covMat[(nLP+nPP+1):(nLP+nPP+nDP), (nLP+nPP+1):(nLP+nPP+nDP)]),
    invlink = "logistic", invlinkGrad = "logistic.grad")

estimateList <- unmarkedEstimateList(list(lambda=lamEstimates, phi=phiEstimates, det=detEstimates))

if(identical(mixture,"NB"))
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")

if(identical(mixture,"ZIP")) {
    estimateList@estimates$psi <- unmarkedEstimate(name="Zero-inflation",
        short.name = "psi", estimates = ests[nP],
        covMat=as.matrix(covMat[nP, nP]), invlink = "logistic",
        invlinkGrad = "logistic.grad")
}

umfit <- new("unmarkedFitGPC", fitType = "gpcount",
    call = match.call(), formlist = formulas,
    data = data, estimates = estimateList, sitesRemoved = D$removed.sites,
    AIC = fmAIC, opt = fm, negLogLike = fm$value, nllFun = nll,
    mixture=mixture, K=K)

return(umfit)
}




