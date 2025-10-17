unmarkedFrameGOccu <- function(y, siteCovs=NULL, obsCovs=NULL, numPrimary,
                             yearlySiteCovs=NULL) {
  y[y > 1] <- 1
  if(numPrimary < 2) stop("numPrimary < 2, use occu instead")
  umf <- unmarkedFrameGPC(y, siteCovs=siteCovs, obsCovs=obsCovs, 
                          numPrimary=numPrimary, yearlySiteCovs=yearlySiteCovs)
  class(umf) <- "unmarkedFrameGOccu"
  umf
}

goccu <- function(psiformula, phiformula, pformula, data,
                  linkPsi = c("logit", "cloglog"), starts, method = "BFGS",
                  se = TRUE, ...){

  linkPsi <- match.arg(linkPsi, c("logit","cloglog"))
  psiLinkFunc <- ifelse(linkPsi=="cloglog", cloglog, plogis)
  psiInvLink <- ifelse(linkPsi=="cloglog", "cloglog", "logistic")
  psiLinkGrad <- ifelse(linkPsi=="cloglog", "cloglog.grad", "logistic.grad")

  formulas <- list(psi=psiformula, phi=phiformula, det=pformula)

  data@y[data@y > 1] <- 1
 
  class(data) <- "unmarkedFrameGOccu"

  gd <- getDesign(data, formulas)
  y <- gd$y

  M <- nrow(y)
  T <- data@numPrimary
  J <- ncol(y) / T

  # Identify entirely missing primary periods at each site
  y_array <- array(t(y), c(J, T, M)) 
  missing_session <- t(apply(y_array, c(2,3), 
                           function(x) as.numeric(all(is.na(x)))))

  # Create possible states in each T
  alpha_potential <- as.matrix(expand.grid(rep(list(c(0, 1)), T)))
  n_possible <- nrow(alpha_potential)

  # Known present at each site
  known_present <- rep(0, M)
  # Known available at each site and T
  known_available <- matrix(0, nrow=M, ncol=T)
  
  for (i in 1:M){
    for (t in 1:T){
      for (j in 1:J){
        if(is.na(y_array[j,t,i])) next
        if(y_array[j, t, i] == 1){
          known_present[i] <- 1
          known_available[i,t] <- 1
        }
      }
    }
  }

  # Determine which configurations of available states should not be
  # included in likelihood because relevant primary periods were missing
  alpha_drop <- matrix(NA, M, nrow(alpha_potential))
  for (i in 1:M){
    dropped <- rep(0, nrow(alpha_potential))
    for (j in 1:nrow(alpha_potential)){
      check_drop <- alpha_potential[j,] * missing_session[i,]
      if(sum(check_drop) > 0) dropped[j] <- 1
    }
    alpha_drop[i,] <- dropped
  }

  # Bundle data for TMB
  dataList <- list(y=y, T=T, link=ifelse(linkPsi=='cloglog', 1, 0), 
                   X_psi=gd$X_psi, X_phi=gd$X_phi, X_det=gd$X_det,
                   n_possible=n_possible,
                   alpha_potential=alpha_potential, alpha_drop = alpha_drop,
                   known_present=known_present, known_available=known_available, 
                   missing_session=missing_session)

  # Provide dimensions and starting values for parameters
  # This part should change to be more like occu() if we add random effects
  psi_ind <- 1:ncol(gd$X_psi)
  phi_ind <- 1:ncol(gd$X_phi) + max(psi_ind)
  det_ind <- 1:ncol(gd$X_det) + max(phi_ind)
  nP <- max(det_ind)
  params <- list(beta_psi = rep(0, length(psi_ind)), 
                 beta_phi = rep(0, length(phi_ind)), 
                 beta_det = rep(0, length(det_ind)))

  # Create TMB object
  tmb_mod <- TMB::MakeADFun(data = c(model = "tmb_goccu", dataList),
                            parameters = params,
                            DLL = "unmarked_TMBExports", silent = TRUE)

  # Optimize TMB object, print and save results
  if(missing(starts) || is.null(starts)) starts <- tmb_mod$par
  opt <- optim(starts, fn = tmb_mod$fn, gr = tmb_mod$gr, method = method,
               hessian = se, ...)

  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  names(ests) <- c(colnames(gd$X_psi), colnames(gd$X_phi), colnames(gd$X_det))
  fmAIC <- 2 * opt$value + 2 * nP


  psi_est <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                              estimates = ests[psi_ind],
                              covMat = covMat[psi_ind, psi_ind, drop=FALSE],
                              fixed = 1:ncol(gd$X_psi),
                              invlink = psiInvLink,
                              invlinkGrad = psiLinkGrad,
                              randomVarInfo=list()
                            )

  phi_est <- unmarkedEstimate(name = "Availability", short.name = "phi",
                              estimates = ests[phi_ind],
                              covMat = covMat[phi_ind, phi_ind, drop=FALSE],
                              fixed = 1:ncol(gd$X_phi),
                              invlink = "logistic",
                              invlinkGrad = "logistic.grad",
                              randomVarInfo=list()
                            )

  det_est <- unmarkedEstimate(name = "Detection", short.name = "p",
                              estimates = ests[det_ind],
                              covMat = covMat[det_ind, det_ind, drop=FALSE],
                              fixed = 1:ncol(gd$X_det),
                              invlink = "logistic",
                              invlinkGrad = "logistic.grad",
                              randomVarInfo=list()
                            )

  estimate_list <- unmarkedEstimateList(list(psi=psi_est, phi=phi_est, det=det_est))

  # Create unmarkedFit object--------------------------------------------------
  umfit <- new("unmarkedFitGOccu", fitType = "goccu", call = match.call(),
                 formlist=formulas, data = data,
                 sitesRemoved = gd$removed.sites,
                 estimates = estimate_list, AIC = fmAIC, opt = opt,
                 negLogLike = opt$value,
                 nllFun = tmb_mod$fn, TMB=tmb_mod)

  return(umfit)

}

# Methods

setMethod("getP_internal", "unmarkedFitGOccu", function(object){
  M <- numSites(object@data)
  J <- ncol(object@data@y)
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  p <- matrix(p, nrow=M, ncol=J, byrow=TRUE)
  p
})

# based on ranef for GPC
setMethod("ranef_internal", "unmarkedFitGOccu", function(object, ...){

  M <- numSites(object@data)
  JT <- obsNum(object@data)
  T <- object@data@numPrimary
  J <- JT / T

  gd <- getDesign(object@data, object@formlist, na.rm=FALSE)
  y_array <- array(t(gd$y), c(J, T, M))

  psi <- drop(plogis(gd$X_psi %*% coef(object, "psi")))
  phi <- drop(plogis(gd$X_phi %*% coef(object, "phi")))
  phi <- matrix(phi, nrow=M, ncol=T, byrow=TRUE)
  p <- getP(object)
  p_array <- array(t(p), c(J, T, M))
  
  Z <- ZZ <- 0:1
  post <- array(NA, c(M, 2, 1))
  colnames(post) <- Z

  for(i in 1:M) {
    if(i %in% object@sitesRemoved) next
    f <- dbinom(Z, 1, psi[i])
    
    ghi <- rep(0, 2)

    for(t in 1:T) {
      gh <- matrix(-Inf, 2, 2)
      for(z in Z) {
        if(z < max(y_array[,,i], na.rm=TRUE)){
          gh[,z+1] <- -Inf
          next
        }
        if(is.na(phi[i,t])) {
          g <- rep(0, 2)
          g[ZZ>z] <- -Inf
        } else{
          g <- dbinom(ZZ, z, phi[i,t], log=TRUE)
        }
        h <- rep(0, 2)
        for(j in 1:J) {
          if(is.na(y_array[j,t,i]) | is.na(p_array[j,t,i])) next
          h <- h + dbinom(y_array[j,t,i], ZZ, p_array[j,t,i], log=TRUE)
        }
        gh[,z+1] <- g + h
      }
      ghi <- ghi + log(colSums(exp(gh)))
    }
    fgh <- exp(f + ghi)
    prM <- fgh/sum(fgh)
    post[i,,1] <- prM
  }

  new("unmarkedRanef", post=post)
})


setMethod("simulate_internal", "unmarkedFitGOccu", 
          function(object, nsim){
 
  y <- object@data@y
  M <- nrow(y)
  T <- object@data@numPrimary
  JT <- ncol(y)
  J <- JT / T
  y_array <- array(t(y), c(J, T, M)) 

  psi <- predict(object, type = "psi", level = NULL, na.rm=FALSE)$Predicted
  phi <- predict(object, type = "phi", level = NULL, na.rm=FALSE)$Predicted
  phi <- matrix(phi, nrow=M, ncol=T, byrow=TRUE)
  p <- getP(object)

  sim_list <- list()

  for (i in 1:nsim){
    z <- suppressWarnings(rbinom(M, 1, psi))
    z <- matrix(z, nrow=M, ncol=T) 
    
    zz <- suppressWarnings(rbinom(M*T, 1, phi*z))
    zz <- matrix(zz, M, T)
    
    colrep <- rep(1:T, each=J)
    zz <- zz[,colrep]

    y <- suppressWarnings(rbinom(M*T*J, 1, zz*p))
    y <- matrix(y, M, JT)
    sim_list[[i]] <- y
  }

  return(sim_list)
})

setMethod("get_fitting_function", "unmarkedFrameGOccu",
          function(object, model, ...){
  goccu
})

setMethod("rebuild_call", "unmarkedFitGOccu", function(object){           
  cl <- object@call
  cl[["data"]] <- quote(object@data)
  cl[["psiformula"]] <- object@formlist$psi
  cl[["phiformula"]] <- object@formlist$phi
  cl[["pformula"]] <- object@formlist$det
  cl
})
