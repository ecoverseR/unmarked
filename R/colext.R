colext <- function(psiformula = ~ 1, gammaformula = ~ 1,
                    epsilonformula = ~ 1, pformula = ~ 1,
                    data, starts, method = "BFGS", se = TRUE, ...){

  ## truncate to 1
  data@y <- truncateToBinary(data@y)
 
  formulas <- list(psi = psiformula, col = gammaformula, ext = epsilonformula,
                   det = pformula)
  check_no_support(formulas)
  dm <- getDesign(data, formulas)

  X_col <- dm$X_col; X_ext <- dm$X_ext
  y <- dm$y
  M <- nrow(y)
  T <- data@numPrimary
  J <- ncol(y) / T
  psiParms <- colnames(dm$X_psi)
  colParms <- colnames(X_col)
  extParms <- colnames(X_ext)
  detParms <- colnames(dm$X_det)

  ## remove final year from transition prob design matrices
  X_col <- as.matrix(X_col[-seq(T,M*T,by=T),])
  X_ext <- as.matrix(X_ext[-seq(T,M*T,by=T),])

  # Determine which periods were sampled at each site
  site_sampled <- matrix(1, M, T)
  Trep <- rep(1:T, each = J)
  for (t in 1:T){
    ind <- which(Trep == t)
    for (i in 1:M){
      if(all(is.na(y[i, ind]))) site_sampled[i,t] <- 0
    }
  }

  # Determine site-periods with no detections
  no_detects <- matrix(1, M, T)
  for (i in 1:M){
    for (t in 1:T){
      ind <- which(Trep == t)
      ysub <- y[i, ind]
      if(all(is.na(ysub))) next
      ysub <- na.omit(ysub)
      if(all(ysub == 0)){
        no_detects[i,t] <- 1
      } else {
        no_detects[i,t] <- 0
      }
    }
  }

  # Parameter indices
  pind_mat <- matrix(0, 4, 2)
  pind_mat[1,] <- c(1, length(psiParms))
  pind_mat[2,] <- max(pind_mat[1,]) + c(1, length(colParms))
  pind_mat[3,] <- max(pind_mat[2,]) + c(1, length(extParms))
  pind_mat[4,] <- max(pind_mat[3,]) + c(1, length(detParms))
  
  tmb_dat <- list(y = as.vector(t(y)), X_psi = dm$X_psi, X_col = X_col, 
                  X_ext = X_ext, X_det = dm$X_det,
                  M = M, T = T, J = J, 
                  site_sampled = site_sampled, nd = no_detects)

  tmb_pars <- list(beta_psi = rep(0, length(psiParms)), 
                   beta_col = rep(0, length(colParms)),
                   beta_ext = rep(0, length(extParms)), 
                   beta_det = rep(0, length(detParms)))

  tmb_obj <- TMB::MakeADFun(data = c(model = "tmb_colext", tmb_dat), 
                            parameters = tmb_pars,
                            DLL = "unmarked_TMBExports", silent=TRUE)
  
  opt <- optim(unlist(tmb_pars), fn=tmb_obj$fn, gr=tmb_obj$gr, 
               method=method, hessian = se, ...)

  fmAIC <- 2 * opt$value + 2 * length(unlist(tmb_pars))

  sdr <- TMB::sdreport(tmb_obj)

  psi_coef <- get_coef_info(sdr, "psi", psiParms, pind_mat[1,1]:pind_mat[1,2])
  col_coef <- get_coef_info(sdr, "col", colParms, pind_mat[2,1]:pind_mat[2,2])
  ext_coef <- get_coef_info(sdr, "ext", extParms, pind_mat[3,1]:pind_mat[3,2])
  det_coef <- get_coef_info(sdr, "det", detParms, pind_mat[4,1]:pind_mat[4,2])
  
  psi <- unmarkedEstimate(name = "Initial", short.name = "psi",
                          estimates = psi_coef$ests,
                          covMat = psi_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  col <- unmarkedEstimate(name = "Colonization", short.name = "col",
                          estimates = col_coef$ests,
                          covMat = col_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  ext <- unmarkedEstimate(name = "Extinction", short.name = "ext",
                          estimates = ext_coef$ests,
                          covMat = ext_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates = det_coef$ests,
                          covMat = det_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(psi = psi, col = col,
                                            ext = ext, det = det))
  
  psis <- plogis(dm$X_psi %*% psi_coef$ests)

  # Compute projected estimates
  phis <- array(NA,c(2,2,T-1,M))
  phis[,1,,] <- plogis(X_col %x% c(-1,1) %*% col_coef$ests)
  phis[,2,,] <- plogis(X_ext %x% c(-1,1) %*% -ext_coef$ests)

  projected <- array(NA, c(2, T, M))
  projected[1,1,] <- 1 - psis
  projected[2,1,] <- psis
  for(i in 1:M) {
    for(t in 2:T) {
      projected[,t,i] <- phis[,,t-1,i] %*% projected[,t-1,i]
    }
  }
  projected.mean <- apply(projected, 1:2, mean)
  rownames(projected.mean) <- c("unoccupied","occupied")
  colnames(projected.mean) <- 1:T

  # Compute smoothed estimates
  smoothed <- calculate_smooth(y = y, psi = psis,
                               col = plogis(X_col %*% col_coef$ests),
                               ext = plogis(X_ext %*% ext_coef$ests),
                               p = plogis(dm$X_det %*% det_coef$ests),
                               M = M, T = T, J = J)
  smoothed.mean <- apply(smoothed, 1:2, mean)
  rownames(smoothed.mean) <- c("unoccupied","occupied")
  colnames(smoothed.mean) <- 1:T

  umfit <- new("unmarkedFitColExt", fitType = "colext",
                call = match.call(),
                formlist = formulas,
                data = data, sitesRemoved = dm$removed.sites,
                estimates = estimateList,
                AIC = fmAIC, opt = opt, negLogLike = opt$value,
                nllFun = tmb_obj$fn,
                projected = projected,
                projected.mean = projected.mean,
                smoothed = smoothed, smoothed.mean = smoothed.mean)

  return(umfit)
}

# Based on Weir, Fiske, Royle 2009 "TRENDS IN ANURAN OCCUPANCY"
# Appendix 1
calculate_smooth <- function(y, psi, col, ext, p, M, T, J){

  smoothed <- array(NA, c(2, T, M))

  # Turn parameters into matrices
  p <- matrix(p, M, T*J, byrow=TRUE)
  col <- matrix(col, M, (T-1), byrow=TRUE)
  ext <- matrix(ext, M, (T-1), byrow = TRUE)

  tind <- rep(1:T, each = J)

  for (i in 1:M){

    # Forward pass
    alpha1 <- matrix(NA, M, T)
    alpha0 <- matrix(NA, M, T)

    # Initialize at t=1
    ysub <- y[i, tind == 1]
    #no_detects <- all(na.omit(ysub) == 0) * 1 
    psub <- p[i, tind == 1]

    if(all(is.na(ysub))){
      # Don't include detection likelihood in calculation if no data
      alpha1[i,1] <- psi[i]
      alpha0[i,1] <- (1-psi[i])
    } else {
      # Case when z = 1
      cp <- prod(na.omit(dbinom(ysub, 1, psub)))
      alpha1[i,1] <- psi[i] * cp 

      # Case when z = 0
      cp <- prod(na.omit(dbinom(ysub, 1, 0)))
      alpha0[i,1] <- (1-psi[i]) * cp
    }

    for (t in 2:T){
      ysub <- y[i, tind == t]
      psub <- p[i, tind == t]

      if(all(is.na(ysub))){
        alpha1[i,t] <- alpha0[i,t-1] * col[i,t-1] + alpha1[i,t-1] * (1 - ext[i,t-1])
        alpha0[i,t] <- alpha0[i,t-1] * (1-col[i,t-1]) + alpha1[i,t-1] * ext[i,t-1]
      } else {
        # Case when z = 1
        cp <- prod(na.omit(dbinom(ysub, 1, psub)))
        alpha1[i,t] <- (alpha0[i,t-1] * col[i,t-1] + alpha1[i,t-1] * (1 - ext[i,t-1])) * cp  

        # Case when z = 0
        cp <- prod(na.omit(dbinom(ysub, 1, 0)))
        alpha0[i,t] <- (alpha0[i,t-1] * (1-col[i,t-1]) + alpha1[i,t-1] * ext[i,t-1]) * cp
      }
    }

    # Backwards pass
    beta1 <- matrix(NA, M, T)
    beta0 <- matrix(NA, M, T)

    # Initialize
    beta1[i, T] <- 1
    beta0[i, T] <- 1
  
    for (t in (T-1):1){
      ysub <- y[i, tind == t+1]
      psub <- p[i, tind == t+1]

      if(all(is.na(ysub))){
        beta1[i, t] <- ext[i,t] * beta0[i, t+1] + (1-ext[i,t]) * beta1[i, t+1]
        beta0[i, t] <- (1-col[i,t]) * beta0[i, t+1] + col[i,t] * beta1[i, t+1]
      } else {
        cp1 <- prod(na.omit(dbinom(ysub, 1, psub)))
        cp0 <- prod(na.omit(dbinom(ysub, 1, 0)))

        # Case when z = 1
        beta1[i, t] <- ext[i,t] * cp0 * beta0[i, t+1] + (1-ext[i,t]) * cp1 * beta1[i, t+1]
        # Case when z = 0
        beta0[i, t] <- (1-col[i,t]) * cp0 * beta0[i, t+1] + col[i,t] * cp1 * beta1[i, t+1]
      }
    }

    out <- rep(0, T)
    for (t in 1:T){
      out[t] <- (alpha1[i,t] * beta1[i,t]) / (alpha0[i,t]*beta0[i,t] + alpha1[i,t] * beta1[i,t])
    }

    smoothed[1, 1:T, i] <- 1 - out
    smoothed[2, 1:T, i] <- out
  }

  smoothed
}
