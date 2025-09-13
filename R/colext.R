colext <- function(psiformula = ~ 1, gammaformula = ~ 1,
                    epsilonformula = ~ 1, pformula = ~ 1,
                    data, starts = NULL, method = "BFGS", se = TRUE, ...){

  ## truncate to 1
  data@y <- truncateToBinary(data@y)
 
  formulas <- list(psi = psiformula, col = gammaformula, ext = epsilonformula,
                   det = pformula)
  check_no_support(formulas)
  dm <- getDesign(data, formulas)

  y <- dm$y
  M <- nrow(y)
  T <- data@numPrimary
  J <- ncol(y) / T

  # remove final year from transition prob design matrices
  dm$X_col <- dm$X_col[-seq(T,M*T,by=T),,drop=FALSE]
  dm$X_ext <- dm$X_ext[-seq(T,M*T,by=T),,drop=FALSE]

  # Set up submodels-----------------------------------------------------------
  estimateList <- unmarkedEstimateList(list(
    psi = unmarkedEstimate("Initial", short.name = "psi", invlink="logistic"),
    col = unmarkedEstimate("Colonization", short.name = "col", invlink="logistic"),
    ext = unmarkedEstimate("Extinction", short.name = "ext", invlink="logistic"),
    det = unmarkedEstimate("Detection", short.name = "p", invlink="logistic")
  ))

  # Set up parameter names and indices-----------------------------------------
  par_inds <- get_parameter_inds(estimateList, dm)

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

  tmb_inputs <- get_TMB_inputs(formulas, dm, par_inds, data,
                               M = M, T = T, J = J,
                               nd = no_detects, site_sampled = site_sampled)
  # Currently the TMB code requires y to be a vector
  tmb_inputs$data$y <- as.vector(t(tmb_inputs$data$y))

  fit <- fit_TMB("tmb_colext", starts, method, estimateList, par_inds,
                  tmb_inputs, data, ...)

  # Compute projected estimates
  psis <- plogis(dm$X_psi %*% coef(fit$estimate_list, "psi"))
  col_coef <- coef(fit$estimate_list, "col") 
  ext_coef <- coef(fit$estimate_list, "ext")
  det_coef <- coef(fit$estimate_list, "det")

  phis <- array(NA,c(2,2,T-1,M))
  phis[,1,,] <- plogis(dm$X_col %x% c(-1,1) %*% col_coef)
  phis[,2,,] <- plogis(dm$X_ext %x% c(-1,1) %*% -ext_coef)

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
                               col = plogis(dm$X_col %*% col_coef),
                               ext = plogis(dm$X_ext %*% ext_coef),
                               p = plogis(dm$X_det %*% det_coef),
                               M = M, T = T, J = J)
  smoothed.mean <- apply(smoothed, 1:2, mean)
  rownames(smoothed.mean) <- c("unoccupied","occupied")
  colnames(smoothed.mean) <- 1:T

  new("unmarkedFitColExt", fitType = "colext", call = match.call(),
      formlist = formulas, data = data, sitesRemoved = dm$removed.sites,
      estimates = fit$estimate_list, AIC = fit$AIC, opt = fit$opt, 
      negLogLike = fit$opt$value, nllFun = fit$nll, TMB = fit$TMB,
      projected = projected, projected.mean = projected.mean,
      smoothed = smoothed, smoothed.mean = smoothed.mean)
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
