
#' Fit the N-mixture point count model

pcount <- function(formula, data, K, mixture = c("P", "NB", "ZIP"), 
                   starts = NULL, method = "BFGS", se = TRUE,
                   engine = c("C", "R", "TMB"), threads = 1, ...)
{

    # Check argument validity--------------------------------------------------
    mixture <- match.arg(mixture, c("P", "NB", "ZIP"))
    mixture_code <- switch(mixture, P = {1}, NB = {2}, ZIP = {3})

    if(!is(data, "unmarkedFramePCount"))
        stop("Data is not an unmarkedFramePCount object.")

    engine <- match.arg(engine, c("C", "R", "TMB"))
    formulas <- split_formula(formula)
    names(formulas) <- c("det", "state")
    if(any(sapply(formulas, has_random))) engine <- "TMB"
    if(identical(mixture, "ZIP") & engine == "R")
        stop("ZIP mixture not available for R engine")

    # Generate design matrices-------------------------------------------------
    dm <- getDesign(data, formulas)
    y <- dm$y

    # Set up submodels-----------------------------------------------------------
    state <- unmarkedEstimate("Abundance", short.name = "lam", invlink="exp")
    det <- unmarkedEstimate("Detection", short.name = "p", invlink="logistic")
    estimateList <- unmarkedEstimateList(list(state=state,det=det))

    if(identical(mixture,"NB")) {
      estimateList@estimates$alpha <- 
        unmarkedEstimate(name="Dispersion", short.name = "alpha", invlink = "exp")
    } else if(identical(mixture,"ZIP")) {
      estimateList@estimates$psi <- 
        unmarkedEstimate(name="Zero-inflation", short.name = "psi", invlink = "logistic")
    }

    # Set up parameter names and indices-----------------------------------------
    par_inds <- get_parameter_inds(estimateList, dm)

    # Handle K (number of possible abundance values to marginalize over)-------
    if(missing(K)) {
        K <- max(y, na.rm = TRUE) + 100
        warning("K was not specified and was set to ", K, ".")
    }
    if(K <= max(y, na.rm = TRUE))
        stop("specified K is too small. Try a value larger than any observation")

    #Minimum observed abundance at each site: Used by C++ and TMB
    Kmin <- apply(y, 1, function(x) max(x, na.rm=TRUE))

    # Specify negative log likelihood functions--------------------------------
    if(identical(engine, "R")) {
        k <- 0:K
        M <- nrow(y)
        J <- ncol(y)
        k.ik <- rep(k, M)
        k.ijk <- rep(k, M*J)
        y.ij <- as.numeric(t(y))
        y.ijk <- rep(y.ij, each = K + 1)
        navec <- is.na(y.ijk)
        ijk <- expand.grid(k = 0:K, j = 1:J, i = 1:M)
        ijk.to.ikj <- with(ijk, order(i, k, j))
        nll <- function(parms) {
            theta.i <- exp(dm$X_state %*% parms[par_inds$state] + dm$offset_state)
            p.ij <- plogis(dm$X_det %*% parms[par_inds$det] + dm$offset_det)
            theta.ik <- rep(theta.i, each = K + 1)
            p.ijk <- rep(p.ij, each = K + 1)

            bin.ijk <- dbinom(y.ijk,k.ijk,p.ijk)
            bin.ijk[which(is.na(bin.ijk))] <- 1
            bin.ik.mat <- matrix(bin.ijk[ijk.to.ikj], M * (K + 1), J,
                                 byrow = TRUE)
            g.ik <- rowProds(bin.ik.mat)

            if(identical(mixture,"P")) {
                f.ik <- dpois(k.ik,theta.ik)
            }
            else if (identical(mixture,"NB")){
                f.ik <- dnbinom(k.ik, mu = theta.ik, size = exp(parms[par_inds$alpha]))
            }
            dens.i.mat <- matrix(f.ik * g.ik, M, K + 1, byrow = TRUE)
            dens.i <- rowSums(dens.i.mat)  # sum over the K

            -sum(log(dens.i))
      }
    } else if(identical(engine, "C")) {
        n_param <- sapply(par_inds, length)
        if(length(n_param) == 2) n_param <- c(n_param, 0)
        nll <- function(parms) {
          nll_pcount(parms, n_param, y, dm$X_state, dm$X_det, dm$offset_state, dm$offset_det, 
                     K, Kmin, mixture_code, threads)
        }
    }

    # Fit model in C or R------------------------------------------------------
    if(engine %in% c("C","R")){

      fit <- fit_optim(nll, starts, method, se, estimateList, par_inds, ...)

    # Fit model in TMB---------------------------------------------------------
    } else if(engine == "TMB"){

      # Set up TMB input data
      tmb_inputs <- get_TMB_inputs(formulas, dm, par_inds, data, 
                                   K=K, Kmin=Kmin, mixture=mixture_code)

      # Fit model with TMB
      fit <- fit_TMB("tmb_pcount", starts, method, estimateList, par_inds,
                      tmb_inputs, data, ...)
    }

    # Create unmarkedFit object------------------------------------------------
    new("unmarkedFitPCount", fitType="pcount", call=match.call(),
        formula = formula, formlist = formulas, data = data,
        sitesRemoved = dm$removed.sites,
        estimates = fit$estimate_list, AIC = fit$AIC, opt = fit$opt,
        negLogLike = fit$opt$value,
        nllFun = fit$nll, K = K, mixture = mixture, TMB=fit$TMB)
}
