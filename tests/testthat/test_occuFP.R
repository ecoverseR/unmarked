context("occuFP fitting function")
skip_on_cran()

test_that("occuFP model can be fit",{
  set.seed(123)
  n = 100
  o = 10
  o1 = 5
  y = matrix(0,n,o)
  p = .7
  r = .5
  fp = 0.05
  y[1:(n*.5),(o-o1+1):o] <- rbinom((n*o1*.5),1,p)
  y[1:(n*.5),1:(o-o1)] <- rbinom((o-o1)*n*.5,1,r)
  y[(n*.5+1):n,(o-o1+1):o] <- rbinom((n*o1*.5),1,fp)
  type <- c((o-o1),o1,0)
  site <- c(rep(1,n*.5*.8),rep(0,n*.5*.2),rep(1,n*.5*.2),rep(0,n*.8*.5))
  occ <- matrix(c(rep(0,n*(o-o1)),rep(1,n*o1)),n,o)
  site <- data.frame(habitat = site)
  occ <- list(METH = occ)
  umf1 <- unmarkedFrameOccuFP(y,site,occ, type = type)

  m1 <- occuFP(detformula = ~ METH, FPformula = ~1,
               stateformula = ~ habitat, data = umf1)
  expect_equal(names(m1), c("state","det","fp"))
  expect_equal(round(unname(coef(m1)), 4),
               c(-1.3750, 2.7602, -0.1954, 1.1562, -3.0884))

  expect_message(fl <- fitList(m1,m1))
  expect_is(fl,"unmarkedFitList")
  expect_equal(length(fl@fits), 2)

  pr <- predict(m1, "fp")
  expect_equal(dim(pr), c(1000, 4))

  # Check update
  new_fit <- update(m1, stateformula=~1, data=umf1[1:10,])
  expect_equal(length(coef(new_fit, "state")), 1)
  expect_equal(nrow(new_fit@data@y), 10)

  # Check parboot (SSE doesn't work)
  pb <- parboot(new_fit, statistic=function(x) x@AIC, nsim=2)
  expect_equal(pb@t.star[1,1], 109.8942, tol=1e-4)

  # Nonparboot
  npb <- nonparboot(new_fit, B=2)
  expect_equal(length(npb@bootstrapSamples), 2)
  expect_equal(npb@bootstrapSamples[[1]]@AIC, 136.29, tol=1e-4)
  expect_equal(numSites(npb@bootstrapSamples[[1]]@data), numSites(npb@data))
  v <- vcov(npb, method='nonparboot')
  expect_equal(nrow(v), length(coef(npb)))

  # getP
  gp <- getP(m1)
  expect_equal(dim(gp), dim(m1@data@y))
  expect_equal(as.vector(gp[5:6, 5:6]), c(0.4513,0.4513,0.7232,0.7232), tol=1e-4)

  # ranef
  expect_error(r <- ranef(m1))

  # Check error when random effect in formula
  expect_error(occuFP(~(1|dummy), ~1, ~1, data=umf1))
})

test_that("occuFP model can be fit with missing values",{
  set.seed(123)
  n = 100
  o = 10
  o1 = 5
  y = matrix(0,n,o)
  p = .7
  r = .5
  fp = 0.05
  y[1:(n*.5),(o-o1+1):o] <- rbinom((n*o1*.5),1,p)
  y[1:(n*.5),1:(o-o1)] <- rbinom((o-o1)*n*.5,1,r)
  y[(n*.5+1):n,(o-o1+1):o] <- rbinom((n*o1*.5),1,fp)
  type <- c((o-o1),o1,0)
  site <- c(rep(1,n*.5*.8),rep(0,n*.5*.2),rep(1,n*.5*.2),rep(0,n*.8*.5))
  occ <- matrix(c(rep(0,n*(o-o1)),rep(1,n*o1)),n,o)
  site <- data.frame(habitat = site)
  occ <- list(METH = occ)

  y[1,1] <- NA
  y[2,] <- NA
  site$habitat[3] <- NA
  occ$METH[4,1] <- NA
  occ$METH[5,] <- NA

  umf1 <- unmarkedFrameOccuFP(y,site,occ, type = type)

  m1 <- expect_warning(occuFP(detformula = ~ METH, FPformula = ~1,
                        stateformula = ~ habitat, data = umf1))
  expect_equal(m1@sitesRemoved, c(2,3,5))
  expect_equal(round(unname(coef(m1)), 4),
               c(-1.3750, 2.6815, -0.2011, 1.1588, -3.0876))

  expect_message(expect_warning(fl <- fitList(m1,m1)))
  expect_is(fl,"unmarkedFitList")
  expect_equal(length(fl@fits), 2)

  pr <- expect_warning(predict(m1, "fp"))
  expect_equal(dim(pr), c(1000-10*length(m1@sitesRemoved), 4))

  # Check update
  new_fit <- expect_warning(update(m1, stateformula=~1, data=umf1[1:10,]))
  expect_equal(length(coef(new_fit, "state")), 1)
  expect_equal(nrow(new_fit@data@y), 10)

  # Check parboot (SSE doesn't work)
  pb <- expect_warning(parboot(new_fit, statistic=function(x) x@AIC, nsim=2))
  expect_equal(pb@t.star[1,1], 96.28742, tol=1e-4)

  # Nonparboot
  npb <- expect_warning(nonparboot(new_fit, B=2))
  expect_equal(length(npb@bootstrapSamples), 2)
  expect_equal(npb@bootstrapSamples[[1]]@AIC, 81.46887, tol=1e-4)
  expect_equal(numSites(npb@bootstrapSamples[[1]]@data), numSites(npb@data))
  v <- vcov(npb, method='nonparboot')
  expect_equal(nrow(v), length(coef(npb)))

  # getP
  gp <- getP(m1)
  expect_equal(dim(gp), dim(m1@data@y))
  expect_true(is.na(gp[4,1]))
  expect_true(all(is.na(gp[5,])))
  expect_equal(gp[1,1], 0.4499, tol=1e-4)
})

test_that("occuFP certain detection/b model can be fit",{
  # Example from AHM2 7.4.1
  set.seed(129)       # RNG seed
  nsites <- 200       # number of sites
  nsurveys1 <- 3      # number of occasions with Type 1 data
  nsurveys2 <- 4      # number of occasions with Type 2 data
  psi <- 0.6          # expected proportion of are occupied
  p <- c(0.7,0.5)     # detection prob of method 1 and method 2
  fp <- 0.05          # false-positive error probability (p_10)
  b <- 0.2            # probability y is recorded as certain

  # Simulate the occupancy states and data
  z <- rbinom(nsites, 1, psi)
  y <- matrix(NA, nrow = nsites, ncol = nsurveys1 + nsurveys2)
  for(i in 1:nsites){
    p1 <- p[1]*z[i]                      # certainly detection (method 1)
    p2 <- p[2]*z[i] + fp*(1-z[i])        # uncertainly detection (method 2)
    y[i,1:3] <- rbinom(nsurveys1, 1, p1) # simulate method 1 data
    y[i,4:7] <- rbinom(nsurveys2, 1, p2) # simulate method 2 data
    # Now introduce certain observations:
    pr.certain <- z[i] * y[i,4:7] * b
    y[i, 4:7] <- y[i, 4:7] + rbinom(4, 1, pr.certain)
  }

  y[1,1] <- NA
  y[2,] <- NA

  # Make a covariate to distinguish between the two methods
  Method <- matrix(c(rep("1", 3), rep("2", 4)), nrow = nsites,
      ncol = nsurveys1 + nsurveys2, byrow = TRUE)

  type <- c(nsurveys1, 0, nsurveys2)

  umf2 <- expect_warning( unmarkedFrameOccuFP(y = y, obsCovs = list(Method = Method),
      type = type))
  m2 <- expect_warning(occuFP(detformula = ~ -1 + Method, FPformula = ~ 1, Bformula = ~ 1,
    stateformula = ~ 1, data = umf2))
  expect_equal(m2@sitesRemoved, 2)
  expect_equal(m2@AIC, 1720.162, tol=1e-4)
})
