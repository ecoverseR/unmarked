context("IDS fitting function")
skip_on_cran()

test_that("IDS can fit models with covariates", {
  set.seed(123)

  # Survey durations, loosely based on real data
  durs <- list(ds = rep(5, 1000), pc=runif(300, 3, 30))

  # Create distance data
  temp_ds <- unmarkedFrameDS(y = matrix(NA, 1000, 6),
                             siteCovs = data.frame(elev = rnorm(1000)),
                             dist.breaks = seq(0, 0.3, length.out=7),
                             survey = "point", unitsIn = "km")

  umf_ds <- simulate(temp_ds, formula = ~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2.5))[[1]]

  # Simulate availability
  avail_ds <- 1 - exp(-1 * durs$ds * plogis(-1.3))
  for (i in 1:nrow(umf_ds@y)){
    umf_ds@y[i,] <- rbinom(6, umf_ds@y[i,], avail_ds[i])
  }

  # Create point count data
  temp_pc <- unmarkedFrameDS(y = matrix(NA, 300, 1),
                             siteCovs = data.frame(elev = rnorm(300)),
                             dist.breaks = c(0, 0.5),
                             survey = "point", unitsIn = "km")

  umf_pc <- simulate(temp_pc, formula=~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2.5))[[1]]
  umf_pc <- unmarkedFramePCount(umf_pc@y, siteCovs = siteCovs(umf_pc))

  # Simulate availability
  avail_pc <- 1 - exp(-1 * durs$pc * exp(-1.3))
  for (i in 1:nrow(umf_pc@y)){
    umf_pc@y[i,] <- rbinom(1, umf_pc@y[i,], avail_pc[i])
  }

  set.seed(123)
  mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1,
                dataDS=umf_ds, dataPC=umf_pc,
                availformula = ~1, durationDS=durs$ds, durationPC=durs$pc,
                maxDistPC=0.5, unitsOut="kmsq")

  expect_equivalent(coef(mod_sim), 
                    c(3.12595,-0.453614,-2.487158,-1.789945), tol=1e-5)

  # Predict
  pr <- predict(mod_sim, type = 'lam')
  expect_equal(nrow(pr), 1000) # predicts only for the distance sampled sites
  pr <- predict(mod_sim, type = 'lam', newdata=umf_pc)

  pr <- predict(mod_sim, type = 'ds')
  expect_equal(nrow(pr), 1000)
  expect_equal(pr$Predicted[1], exp(-2.5), tol=0.05)
  pr <- predict(mod_sim, type = 'pc')
  expect_equal(nrow(pr), 300)
  expect_equal(pr$Predicted[1], exp(-2.5), tol=0.05)
  pr <- predict(mod_sim, type = 'phi')
  #expect_equal(pr$ds$Predicted[1], exp(-1.37), tol=0.05)
 
  # fitted
  ft <- fitted(mod_sim)
  expect_equal(lapply(ft, dim), list(ds=c(1000,6), pc = c(300,1)))
  expect_equal(round(ft[[1]],4)[1:2,1:2],
    structure(c(0.0648, 0.0654, 0.1369, 0.1381), dim = c(2L, 2L)))

  # residuals
  r <- residuals(mod_sim)
  expect_equal(lapply(r, dim), list(ds=c(1000,6), pc = c(300,1)))

  pdf(NULL)
  plot(mod_sim)
  hist(mod_sim)
  dev.off()

  # unsupported methods
  expect_error(parboot(mod_sim, nsim=2))
  expect_error(nonparboot(mod_sim))
  expect_error(ranef(mod_sim))

  # Separate detection models
  mod_sep <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaPC = ~1,
                dataDS=umf_ds[1:100,], dataPC=umf_pc[1:100,],
                availformula = ~1, durationDS=durs$ds[1:100], durationPC=durs$pc[1:100],
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")
  expect_equal(length(coef(mod_sim)), 4)
  expect_equal(length(coef(mod_sep)), 5)
})

test_that("IDS can fit models with occupancy data", {

  set.seed(123)
  # Create distance data
  temp_ds <- unmarkedFrameDS(y = matrix(NA, 100, 6),
                             siteCovs = data.frame(elev = rnorm(100)),
                             dist.breaks = seq(0, 0.3, length.out=7),
                             survey = "point", unitsIn = "km")

  umf_ds <- simulate(temp_ds, formula = ~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2.5))[[1]]

  # Create point count data
  temp_pc <- unmarkedFrameDS(y = matrix(NA, 100, 1),
                             siteCovs = data.frame(elev = rnorm(100)),
                             dist.breaks = c(0, 0.5),
                             survey = "point", unitsIn = "km")

  umf_pc <- simulate(temp_pc, formula=~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2.5))[[1]]
  umf_pc <- unmarkedFramePCount(umf_pc@y, siteCovs = siteCovs(umf_pc))

  # Create occupancy data
  temp_oc <- unmarkedFrameDS(y = matrix(NA, 100, 1),
                             siteCovs = data.frame(elev = rnorm(100)),
                             dist.breaks = c(0, 0.5),
                             survey = "point", unitsIn = "km")

  umf_oc <- simulate(temp_oc, formula=~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2))[[1]]
  umf_oc@y[,1] <- ifelse(umf_oc@y[,1] > 0, 1, 0)
  umf_oc <- unmarkedFrameOccu(umf_oc@y, siteCovs = siteCovs(umf_oc))

  mod_oc <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=umf_ds, dataPC=umf_pc, dataOC=umf_oc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")

  expect_equivalent(coef(mod_oc), c(2.87601, -0.578772, -2.499710, -1.9705832),
                    tol=1e-5)
  
  pr <- predict(mod_oc, type='oc')
  expect_equal(pr$Predicted[1], 0.1393755, tol=1e-5)
  
  res <- residuals(mod_oc)
  expect_equal(lapply(res, dim), list(ds=c(100,6), pc = c(100,1), oc=c(100,1)))

  # Don't estimate availability if OC data
  expect_error(
  mod_oc <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                availformula = ~1,
                dataDS=umf_ds, dataPC=umf_pc, dataOC=umf_oc,
                maxDistPC=0.5, maxDistOC=0.5,
                durationDS=durs$ds, durationPC=durs$pc, durationOC=durs$pc,
                unitsOut="kmsq")
  )

  # Just occupancy data
  mod_oc <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=umf_ds, dataOC=umf_oc,
                maxDistOC=0.5,
                unitsOut="kmsq")
  
  expect_equal(names(mod_oc), c("lam","ds","oc"))

})

test_that("IDS handles missing values", {

  set.seed(123)
  # Create distance data
  temp_ds <- unmarkedFrameDS(y = matrix(NA, 100, 6),
                             siteCovs = data.frame(elev = rnorm(100)),
                             dist.breaks = seq(0, 0.3, length.out=7),
                             survey = "point", unitsIn = "km")

  umf_ds <- simulate(temp_ds, formula = ~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2.5))[[1]]

  # Create point count data
  temp_pc <- unmarkedFrameDS(y = matrix(NA, 100, 1),
                             siteCovs = data.frame(elev = rnorm(100)),
                             dist.breaks = c(0, 0.5),
                             survey = "point", unitsIn = "km")

  umf_pc <- simulate(temp_pc, formula=~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2.5))[[1]]
  umf_pc <- unmarkedFramePCount(umf_pc@y, siteCovs = siteCovs(umf_pc))

  # Create occupancy data
  temp_oc <- unmarkedFrameDS(y = matrix(NA, 100, 1),
                             siteCovs = data.frame(elev = rnorm(100)),
                             dist.breaks = c(0, 0.5),
                             survey = "point", unitsIn = "km")

  umf_oc <- simulate(temp_oc, formula=~1~elev, unitsOut = "kmsq",
                     coefs = list(state = c(3, -0.5), det = -2))[[1]]
  umf_oc@y[,1] <- ifelse(umf_oc@y[,1] > 0, 1, 0)
  umf_oc <- unmarkedFrameOccu(umf_oc@y, siteCovs = siteCovs(umf_oc))

  umf_pc@y[1,1] <- NA
  umf_pc@y[2,] <- NA

  umf_oc@y[1,1] <- NA
  umf_oc@y[2,] <- NA

  expect_warning(
  mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=umf_ds, dataPC=umf_pc, dataOC=umf_oc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")
  )

  expect_equivalent(coef(mod_sim), c(2.881246,-0.57050,-2.4998305,-1.9608889),
                    tol=1e-5)

  ft <- fitted(mod_sim)
  expect_equal(dim(ft$ds), dim(umf_ds@y))
  expect_equal(dim(ft$pc), dim(umf_pc@y))
  expect_equal(dim(ft$oc), dim(umf_oc@y))

  umf_ds@siteCovs$elev[3] <- NA
  expect_warning(
  mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=umf_ds, dataPC=umf_pc, dataOC=umf_oc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq"))

  ft <- fitted(mod_sim)
  expect_equal(dim(ft$ds), dim(umf_ds@y))
  expect_true(all(is.na(ft$ds[3,]))) # missing site covariate

  umf_ds@y[1,1] <- NA
  umf_ds@y[2,] <- NA

  expect_error(
  expect_warning(
  mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=umf_ds, dataPC=umf_pc, dataOC=umf_oc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")
  ))

})

test_that("Hazard-rate works with IDS", {

  # Simulate IDS dataset
  set.seed(123)

  # First simulate distance sampling data
  # Distance sampling sites
  Mds <- 300
  elev <- rnorm(Mds)
  # Distance breaks
  db <- seq(0, 0.3, length.out=7)
  # blank unmarked frame
  umf_ds <- unmarkedFrameDS(y = matrix(NA, Mds, length(db)-1),
                            siteCovs = data.frame(elev = elev), 
                            survey="point", dist.breaks = db, unitsIn='km')

  # True coefficient values
  cf <- list(state = c(3, -0.5), # abundance intercept and effect of elevation 
            det = -2.5,         # detection sigma - exp(-2.5) = 0.08 
            scale = log(2))     # detection hazard rate scale = 2

  # Simulate response
  umf_ds <- simulate(umf_ds, formula=~1~elev, coefs=cf,
                     keyfun = "hazard", output = "density", unitsOut = "kmsq")[[1]]

  # Fit regular distance sampling model as test
  mod <- distsamp(~1~elev, umf_ds, keyfun="hazard", output="density", unitsOut = "kmsq")
  mod

  # good correspondance
  expect_equivalent(coef(mod), c(2.8913, -0.4381, -2.4697, 0.60489), tol=1e-4)

  # "Point count" sites
  # simulate these as 1-bin distance
  Mpc <- 500
  elev2 <- rnorm(Mpc)
  max_dist <- 0.5
  db2 <- c(0, max_dist) # single bin

  # blank unmarked frame
  umf_pc <- unmarkedFrameDS(y = matrix(NA, Mpc, 1), # single column in y
                          siteCovs = data.frame(elev = elev2),
                          survey = "point", dist.breaks=db2, unitsIn="km")

  # simulate response (same coefficients)
  umf_pc <- simulate(umf_pc, formula=~1~elev, coefs=cf,
                   keyfun = "hazard", output = "density", unitsOut = "kmsq")[[1]]


  # convert to unmarkedFramePcount
  umf_pc2 <- unmarkedFramePCount(y = umf_pc@y,
                                siteCovs = umf_pc@siteCovs)

  # Fit IDS model (in this case: common abundance and detection processes,
  # and common duration)
  mod2 <- expect_warning(IDS(lambdaformula = ~elev, detformulaDS = ~1,
              dataDS = umf_ds, dataPC = umf_pc2, keyfun = "hazard", 
              unitsOut = "kmsq", maxDistPC = 0.5))
 
  # similar to just distance sampling data
  expect_equivalent(coef(mod2),
                    c(2.8419,-0.5070, -2.4420, 0.6313), tol=1e-5)
  #cbind(true=unlist(cf), est=coef(mod2))

  # Check with separate formulas
  mod3 <- expect_warning(IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaPC = ~1,
              dataDS = umf_ds, dataPC = umf_pc2, keyfun = "hazard", 
              unitsOut = "kmsq", maxDistPC = 0.5))
  
  # Also make sure the coefs are in the right order: lam, then all sigma, then all haz scales
  expect_equal(coef(mod3),
                    c(`lam(Int)` = 2.86807, `lam(elev)` = -0.50746, 
                      `ds(Int)` = -2.47106, `pc(Int)` = -3.37866, 
                      `ds_scale(Int)` = 0.60447, `pc_scale(Int)` = -0.02817), tol=1e-5)
  
  # estimation doesn't really work for this one
  expect_warning(se <- SE(mod3))
  expect_true(is.nan(se["pc(Int)"]))
  expect_true(is.nan(se["pc_scale(Int)"]))

  # Make sure different formulas for sigma work
  mod4 <- expect_warning(IDS(lambdaformula = ~1, detformulaDS = ~elev, detformulaPC = ~1,
              dataDS = umf_ds, dataPC = umf_pc2, keyfun = "hazard", 
              unitsOut = "kmsq", maxDistPC = 0.5))
})
