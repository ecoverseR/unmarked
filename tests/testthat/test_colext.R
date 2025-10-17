context("colext fitting function")
skip_on_cran()

# Simulate data
set.seed(123)
nsites <- 6
nyr <- 4
nrep <- 2
y <- matrix(c(
        1,0, 1,1, 0,0, 0,0,
        1,1, 0,0, 0,0, 0,0,
        0,0, 0,0, 0,0, 0,0,
        0,0, 1,1, 0,0, 0,0,
        1,1, 1,0, 0,1, 0,0,
        0,0, 0,0, 0,0, 1,1), nrow=nsites, ncol=nyr*nrep, byrow=TRUE)

sc <- data.frame(sc1 = rnorm(nsites))
oc <- matrix(rnorm(nsites*nyr*nrep), nsites, nyr*nrep)
ysc <- matrix(rnorm(nsites*nyr), nsites, nyr)

test_that("unmarkedMultFrame construction works",{

  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, yearlySiteCovs=list(ysc=ysc),
                            obsCovs=list(oc=oc), numPrimary=4)
  expect_is(umf1, "unmarkedMultFrame")

  expect_error(unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc[1:5,]), numPrimary=4))

  expect_error(unmarkedMultFrame(y=y, siteCovs=sc, yearlySiteCovs=list(ysc=ysc),
                            obsCovs=list(oc=oc), numPrimary=3))

  plot(umf1)

  # Subsetting
  umf_sub <- umf1[2:3,]
  expect_equal(numSites(umf_sub), 2)
  expect_equivalent(umf_sub[1,], umf1[2,])
  expect_equivalent(umf_sub[2,], umf1[3,])
  umf_sub <- umf1[c(2,2,4),]
  expect_equivalent(umf_sub[1,], umf1[2,])
  expect_equivalent(umf_sub[2,], umf1[2,])
  expect_equivalent(umf_sub[3,], umf1[4,])
})


test_that("colext model fitting works", {

  set.seed(123)
  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)

  fm1 <- colext(~1, ~1, ~1, ~1, umf1)
    expect_equivalent(coef(fm1),
        c(0.14227, -1.49561,  0.210121,  1.200177),
        tol=1e-5)

  # With site covs
  fm2 <- colext(~sc1, ~1, ~1, ~1, umf1)
  expect_equivalent(coef(fm2), c(1.3473, -6.3828,-1.5830,0.1416,1.1643),
                    tol=1e-4)

  ft <- fitted(fm2)
  expect_equal(ft[1:3,1:3],
    structure(c(0.7566, 0.7191, 0.00014, 0.7566, 0.7191, 0.00014,
    0.3525, 0.3415, 0.1299), dim = c(3L, 3L)), tol = 1e-4)

  # With obs covs
  fm3 <- colext(~1, ~1, ~1, ~oc, umf1)
  expect_equivalent(coef(fm3),
        c(0.1434,-1.4980,0.2083,1.2006,-0.03827),
        tol=1e-4)

  # With yearly site covs
  fm4 <- colext(~1, ~ysc, ~ysc, ~1, umf1)
  expect_equivalent(coef(fm4),
                    c(0.2677,-2.0574,-1.0604,0.2156,0.6871,1.1025), tol=1e-4)

  # ranef
  r <- ranef(fm4)
  expect_equal(dim(r@post), c(nsites, nrep, nyr))
  expect_equal(dim(bup(r)), c(nsites, nyr))
  expect_equal(r@post[1,1,], c(0,0,0.9865,0.99098), tol=1e-4)

  # nonparboot
  expect_true(is.null(fm4@projected.mean.bsse))
  expect_true(is.null(fm4@smoothed.mean.bsse))
  npb <- nonparboot(fm4, B=2)
  expect_equal(length(npb@bootstrapSamples), 2)
  expect_equal(npb@bootstrapSamples[[1]]@AIC, 19.6418, tol=1e-4)
  v <- vcov(npb, method='nonparboot')
  expect_equal(nrow(v), length(coef(npb)))
  expect_is(npb@projected.mean.bsse, "matrix")
  expect_is(npb@smoothed.mean.bsse, "matrix")
   
  # parboot
  pb <- parboot(fm4, nsim=2)
  expect_equal(pb@t.star[1,1], 10.2837, tol=1e-4)

  # getP
  gp <- getP(fm4)
  expect_equal(dim(gp), c(6,8))
  expect_equal(gp[1,1], 0.7507, tol=1e-4)
})

test_that("colext handles missing values",{

  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)

  umf2 <- umf1
  umf2@y[1,3] <- NA

  fm1 <- colext(~1, ~1, ~1, ~1, umf2)
  expect_equal(fm1@AIC, 50.48498, tol=1e-4)

  umf3 <- umf1
  umf3@y[1,] <- NA
  expect_warning(fm2 <- colext(~1, ~1, ~1, ~1, umf3))
  expect_equal(fm2@AIC, 41.87926, tol=1e-4)
  expect_equal(fm2@sitesRemoved, 1)

  umf4 <- umf1
  umf4@y[1,3:4] <- NA
  fm3 <- colext(~1, ~1, ~1, ~1, umf4)
  expect_equal(fm3@AIC, 47.74454, tol=1e-4)

  umf5 <- umf1
  umf5@siteCovs$sc1[2] <- NA
  umf5@obsCovs$oc[1] <- NA
  expect_warning(fm4 <- colext(~sc1, ~1, ~1, ~oc, umf5))
  expect_warning(pr <- predict(fm4, 'det'))
  expect_equal(nrow(pr), (nsites-1)*nyr*nrep)
  expect_true(all(is.na(pr[1,])))
  ft <- fitted(fm4)
  expect_equal(dim(ft), dim(umf5@y))
  expect_true(is.na(ft[1,1]))
  expect_true(all(is.na(ft[2,])))

  gp <- getP(fm4)
  expect_equal(dim(gp), dim(umf5@y))
  expect_true(all(!is.na(gp[2,])))
  expect_equal(as.vector(gp[1:2,1:2]), c(NA, 0.7318,0.8017,0.7906), tol=1e-4)

  r <- ranef(fm4)
  expect_true(all(is.na(r@post[fm4@sitesRemoved,,1])))
  expect_equal(nrow(r@post), numSites(fm4@data))

  umf5 <- umf1
  umf5@yearlySiteCovs$ysc[1] <- NA
  # This should work, right?
  expect_error(expect_warning(fm4 <- colext(~1, ~1, ~ysc, ~1, umf5)))

})

test_that("colext errors when random effects are in formula",{
  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)
  expect_error(colext(~(1|dummy), ~1, ~ysc, ~1, umf1))
})

test_that("colext methods work",{

  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)
  fm1 <- colext(~sc1, ~1, ~ysc, ~oc, umf1)

  pdf(NULL)
  plot(fm1)
  dev.off()
  res <- residuals(fm1)
  expect_equal(dim(res), c(6,8))
  r <- ranef(fm1)
  expect_equal(dim(r@post), c(nsites, nrep, nyr))
  pr1 <- predict(fm1, 'psi')
  expect_equal(nrow(pr1), 6)
  pr2 <- predict(fm1, 'col')
  expect_equal(nrow(pr2), nsites*nyr)
  pr3 <- predict(fm1, 'det')
  expect_equal(nrow(pr3), nsites*nyr*nrep)

  nd <- data.frame(sc1=c(0,1))
  pr4 <- predict(fm1, 'psi', newdata=nd)
  expect_equal(nrow(pr4), 2)
})

test_that("Missing primary periods are handled correctly", {
  #https://github.com/ecoverseR/unmarked/pull/40
  y <- readRDS("data/colext_missing_periods.Rds")
  umf <- unmarkedMultFrame(y = y, numPrimary = 9)
  fit <- colext(data = umf)
  expect_equal(round(unname(coef(fit)), 4),
               c(-0.0254, -2.1375, -2.0540, 0.8759))
})
