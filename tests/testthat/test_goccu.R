context("goccu fitting function")
skip_on_cran()

set.seed(123)
M <- 100
T <- 5
J <- 4

psi <- 0.5
phi <- 0.3
p <- 0.4

z <- rbinom(M, 1, psi)
zmat <- matrix(z, nrow=M, ncol=T)

zz <- rbinom(M*T, 1, zmat*phi)
zz <- matrix(zz, nrow=M, ncol=T)

zzmat <- zz[,rep(1:T, each=J)]
y <- rbinom(M*T*J, 1, zzmat*p)
y <- matrix(y, M, J*T)
umf <- unmarkedMultFrame(y=y, numPrimary=T)

test_that("unmarkedFrameGOccu can be constructed", {
  set.seed(123)
  sc <- data.frame(x=rnorm(M))
  ysc <- matrix(rnorm(M*T), M, T)
  oc <- matrix(rnorm(M*T*J), M, T*J)

  umf2 <- unmarkedFrameGOccu(y, siteCovs=sc, obsCovs=list(x2=oc),
                           yearlySiteCovs=list(x3=ysc), numPrimary=T)
  expect_is(umf2, "unmarkedFrameGOccu")
  expect_equal(names(umf2@yearlySiteCovs), "x3")

  umf3 <- umf2[c(2,2,4),]
  expect_equal(numSites(umf3), 3)
  expect_equivalent(umf3[1,], umf2[2,])
  expect_equivalent(umf3[2,], umf2[2,])
  expect_equivalent(umf3[3,], umf2[4,])

  umf4 <- umf2[2:3,]
  expect_equal(numSites(umf4), 2)
  expect_equivalent(umf4[1,], umf2[2,])
  expect_equivalent(umf4[2,], umf2[3,])
})

test_that("goccu can fit models", {

  # Without covariates
  mod <- goccu(~1, ~1, ~1, umf)
  expect_equivalent(coef(mod), c(0.16129, -0.97041, -0.61784), tol=1e-5)
  
  # With covariates
  set.seed(123)
  sc <- data.frame(x=rnorm(M))
  ysc <- matrix(rnorm(M*T), M, T)
  oc <- matrix(rnorm(M*T*J), M, T*J)

  umf2 <- unmarkedMultFrame(y=y, siteCovs=sc, yearlySiteCovs=list(x2=ysc),
                            obsCovs=list(x3=oc), numPrimary=T)

  mod2 <- goccu(~x, ~x2, ~x3, umf2)
  expect_equivalent(coef(mod2), c(0.18895, -0.23629,-0.97246,-0.094335,-0.61808,
                                  -0.0040056), tol=1e-5)

  # predict
  pr <- predict(mod2, 'psi')
  expect_equal(dim(pr), c(M, 4))
  expect_equal(pr$Predicted[1], 0.5796617, tol=1e-5)

  # phi should not drop last level
  pr2 <- predict(mod2, 'phi')
  expect_equal(dim(pr2), c(M*T, 4))

  nd <- data.frame(x=1)
  pr3 <- predict(mod2, 'psi', newdata=nd)
  expect_true(nrow(pr3) == 1)
  expect_equal(pr3$Predicted[1], 0.488168, tol=1e-5)

  # Other methods
  ft <- fitted(mod2)
  expect_equal(dim(ft), dim(umf2@y))
  expect_true(all(ft >=0 & ft <= 1))
  expect_equal(round(ft,4)[1:2,1:2],
    structure(c(0.0583, 0.0529, 0.0586, 0.0531), dim = c(2L, 2L)))

  res <- residuals(mod2)
  expect_equal(dim(res), dim(umf2@y))

  gp <- getP(mod2)
  expect_equal(dim(gp), dim(umf2@y))
  expect_equal(as.vector(gp[1:2,1:2]), c(0.34923,0.35024,0.35088,0.35162), tol=1e-5) 

  set.seed(123)
  s <- simulate(mod2, nsim=2)
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(mod2@data@y))
  simumf <- umf2
  simumf@y <- s[[1]]
  simmod <- update(mod2, data=simumf)
  expect_equivalent(coef(simmod),
               c(0.174991, -0.27161, -1.32766, 0.054459,-0.41610,-0.073922), tol=1e-5)
  
  r <- ranef(mod2)
  expect_equal(dim(r@post), c(M, 2, 1))
  expect_equal(sum(bup(r)), 53.13565, tol=1e-4)
 
  pb <- parboot(mod2, nsim=2)
  expect_is(pb, "parboot")
  expect_equal(pb@t.star[1,1], 117.2043, tol=1e-4)

  npb <- nonparboot(mod2, B=2)
  expect_equal(length(npb@bootstrapSamples), 2)
  expect_equal(npb@bootstrapSamples[[1]]@AIC, 684.094, tol=1e-4)
  expect_true(npb@bootstrapSamples[[1]]@AIC != mod2@AIC)
  v <- vcov(npb, method='nonparboot')
  expect_equal(nrow(v), length(coef(npb)))
})

test_that("goccu handles missing values", {

  set.seed(123)
  y2 <- y
  y2[1,1] <- NA
  y2[2,1:J] <- NA

  sc <- data.frame(x=rnorm(M))
  sc$x[3] <- NA
  ysc <- matrix(rnorm(M*T), M, T)
  ysc[4,1] <- NA
  oc <- matrix(rnorm(M*T*J), M, T*J)
  oc[5,1] <- NA
  oc[6,1:J] <- NA

  umf_na <- unmarkedMultFrame(y=y2, siteCovs=sc, yearlySiteCovs=list(x2=ysc),
                            obsCovs=list(x3=oc), numPrimary=T)
 
  mod_na <- expect_warning(goccu(~x, ~x2, ~x3, umf_na))
  
  pr <- expect_warning(predict(mod_na, 'psi'))
  expect_equal(nrow(pr), M-1)

  # Need to re-write these to use the design matrix instead of predict
  gp <- getP(mod_na)
  expect_equal(dim(gp), c(100, 20))
  expect_true(is.na(gp[5,1]))
  expect_true(all(is.na(gp[6, 1:4])))
  s <- simulate(mod_na)
  expect_equal(dim(s[[1]]), dim(mod_na@data@y))
  ft <- fitted(mod_na)
  expect_equal(dim(ft), dim(mod_na@data@y))

  expect_true(all(is.na(ft[3,]))) # site covariate for site 3 missing
  expect_true(all(is.na(ft[4,1:4]))) # ysc covariate for site 4 per 1 missing
  expect_false(is.na(ft[4,5]))
  expect_true(is.na(ft[5,1])) # missing obs cov
  expect_true(all(is.na(ft[6,1:J]))) # missing obs cov
  
  r <- ranef(mod_na)
  expect_equal(dim(r@post), c(100, 2, 1))
  expect_true(is.na(bup(r)[3]))
 
  pb <- expect_warning(parboot(mod_na, nsim=2))
  expect_is(pb, "parboot")
})

test_that("larksparrow example from RMark", {
  # good test of missing primary periods
  umf_lasp <- readRDS("data/umf_lasp.Rds")
  fit1 <- goccu(psiformula = ~1, phiformula = ~1, pformula = ~1, umf_lasp)
  expect_equal(fit1@AIC, 331.2776, tol=1e-4)
  fit2 <- goccu(psiformula = ~1, phiformula = ~cw, pformula = ~td, umf_lasp)
  expect_equal(fit2@AIC, 323.4963, tol=1e-4)
  fit3 <- goccu(psiformula = ~ceap, phiformula = ~cw, pformula = ~td, umf_lasp)
  expect_equal(fit3@AIC, 323.2437, tol=1e-4)
})
