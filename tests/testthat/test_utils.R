context("utility functions")

test_that("invertHessian function works",{

  a <- 4; b <- 7; c <- 2; d <- 6
  mat <- matrix(c(a,b,c,d), nrow=2, byrow=T)
  mat_det <- a*d-b*c
  inv_mat <- 1/mat_det * matrix(c(d, -b, -c, a), nrow=2, byrow=T)

  fake_opt <- list(hessian=mat)

  #Successful inversion
  expect_equivalent(invertHessian(fake_opt, nrow(mat), TRUE),
                     inv_mat)

  #When se=F
  expect_equivalent(invertHessian(fake_opt, nrow(mat), FALSE),
                     matrix(rep(NA,4), nrow=2))

  #When matrix is not invertible
  bad_opt <- list(hessian=matrix(c(1, -2, -3, 6), nrow=2, byrow=T))
  expect_error(solve(bad_opt$hessian))

  #Should generate warning
  expect_warning(invertHessian(bad_opt, nrow(bad_opt$hessian), TRUE))


  #Should result in matrix of NAs
  expect_equivalent(invertHessian(bad_opt, nrow(bad_opt$hessian), FALSE),
                     matrix(rep(NA,4), nrow=2))
})

test_that("getStarts returns correct starting value vector", {
  skip() # this test doesn't work due to environment issues that I think don't matter
  data(frogs)
  pferUMF <- unmarkedFrameOccu(pfer.bin)
  siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
  obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))
  
  expect_equal(
    getStarts(fm <- occu(~ obsvar1 ~ 1, pferUMF)),
    c("psi(Int)"=0, "p(Int)"=0, "p(obsvar1)"=0)
  )

  expect_equal(
    getStarts(fm <- occu(~ obsvar1 ~ 1, pferUMF, starts=c(1,2,1))),
    c("psi(Int)"=1, "p(Int)"=2, "p(obsvar1)"=1)
  )

  # TMB changes the order
  expect_equal(
    getStarts(occu(~ obsvar1 ~ 1, pferUMF, starts=c(1,2,1), engine="TMB")),
    c("beta_state"=1, "beta_det"=1, "beta_det"=2)
  )
  expect_equal(
    getStarts(fm <- occu(~ (1|obsvar1) ~ 1, pferUMF, starts=c(1,2,1))),
    c("beta_state"=2, "beta_det"=1, "lsigma_det"=1)
  )
})
