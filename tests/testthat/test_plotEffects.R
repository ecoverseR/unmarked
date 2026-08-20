context("plotEffects")

skip_on_cran()

set.seed(123)
dat_occ <- data.frame(x1=rnorm(500))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
z <- rep(NA, 500)

b <- c(0.4, -0.5, 0.3, 0.5)

re_fac <- factor(sample(letters[1:5], 500, replace=T))
dat_occ$group <- re_fac
re <- rnorm(5, 0, 1.2)
re_idx <- as.numeric(re_fac)

idx <- 1
for (i in 1:500){
  z[i] <- rbinom(1,1, plogis(b[1] + b[2]*dat_occ$x1[i] + re[re_idx[i]]))
  for (j in 1:5){
    y[i,j] <- z[i]*rbinom(1,1,
                    plogis(b[3] + b[4]*dat_p$x2[idx]))
    idx <- idx + 1
  }
}

umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)

fm <- occu(~x2 ~x1 + group, umf)

test_that("plotEffects works", {
  
  pdf(NULL)
  plotEffects(fm, "state", "x1")
  plotEffects(fm, "state", "group")
  plotEffects(fm, "det", "x2")
  dev.off()

  expect_error(plotEffects(fm, "state", "x2"))

  dat <- plotEffectsData(fm, "state", "group")

  expect_true(inherits(dat, "data.frame"))
  expect_equal(nrow(dat), 5)

})

test_that("plotEffects works when there are non numeric/factor columns in covs", {

  # #55 https://github.com/ecoverseR/unmarked/issues/55
  data(frogs)
  
  # 1. Create a basic UMF and fit model 
  pferUMF <- unmarkedFrameOccu(pfer.bin)
  siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
  fm1 <- occu(~ 1 ~ sitevar1, pferUMF)
  
  # This works as expected
  expect_is(plotEffectsData(fm1, type = "state", covariate = "sitevar1"),
            "data.frame")
  
  # 2. Add a Date column
  # Note: 'date' is NOT used in the model formula
  siteCovs(pferUMF)$date <- as.Date("2026-01-01")
  fm2 <- occu(~ 1 ~ sitevar1, pferUMF)
  
  expect_is(plotEffectsData(fm2, type = "state", covariate = "sitevar1"),
            "data.frame")
  
  # 3. Add a Logical column
  siteCovs(pferUMF)$date <- NULL
  siteCovs(pferUMF)$logic_col <- TRUE
  fm3 <- occu(~ 1 ~ sitevar1, pferUMF)
  
  expect_is(plotEffectsData(fm3, type = "state", covariate = "sitevar1"),
            "data.frame")

})
