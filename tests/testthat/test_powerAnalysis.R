context("powerAnalysis method")
skip_on_cran()

temp <- unmarkedFrameOccu(y=matrix(NA, 300, 8), 
          siteCovs=data.frame(elev=rnorm(300),
                              group=factor(sample(letters[1:20], 300, replace=TRUE))))

test_that("powerAnalysis method works",{

  # When no effect sizes provided
  nul <- capture.output(expect_error(powerAnalysis(temp, model=occu, formula=~1~elev)))

  ef <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))

  nul <- capture_output({

  set.seed(123)
  pa <- powerAnalysis(temp, model=occu, formula=~1~elev, effects=ef, nsim=10)
  expect_is(pa, "unmarkedPower")
  s <- summary(pa)$Power
  expect_true(s>0.5)

  # output printout
  out <- capture.output(pa)
  expect_equal(out[7], "Power Statistics:")

  # set null
  set.seed(123)
  nul <- list(state=c(intercept=5, elev=0), det=c(intercept=0))
  pa2 <- powerAnalysis(temp, model=occu, formula=~1~elev, effects=ef, nulls=nul, nsim=10)
  s <- summary(pa2, showIntercepts=TRUE)
  expect_equivalent(s$Null, c(5,0,0))

  # list
  pl <- unmarkedPowerList(pa, pa2)
  expect_is(pl, "unmarkedPowerList")
  s <- summary(pl)
  expect_is(s, "data.frame")

  pdf(NULL)
  pl_plot <- plot(pl)
  expect_is(pl_plot,"list")
  dev.off()

  # With random effect
  set.seed(123) 
  ef <- list(state=c(0, -0.4, 1), det=0)
  pa3 <- powerAnalysis(temp, model=occu, formula=~1~elev+(1|group), effects=ef, nsim=10)
  s <- summary(pa3, showIntercepts=TRUE)
  expect_equal(nrow(s), 3)

  # Only one random effect allowed
  expect_error(powerAnalysis(temp, model=occu, formula=~1~elev+(elev||group), effects=ef, nsim=10))
  })
})

test_that("custom datasets can be passed to powerAnalysis",{
  set.seed(123)
  ef <- list(state=c(0, -0.4), det=0)
  s <- simulate(temp, model=occu, formula=~1~elev, nsim=10, coefs=ef, quiet=TRUE)
  
  nul <- capture_output({
    pa <- powerAnalysis(s, model=occu, formula=~1~elev, effects=ef)
    expect_is(pa, "unmarkedPower")
    expect_equal(length(pa@estimates), length(s))
  })
})

test_that("powerAnalysis can be run in parallel",{
  skip_on_cran()
  skip_on_ci()
  #skip_if(!requireNamespace("pbapply", quietly=TRUE), 
  #        "pbapply package unavailable")

  set.seed(123)
  ef <- list(state=c(0, -0.4), det=0)

  nul <- capture_output({
    pa <- powerAnalysis(temp, model=occu, formula=~1~elev, effects=ef, parallel=TRUE, nsim=2)
    expect_is(pa, "unmarkedPower")
    expect_equal(length(pa@estimates), 2)
  })
})

test_that("Error returned for unsupported unmarkedFit-based power analysis",{
  temp2 <- temp
  temp2@y[] <- 0
  fit <- suppressWarnings(occu(~1~1, temp2))
  expect_error(powerAnalysis(fit), "no longer supported")

  expect_error(shinyPower(fit), "no longer supported")
})
