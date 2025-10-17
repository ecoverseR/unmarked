context("simulate method")
skip_on_cran()

test_that("simulate can generate new datasets from scratch",{

  set.seed(123)

  y <- matrix(NA, 300, 5)
  sc <- data.frame(elev=rnorm(300))

  umf <- unmarkedFrameOccu(y=y, siteCovs=sc)

  # When no coefficients list provided
  nul <- capture_output(expect_error(simulate(umf, model=occu, formula=~1~elev)))

  cf <- list(state=c(0, -0.4), det=0)

  # When model not provided
  expect_error(simulate(umf, formula=~1~elev, coefs=cf))

  s <- expect_message(simulate(umf, model=occu, formula=~1~elev, coefs=cf)[[1]])

  expect_equivalent(dim(s@y), c(300,5))
  expect_equal(names(s@siteCovs), "elev")

  fm <- occu(~1~elev, s)
  expect_equivalent(coef(fm), c(-0.1361,-0.5984,-0.002639), tol=1e-4)

  # With random effect
  set.seed(123)
  umf@siteCovs$group <- factor(sample(letters[1:20], 300, replace=TRUE))

  cf2 <- list(state=c(0, -0.4, 1), det=0)

  s <- expect_message(simulate(umf, model=occu, formula=~1~elev+(1|group), coefs=cf2)[[1]])

  fm <- occu(~1~elev+(1|group), s)
  expect_equal(sigma(fm)$sigma, 1.04565, tol=1e-4)

  # pcount
  set.seed(123)
  temp <- unmarkedFramePCount(y=y, siteCovs=sc)
  cf$alpha <- c(alpha=0.5)
  s <- expect_message(simulate(temp, formula=~1~elev, K=10, mixture="NB", coefs=cf)[[1]])

  fm2 <- pcount(~1~elev, s, mixture="NB", K=10)
  expect_equivalent(coef(fm2), c(-0.0843,-0.3777,-0.0505,0.666), tol=1e-3)

  # distsamp
  set.seed(123)
  cf$alpha <- NULL
  cf$det[1] <- log(30)
  cf$state <- c(intercept=2, elev=0.5)
  temp <- unmarkedFrameDS(y=y, siteCovs=sc, dist.breaks=c(0,10,20,30,40,50), survey='point', unitsIn='m')
  s <- expect_message(simulate(temp, formula=~1~elev, coefs=cf)[[1]])
  fm <- distsamp(~1~elev, s)
  expect_equivalent(coef(fm), c(1.9734, 0.5283, 3.403), tol=1e-4)

  # Mpois
  set.seed(123)
  cf$dist[1] <- 0
  cf$state <- c(intercept=1, elev=0.5)
  cf$det <- 0
  temp <- unmarkedFrameMPois(y=y, siteCovs=sc, type='removal')
  s <- expect_message(simulate(temp, formula=~1~elev, coefs=cf)[[1]])
  fm <- multinomPois(~1~elev, s)
  expect_equivalent(coef(fm), c(0.975,0.513,0.112), tol=1e-3)

  #colext
  set.seed(123)
  cf_colext <- list(psi=c(intercept=0, elev=0.5), col=c(intercept=0),
                    ext=c(intercept=0), det=c(intercept=0))
  y_ce <- matrix(NA, 300, 15)
  temp <- unmarkedMultFrame(y_ce, siteCovs=sc, numPrimary=3)

  s <- expect_message(simulate(temp, psiformula=~elev, coefs=cf_colext)[[1]])
  fm <- colext(~elev, ~1, ~1, ~1, s)
  expect_equivalent(coef(fm), c(-0.132,0.381,0.0701,0.158,0.015),
                     tol=1e-3)

  #occuTTD
  set.seed(123)
  cf_ttd <- cf_colext
  cf_ttd$det <- c(intercept=log(0.5))
  temp <- unmarkedFrameOccuTTD(y_ce, siteCovs=sc, numPrimary=3, surveyLength=3)

  s <- expect_message(simulate(temp, model=occuTTD, psiformula=~elev, coefs=cf_ttd)[[1]])
  fm <- occuTTD(~elev, ~1, ~1, ~1, s)
  expect_equivalent(coef(fm), c(0.115,0.642,-0.065,-0.095,-0.693),
                     tol=1e-3)

  #gdistsamp
  set.seed(123)
  cf_gds <- list(det=c(intercept=log(30)), lambda=c(intercept=1, elev=0.5),
                 phi=c(intercept=0))

  temp <- unmarkedFrameGDS(y_ce, siteCovs=sc, numPrimary=3, dist.breaks=c(0,10,20,30,40,50), survey='line',
                    tlength=rep(100,300), unitsIn='m')

  s <- expect_message(simulate(temp, lambdaformula=~elev, phiformula=~1, pformula=~1, coefs=cf_gds)[[1]])
  fm <- gdistsamp(~elev, ~1, ~1, data=s, K=15)
  expect_equivalent(coef(fm), c(0.999, 0.451, -0.010, 3.37),
                     tol=1e-3)

  #gmultmix
  set.seed(123)
  cf_gmm <- list(det=c(0), lambda=c(intercept=1, elev=0.5),
                 phi=c(intercept=0))
  temp <- unmarkedFrameGMM(y_ce, siteCovs=sc, numPrimary=3, type='removal')
  s <- expect_message(simulate(temp, lambdaformula=~elev, phiformula=~1, pformula=~1, coefs=cf_gmm, K=15)[[1]])
  fm <- gmultmix(~elev, ~1, ~1, data=s, K=15)
  expect_equivalent(coef(fm), c(1.0025,0.4762,0.022,-0.04318),
                     tol=1e-3)

  #gpcount
  set.seed(123)
  cf_gmm$lambda <- c(0.5, 0.5)
  temp <- unmarkedFrameGPC(y_ce, siteCovs=sc, numPrimary=3)
  s <- expect_message(simulate(temp, lambdaformula=~elev, phiformula=~1, pformula=~1, coefs=cf_gmm, K=10)[[1]])
  fm <- gpcount(~elev, ~1, ~1, data=s, K=10)
  expect_equivalent(coef(fm), c(0.428,0.525,0.0885,-0.040),
                     tol=1e-3)

  #pcountOpen
  set.seed(123)
  cf_pco <- list(lambda=c(intercept=2, elev=0.5), det=c(intercept=0),
                 gamma=c(intercept=0), omega=c(intercept=0))
  y_pco <- matrix(NA, 100, 15)
  temp <- unmarkedFramePCO(y_pco, siteCovs=data.frame(elev=rnorm(100)), numPrimary=3)

  s <- expect_message(simulate(temp, lambdaformula=~elev, gammaformula=~1, 
                               omegaformula=~1, pformula=~1, coefs=cf_pco, K=30)[[1]])

  fm <- pcountOpen(~elev, ~1, ~1, ~1, data=s, K=30)
  expect_equivalent(coef(fm), c(1.9802,0.4691,-0.0366,-0.0054,0.0323), tol=1e-4)

  #multmixOpen
  set.seed(123)
  temp <- unmarkedFrameMMO(y_pco, siteCovs=data.frame(elev=rnorm(100)), numPrimary=3, type='removal')
  s <- expect_message(simulate(temp, lambdaformula=~elev, gammaformula=~1, 
                               omegaformula=~1, pformula=~1, coefs=cf_pco, K=15)[[1]])
  expect_is(s, "unmarkedFrameMMO")

  #distsampOpen
  set.seed(123)
  cf_dso <- cf_pco
  cf_dso$det <- c(intercept=log(30))
  temp <- unmarkedFrameDSO(y_pco, siteCovs=data.frame(elev=rnorm(100)), numPrimary=3,
                          dist.breaks=c(0,10,20,30,40,50),unitsIn='m', survey='point')
  s <- expect_message(simulate(temp, lambdaformula=~elev, gammaformula=~1, 
                               omegaformula=~1, pformula=~1, coefs=cf_dso, K=15)[[1]])
  expect_is(s, "unmarkedFrameDSO")

  # occuMulti
  set.seed(123)
  occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3','~1','~1','~1','~1')
  detFormulas <- c('~1','~1','~1')
  beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
  p_true <- c(0.6,0.7,0.5)
  cf <- list(state=beta, det=log(p_true/(1-p_true)))
  sc <- data.frame(occ_cov1=rnorm(300), occ_cov2=rnorm(300), occ_cov3=rnorm(300))
  temp <- unmarkedFrameOccuMulti(list(sp1=y, sp2=y, sp3=y), siteCovs=sc)
  s <- expect_message(simulate(temp, stateformulas=occFormulas, detformulas=detFormulas,
                               coefs=cf)[[1]])
  fm <- occuMulti(stateformulas=occFormulas, detformulas=detFormulas, s)
  expect_equivalent(coef(fm, 'det'), c(0.2982,0.8416,-0.01816), tol=1e-4)

  # occuMS
  set.seed(123)
  bstate <- c(-0.5, 1, -0.6, -0.7)
  bdet <- c(-0.4, 0, -1.09, -0.84)
  detformulas <- c('~V1','~1','~1')
  stateformulas <- c('~V1','~V2')
  cf <- list(state=bstate, det=bdet)
  sc <- data.frame(V1=rnorm(300), V2=rnorm(300))
  y_ms <- y ## FIX THIS
  y_ms[] <- 2
  temp <- unmarkedFrameOccuMS(y_ms, siteCovs=sc)
  s <- expect_message(simulate(temp, psiformulas=stateformulas, detformulas=detformulas,
                               coefs=cf)[[1]])
  fm <- occuMS(detformulas, stateformulas, data=s, parameterization="multinomial")
  expect_equivalent(coef(fm, 'state'), c(-0.3121,0.8289,-0.4710,-0.8786), tol=1e-3)

  # gdistremoval
  set.seed(123)
  cf <- list(lambda=c(intercept=log(5), sc1=0.7), dist=c(intercept=log(50)),
           rem=c(intercept=log(0.2/(1-0.2)), oc1=0.4))
  sc <- data.frame(sc1=rnorm(200))
  oc <- data.frame(oc1=rnorm(200*5))
  temp <- unmarkedFrameGDR(yDistance=matrix(NA, 200, 4), yRemoval=matrix(NA, 200, 5),
                           siteCovs=sc, obsCovs=oc, dist.breaks=c(0,25,50,75,100), 
                           unitsIn='m')
  s <- expect_message(simulate(temp, lambdaformula=~sc1, removalformula=~oc1, distanceformula=~1,
                               output='abund', coefs=cf)[[1]])

  fm <- gdistremoval(lambdaformula=~sc1, removalformula=~oc1, distanceformula=~1,
                     output='abund', K=50, data=s)
  expect_is(fm, "unmarkedFitGDS")
  expect_equivalent(coef(fm, 'lambda'), c(1.4914,0.6454), tol=1e-4)

})

test_that("Random effects are re-generated for each simulated datset", {

  set.seed(123)

  M <- 100
  J <- 3
  x <- rnorm(M)
  group <- sample(letters[1:10], M, replace=TRUE)
  coefs <- list(state = c(0, 0.3, exp(0.5)), det = 0)
  y <- matrix(NA, M, J)

  umf <- unmarkedFrameOccu(y, siteCovs=data.frame(x=x, group=factor(group)))

  sims <- simulate(umf, model = occu, formula = ~1~x+(1|group),
                   coefs = coefs, nsim=30)

  fits <- lapply(sims, function(x) occu(~1~x+(1|group), x))

  sigs <- sapply(fits, function(x) sigma(x)$sigma)
  ints <- sapply(fits, function(x) coef(x)[1])
  
  # Check that true values are roughly in center of distribution
  # this was generally not true when random effects were only simulated once
  # and then used in every subsequent simulated dataset
  pct <- ecdf(sigs)
  expect_equal(pct(exp(0.5)), 0.47, tol=1e-2)
  pct <- ecdf(ints)
  expect_equal(pct(0), 0.5, tol=1e-2)
})

test_that("Error returned for unsupported model-name-based simulation",{
  expect_error(simulate("occu", formula=~1~1), "no longer supported")
})
