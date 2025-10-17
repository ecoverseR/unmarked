context("predict-related functions")

skip_on_cran()

set.seed(123)
umf <- unmarkedFrameOccu(
  y = matrix(NA, 100, 5),
  siteCovs = data.frame(elev = rnorm(100),
                        group = factor(sample(c("A","B","C"), 100, replace=TRUE)))
)
umf <- simulate(umf, model = occu, formula=~1~elev+group,
                coefs = list(state = c(0, 0.4, -0.5, 0.6), det = 0))[[1]]
mod <- occu(~1~elev+group, umf)

test_that("clean_up_covs works with dynamic model data",{

  # Dynamic data
  y <- matrix(c(
      3, 2, 1, 4,
      3, 4, 2, 1,
      0, 1, 2, 3
      ), 3, 4, byrow=TRUE)
  siteCovs <- data.frame(sc1 = 1:3)
  obsCovs <- data.frame(oc1 = 1:12)
  ysc <- data.frame(ysc1 = 1:6, ysc2=factor(rep(c("a","b"), 3)))
  #ysc <- data.frame(ysc1 = 1:6, ysc2=factor(rep(c("a","b"), each=3)))
  umf <- unmarkedFramePCO(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
                          yearlySiteCovs=ysc, numPrimary=2)

  dr <- unmarked:::clean_up_covs(umf, drop_final=TRUE)
  expect_equal(dr$site_covs, data.frame(sc1=1:3))
  expect_equal(dr$yearly_site_covs, data.frame(ysc1=c(1,NA,3,NA,5,NA),
                                                 ysc2=factor(c("a",NA,"a",NA,"a",NA)),
                                                 sc1=c(1,NA,2,NA,3,NA)))
  expect_equivalent(dr$obs_covs, data.frame(oc1=1:12, ysc1=rep(1:6, each=2),
                                       ysc2=factor(rep(c("a","b"), each=2)),
                                       sc1=rep(1:3, each=4),
                                       obsNum=factor(rep(1:4,3))))

  no_drop <- unmarked:::clean_up_covs(umf)
  expect_equivalent(no_drop$yearly_site_covs, data.frame(ysc1=1:6,
                                                    ysc2=factor(rep(c("a","b"),3)),
                                                    sc1=rep(1:3, each=2)))

  umf <- unmarkedFramePCO(y=y, numPrimary=2)

  cc <- unmarked:::clean_up_covs(umf, drop_final=TRUE)
  expect_equivalent(cc$obs_covs,
                    data.frame(.dummy3=rep(1,12), .dummy2=rep(1,12), .dummy1=rep(1,12),
                               obsNum=factor(rep(1:4, 3))))
})

test_that("clean_up_covs works with single-season models",{
  y <- matrix(c(0,1,0,0,1,1), nrow=3)
  umf <- unmarkedFrameOccu(y=y, siteCovs=data.frame(sc1=1:3),
                           obsCovs=data.frame(oc1=1:6))
  cc <- unmarked:::clean_up_covs(umf)
  expect_equal(names(cc), c("site_covs","obs_covs"))
  expect_equivalent(cc$site_covs, data.frame(sc1=1:3))
  expect_equivalent(cc$obs_covs, data.frame(oc1=1:6, sc1=rep(1:3, each=2),
                                            obsNum = factor(rep(1:2, 3))))
  cc2 <- unmarked:::clean_up_covs(umf, drop_final=TRUE)
  expect_equal(cc, cc2)
})

test_that("clean_up_covs works with models with no obs covs",{
  # single season
  ltUMF <- with(linetran, {
    unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
      siteCovs = data.frame(Length, area, habitat),
      dist.breaks = c(0, 5, 10, 15, 20),
      tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
    })
  ltUMF

  cc <- unmarked:::clean_up_covs(ltUMF)
  expect_equal(names(cc), c("site_covs", "obs_covs"))
  expect_equal(dim(cc$obs_covs), c(12,5))
})

test_that("clean_up_covs works with models where length(y) != length(p)",{
  # double observer, etc
  nSites <- 3
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }

  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
  expect_warning(umf <- unmarkedFrameMPois(y=y,
         siteCovs <- data.frame(sc1=1:3),
         obsCovs=list(observer=observer),
         type="double"))

  cc <- unmarked:::clean_up_covs(umf)
  expect_equivalent(cc$site_covs, data.frame(sc=1:3))
  expect_equivalent(cc$obs_covs, data.frame(observer=factor(c(rep(c("A","B"), 3))),
                                            sc1=rep(1:3, each=2),
                                            obsNum=factor(rep(1:2, 3))))

})

test_that("predicting from raster works",{

  skip_if(!requireNamespace("raster", quietly=TRUE), 
          "raster package unavailable")

  set.seed(123)
  # Create rasters
  # Elevation
  r_elev <- data.frame(x=rep(1:10, 10), y=rep(1:10, each=10), z=rnorm(100))
  r_elev <- raster::rasterFromXYZ(r_elev)

  #Group
  r_group <- data.frame(x=rep(1:10, 10), y=rep(1:10, each=10),
                      z=sample(1:length(levels(umf@siteCovs$group)), 100, replace=T))
  # Convert to 'factor' raster
  r_group <- raster::as.factor(raster::rasterFromXYZ(r_group))
  r_group@data@attributes <- data.frame(ID=raster::levels(r_group)[[1]], group=levels(umf@siteCovs$group))

  # Stack
  nd_raster <- raster::stack(r_elev, r_group)
  names(nd_raster) <- c("elev", "group")
  raster::crs(nd_raster) <- 32616

  pr <- predict(mod, 'state', newdata=nd_raster)
  expect_is(pr, 'RasterStack')
  expect_equal(names(pr), c("Predicted","SE","lower","upper"))
  expect_equal(pr[1,1][1], 0.3675313, tol=1e-5)
  expect_equal(raster::crs(pr), raster::crs(nd_raster))

  #append data
  pr <- predict(mod, 'state', newdata=nd_raster, appendData=TRUE)
  expect_is(pr, 'RasterStack')
  expect_equal(names(pr)[5:6], c("elev","group"))

  # Missing levels are handled
  nd_2 <- nd_raster[[1]]
  expect_error(predict(mod, 'state', newdata=nd_2))
})

test_that("predicting from terra::rast works",{

  skip_if(!requireNamespace("terra", quietly=TRUE), 
          "terra package unavailable")

  set.seed(123)
  # Create rasters
  # Elevation
  r_elev <- data.frame(x=rep(1:10, 10), y=rep(1:10, each=10), z=rnorm(100))
  r_elev <- terra::rast(r_elev, type="xyz")

  #Group
  r_group <- data.frame(x=rep(1:10, 10), y=rep(1:10, each=10),
                      z=sample(1:length(levels(umf@siteCovs$group)), 100, replace=T))
  # Convert to 'factor' raster
  r_group <- terra::as.factor(terra::rast(r_group, type="xyz"))
  levels(r_group) <- data.frame(ID=terra::levels(r_group)[[1]]$ID, group=levels(umf@siteCovs$group))

  # Stack
  nd_raster <- c(r_elev, r_group)
  names(nd_raster) <- c("elev", "group")
  terra::crs(nd_raster) <- "epsg:32616"

  pr <- predict(mod, 'state', newdata=nd_raster)
  expect_is(pr, 'SpatRaster')
  expect_equal(names(pr), c("Predicted","SE","lower","upper"))
  expect_equivalent(pr[1,1][1], 0.3675313, tol=1e-5)
  expect_equal(terra::crs(pr), terra::crs(nd_raster))

  #append data
  pr <- predict(mod, 'state', newdata=nd_raster, appendData=TRUE)
  expect_is(pr, 'SpatRaster')
  expect_equal(names(pr)[5:6], c("elev","group"))

  # Missing levels are handled
  nd_2 <- nd_raster[[1]]
  expect_error(predict(mod, 'state', newdata=nd_2))

  # Missing data in raster handled
  set.seed(123)
  # Create rasters
  # Elevation
  r_elev <- data.frame(x=rep(1:10, 10), y=rep(1:10, each=10), z=rnorm(100))
  r_elev$z[1] <- NA
  r_elev <- terra::rast(r_elev, type="xyz")
  expect_true(is.na(r_elev[91]))

  #Group
  r_group <- data.frame(x=rep(1:10, 10), y=rep(1:10, each=10),
                      z=sample(1:length(levels(umf@siteCovs$group)), 100, replace=T))
  # Convert to 'factor' raster
  r_group <- terra::as.factor(terra::rast(r_group, type="xyz"))
  levels(r_group) <- data.frame(ID=terra::levels(r_group)[[1]]$ID, group=levels(umf@siteCovs$group))

  # Stack
  nd_raster <- c(r_elev, r_group)
  names(nd_raster) <- c("elev", "group")
  terra::crs(nd_raster) <- "epsg:32616"

  pr_na <- predict(mod, 'state', newdata=nd_raster)
  expect_true(is.na(pr_na$Predicted[91]))

  # Random effect in formula
  mod2 <- update(mod, ~1 ~elev + (1|group))

  expect_is(predict(mod2, 'state', newdata=nd_raster), 'SpatRaster')
  expect_is(predict(mod2, 'state', newdata=nd_raster, re.form=NA), 'SpatRaster')
})

test_that("Warning when predicting and there are nested functions in formula", { 
  nd <- data.frame(elev=c(0,1), group="B")
  mod2 <- occu(~1~scale(elev)+group, umf)
  expect_no_warning(predict(mod2, type='state', newdata=nd))

  mod3 <- occu(~1~I(elev^2)+group, umf)
  expect_no_warning(predict(mod3, type='state', newdata=nd))

  mod4 <- occu(~1~I(scale(elev)^2)+group, umf)
  expect_warning(predict(mod4, type='state', newdata=nd), "probably incorrect")

  mod5 <- occu(~1~scale(I(elev^2))+group, umf)
  expect_warning(predict(mod5, type='state', newdata=nd), "probably incorrect")
})
