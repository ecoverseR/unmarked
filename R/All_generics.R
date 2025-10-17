# This filename starts with All_ so it is loaded by R first

# Converting existing R functions to generics----------------------------------
setGeneric("coef")
setGeneric("confint")
setGeneric("fitted")
setGeneric("head")
setGeneric("hist")
setGeneric("logLik")
setGeneric("plot")
setGeneric("predict")
setGeneric("profile")
setGeneric("residuals")
setGeneric("simulate")
setGeneric("summary")
setGeneric("update")
setGeneric("vcov")


# Generics for exported unmarked-specific functions----------------------------

# Transform an object to its natural scale.
setGeneric("backTransform", function(obj, ...) {
  standardGeneric("backTransform")
})

# Best unbiased predictors
# Method in ranef.R
setGeneric("bup", function(object, stat=c("mean", "mode"), ...){
  standardGeneric("bup")
})

# Cross validation for unmarkedFit
# Methods in unmarkedCrossVal.R
setGeneric("crossVal", function(object,
    method=c("Kfold","holdout","leaveOneOut"), folds=10, holdoutPct=0.25,
    statistic=RMSE_MAE, ...) standardGeneric("crossVal"))

# Get availability probability matrix, sites x primary periods
# Only implemented for unmarkedFitIDS, in IDS.R
setGeneric("getAvail", function(object, ...) standardGeneric("getAvail"))

# Get false positive detection probability from unmarkedFitOccuFP object
# Method in unmarkedFit.R
setGeneric("getB", function(object, ...) standardGeneric("getB"))

# Get data (unmarkedFrame) from an unmarkedFit object
# Method in unmarkedFit.R
setGeneric("getData", function(object) standardGeneric("getData"))

# Get false positive detection probability from unmarkedFitOccuFP object
# Method in unmarkedFit.R
setGeneric("getFP", function(object, ...) standardGeneric("getFP"))

# Get observation duration
# Only used by occuCOP
# Method in occuCOP.R
setGeneric("getL", function(object) standardGeneric("getL"))

# Get a matrix of detection probabilities for each observation
# Method for unmarkedFit in getP.R
setGeneric("getP", function(object, ...) standardGeneric("getP"))

# Get observation data from unmarkedFrame
# Method in unmarkedFrame.R
setGeneric("getY", function(object) standardGeneric("getY"))

# Extract hessian matrix
setGeneric("hessian",	function(object) standardGeneric("hessian"))

# Compute linear combinations of parameters.
setGeneric("linearComb", function(obj, coefficients, ...) {
  standardGeneric("linearComb")
})

# Likelihood ratio test
setGeneric("LRT", function(m1, m2) standardGeneric("LRT"))

# Model selection
# Method in unmarkedFitList.R
setGeneric("modSel", function(object, ...) standardGeneric("modSel"))

# Nonparametric bootstrap
# Method in nonparboot.R
setGeneric("nonparboot", function(object, ...) standardGeneric("nonparboot"))

# Get number of sites from an unmarkedFrame
# Method in unmarkedFrame.R
setGeneric("numSites", function(object) standardGeneric("numSites"))

# Get number of observations from an unmarkedFrame
# Method in unmarkedFrame.R
setGeneric("numY", function(object) standardGeneric("numY"))

# Extract observation covs from an unmarkedFrame
# Method in unmarkedFrame.R
setGeneric("obsCovs", function(object,...) standardGeneric("obsCovs"))

# Assign value to obsCovs slot
# Method in unmarkedFrame.R
setGeneric("obsCovs<-", function(object, value) standardGeneric("obsCovs<-"))

# Get number of independent observations from an unmarkedFrame
# Method in unmarkedFrame.R
setGeneric("obsNum", function(object) standardGeneric("obsNum"))

# Extract obsToY matrix
# Method in unmarkedFrame.R
setGeneric("obsToY", function(object) standardGeneric("obsToY"))

# Assign value to obsToY slot
# Method in unmarkedFrame.R
setGeneric("obsToY<-", function(object, value) standardGeneric("obsToY<-"))

# Optimize penalize likelihood penalty
# Only used by occuMulti
# Method in occuMulti.R
setGeneric("optimizePenalty", 
  function(object, penalties=c(0,2^seq(-4,4)), k = 5, boot = 30, ...){
    standardGeneric("optimizePenalty")
})

# Parametric bootstrap
# Methods in parboot.R
setGeneric("parboot", function(object, ...) standardGeneric("parboot"))

# Plot covariate effects
# Method in plotEffects.R
setGeneric("plotEffects", function(object, ...) standardGeneric("plotEffects"))

# Get data for plotting effects
# Method in plotEffects.R
setGeneric("plotEffectsData", function(object, ...) standardGeneric("plotEffectsData"))

# Get posterior samples from ranef
# Method in posteriorSamples.R
setGeneric("posteriorSamples", function(object, nsims, ...){
  standardGeneric("posteriorSamples")
})

# Power analysis
# Method in powerAnalysis.R
setGeneric("powerAnalysis", function(object, ...) standardGeneric("powerAnalysis"))

# Get projected values from an unmarkedFitColExt object
# Method in unmarkedFit.R
setGeneric("projected", function(object, mean=TRUE) standardGeneric("projected"))

# Get data frame of random effect estimates
# Method in mixedModelTools.R
setGeneric("randomTerms", function(object, ...) standardGeneric("randomTerms"))

# Get posteriors of latent effects
# Method in ranef.R
setGeneric("ranef", function(object, ...) standardGeneric("ranef"))

# Get richness
# Only used by occuComm
# Method in occuComm.R
setGeneric("richness", function(object, ...) standardGeneric("richness"))

# Get sample size (number of sites after dropping sites with missing values)
# Method in unmarkedFit.R
setGeneric("sampleSize", function(object) standardGeneric("sampleSize"))

# Get standard error values
setGeneric("SE", function(obj, ...) standardGeneric("SE"))

# Extract site covariates from an unmarkedFrame
# Method in unmarkedFrame.R
setGeneric("siteCovs", function(object,...) standardGeneric("siteCovs"))

# Assign value to siteCovs slot
# Method in unmarkedFrame.R
setGeneric("siteCovs<-", function(object, value) standardGeneric("siteCovs<-"))

# Get smoothed values from an unmarkedFitColExt object
# Method in unmarkedFit.R
setGeneric("smoothed", function(object, mean=TRUE) standardGeneric("smoothed"))

# Sum of squared errors statistic
# Methods in utils.R
setGeneric("SSE", function(fit, ...) standardGeneric("SSE"))

# Create list of unmarkedPower objects
# Method in power.R
setGeneric("unmarkedPowerList", function(object, ...){
  standardGeneric("unmarkedPowerList")
})

# Extract yearly site covs (primary period covs) from an unmarkedFrame
# Method in unmarkedFrame.R
setGeneric("yearlySiteCovs", function(object,...){
  standardGeneric("yearlySiteCovs")
})

# Assign value to yearlySiteCovs slot
# Method in unmarkedFrame.R
setGeneric("yearlySiteCovs<-",
	function(object, value) standardGeneric("yearlySiteCovs<-"))

# Generics for internal unmarked-specific functions----------------------------

# Internal method for fitted
# Methods for different unmarkedFit types in fitted.R
setGeneric("fitted_internal", function(object){
  standardGeneric("fitted_internal")
})

# Get y matrix from unmarkedFit for fitList
# Methods in unmarkedFitList.R
setGeneric("fl_getY", function(fit, ...) standardGeneric("fl_getY"))

# Create design matrices and handle missing values
# Methods for various unmarkedFrame types in getDesign.R
setGeneric("getDesign", function(umf, ...) standardGeneric("getDesign"))

# Get fitting function associated with an unmarkedFrame type
# Methods in simulate.R
setGeneric("get_fitting_function", function(object, model, ...){
  standardGeneric("get_fitting_function")
})

# Get correct formula for predict based on type
# Methods in predict.R
setGeneric("get_formula", function(object, type, ...){
  standardGeneric("get_formula")
})

# Internal method for getP
# Methods for different unmarkedFit types in getP.R
setGeneric("getP_internal", function(object) standardGeneric("getP_internal"))

# Internal method for getY
# Extracts observation data from an unmarkedFit object
# Methods for various unmarkedFit types in unmarkedFit.R
setGeneric("getY_internal", function(object) standardGeneric("getY_internal"))

# Get data frame of covariates for a given submodel type
# Methods in unmarkedFit.R
setGeneric("get_covariates", function(object, type, ...){
  standardGeneric("get_covariates")
})

# Get estimates from unmarkedEstimateList
# Method in unmarkedEstimate.R
setGeneric("estimates", function(object) standardGeneric("estimates"))

# Check if object has random effect
# Methods in mixedModelTools.R
setGeneric("has_random", function(object) standardGeneric("has_random"))

# Get maximum likelihood estimate from an unmarkedFit object
# Method in unmarkedFit.R
setGeneric("mle", function(object) standardGeneric("mle"))

# Get NLL function from an unmarkedFit object
# Method in unmarkedFit.R
setGeneric("nllFun", function(object) standardGeneric("nllFun"))

# Internal nonparboot method for specific unmarkedFit types
# Methods in nonparboot.R
setGeneric("nonparboot_internal", function(object, B, keepOldSamples){
  standardGeneric("nonparboot_internal")
})

# Internal nonparboot update method for specific unmarkedFit types
# Methods in nonparboot.R
setGeneric("nonparboot_update", function(object, data){
  standardGeneric("nonparboot_update")
})

# Run prediction calculations in chunks for speed
# Method in predict.R
setGeneric("predict_by_chunk", function(object, ...){
  standardGeneric("predict_by_chunk")
})

# Internal method for ranef
# Methods for specific unmarkedFit types in ranef.R
setGeneric("ranef_internal", function(object, ...) standardGeneric("ranef_internal"))

# Internal method for predict
# Methods for specific unmarkedFit types in predict.R
setGeneric("predict_internal", function(object, ...) standardGeneric("predict_internal"))

# Rebuild code used to call model
# Methods in update.R
setGeneric("rebuild_call", function(object) standardGeneric("rebuild_call"))

# Internal method for residuals
# Methods for specific unmarkedFit types in residuals.R
setGeneric("residuals_internal", function(object){
  standardGeneric("residuals_internal")
})

# Replace y-matrix in unmarkedFrame, used by parboot
# Methods in parboot.R
setGeneric("replaceY", function(object, newY, replNA = TRUE, ...){
  standardGeneric("replaceY")
})

# Create residual plot from an unmarkedFit object
# Method in unmarkedFit.R
setGeneric("residual_plot", function(x, ...) standardGeneric("residual_plot"))

# Internal method for simulate
# Methods for specific unmarkedFit types in simulate.R
setGeneric("simulate_internal", function(object, nsim) standardGeneric("simulate_internal"))

# Subset unmarkedFrame observations
# Method in square_brackets.R
setGeneric("subset_obs", function(umf, j) standardGeneric("subset_obs"))

# Subset unmarkedFrame sites
# Method in square_brackets.R
setGeneric("subset_sites", function(umf, i) standardGeneric("subset_sites"))

# Internal summary method
# Methods for different unmarkedFit types in unmarkedFit.R
setGeneric("summary_internal", function(object) standardGeneric("summary_internal"))

# Set all observation data values to 0
# Method in simulate.R
setGeneric("y_to_zeros", function(object, ...){
  standardGeneric("y_to_zeros")
})
