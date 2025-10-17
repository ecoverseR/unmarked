# This file starts with All_ so it is loaded first by R

# Class unions-----------------------------------------------------------------
setClassUnion("optionalDataFrame", c("data.frame","NULL"))
setClassUnion("optionalMatrix", c("matrix","NULL"))
setClassUnion("optionalNumeric", c("numeric","NULL"))
setClassUnion("optionalCharacter", c("character","NULL"))
setClassUnion("optionalList", c("list","NULL"))
setClassUnion("numericOrLogical", c("numeric", "logical"))
setClassUnion("matrixOrVector", c("matrix","numeric"))

# unmarkedFrame classes--------------------------------------------------------

# Basic unmarkedFrame
validunmarkedFrame <- function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    if(!is.null(object@siteCovs))
        if(nrow(object@siteCovs) != M)
            errors <- c(errors,
               "siteCovData does not have same size number of sites as y.")
    if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
        if(nrow(object@obsCovs) != M*obsNum(object))
            errors <- c(errors, "obsCovData does not have M*obsNum rows.")
    if(length(errors) == 0)
        TRUE
    else
        errors
}

setClass("unmarkedFrame",
    representation(y = "matrix",
        obsCovs = "optionalDataFrame",
        siteCovs = "optionalDataFrame",
        obsToY = "optionalMatrix"),
    validity = validunmarkedFrame)

### Single-season models ###

# Count-based occupancy
setClass("unmarkedFrameOccuCOP",
  representation(L = "matrix"),
  contains = "unmarkedFrame",
  validity = function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    y_integers = sapply(object@y, check.integer, na.ignore = T)
    if (!all(y_integers)) {
      errors <- c(errors,
                  paste(
                    "Count detection should be integers. Non-integer values:",
                    paste(object@y[which(!y_integers)], collapse = ', ')
                  ))
    }
    if (!all(all(dim(object@L) == dim(object@y)))){
      errors <- c( errors, paste(
        "L should be a matrix of the same dimension as y, with M =", M,
        "rows (sites) and J =", J, "columns (sampling occasions)."
      ))}
    if (length(errors) == 0) TRUE
    else errors
  }
)

# Distance sampling
setClass("unmarkedFrameDS",
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrame",
    validity = function(object) {
        errors <- character(0)
        J <- numY(object)
        db <- object@dist.breaks
        if(J != length(db) - 1)
            errors <- c(errors, "ncol(y) must equal length(dist.breaks)-1")
        if(db[1] != 0)
            errors <- c(errors, "dist.breaks[1] must equal 0")
        if(!is.null(obsCovs(object)))
            "obsCovs cannot be used with distsamp"
        if(length(errors) == 0) TRUE
        else errors
        })

# Multinomial mixture model
setClass("unmarkedFrameMPois",
	representation(samplingMethod = "character", piFun = "character"),
	contains = "unmarkedFrame"
)

# Single-season occupancy model
setClass("unmarkedFrameOccu", contains = "unmarkedFrame")

# Community occupancy model
setClass("unmarkedFrameOccuComm",
  representation(ylist = "list", speciesCovs="optionalList"),
  contains = "unmarkedFrame"
)

# False-positive occupancy model
setClass("unmarkedFrameOccuFP",
  representation(type = "numeric"),
  contains = "unmarkedFrame"
)

# Rota multiple species occupancy model
setClass("unmarkedFrameOccuMulti",
  representation(ylist = "list", fDesign = "matrix"),
  contains = "unmarkedFrame",
  validity = function(object) {
      errors <- character(0)
      M <- nrow(object@y)
      J <- ncol(object@y)
      Ms <- sapply(object@ylist,nrow)
      Js <- sapply(object@ylist,ncol)
      if(length(unique(Ms)) != 1)
        errors <- c(errors, "All species must have same number of sites")
      if(length(unique(Js)) != 1)
        errors <- c(errors, "All species must have same number of observations")
      if(!is.null(object@siteCovs))
        if(nrow(object@siteCovs) != M)
            errors <- c(errors,
               "siteCovData does not have same size number of sites as y.")
      if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
        if(nrow(object@obsCovs) != M*obsNum(object))
            errors <- c(errors, "obsCovData does not have M*obsNum rows.")
      if(length(errors) == 0)
        TRUE
      else
        errors
  }
)

# Twining multispecies occupancy-abundance model
setClass("unmarkedFrameOccuRNMulti",
         contains="unmarkedFrameOccuMulti")

# N-mixture model
setClass("unmarkedFramePCount", contains = "unmarkedFrame")


### Multi-season models ###

# Basic multi-season class used by dynamic occupancy model
setClass("unmarkedMultFrame",
  representation(numPrimary = "numeric", yearlySiteCovs = "optionalDataFrame"),
  contains="unmarkedFrame"
)

# Multi-state occupancy
setClass('unmarkedFrameOccuMS',
  representation(numStates = "numeric", phiOrder = "list"),
  contains = "unmarkedMultFrame"
)

# Time-to-detection occupancy
setClass("unmarkedFrameOccuTTD",
         representation(
            surveyLength = "matrix"),
          contains = "unmarkedMultFrame")


### Temporary emigration (TE) models

# Basic class for temporary emigration data
setClass("unmarkedFrameG3", contains = "unmarkedMultFrame")

# TE multinomial mixture
setClass("unmarkedFrameGMM",
  representation(piFun = "character", samplingMethod = "character"),
  contains = "unmarkedFrameG3"
)

# TE distance sampling
setClass("unmarkedFrameGDS",
  representation(dist.breaks = "numeric", tlength = "numeric",
    survey = "character", unitsIn = "character"),
  contains = "unmarkedFrameG3"
)

# TE N-mixture
setClass("unmarkedFrameGPC", contains = "unmarkedFrameG3")

# TE occupancy
setClass("unmarkedFrameGOccu", contains = "unmarkedFrameG3")

# TE distance-removal model
setClass("unmarkedFrameGDR",
  representation(
    yDistance = "matrix",
    yRemoval = "matrix",
    survey = "character",
    dist.breaks = "numeric",
    unitsIn = "character",
    period.lengths = "numeric"
  ),
  contains="unmarkedFrameG3"
)


### Open population count models ###

# Basic class for all open pop models
setClass("unmarkedFrameDailMadsen",
  representation(primaryPeriod = "matrix"),
  contains = "unmarkedMultFrame"
)

# Open-pop distance sampling
setClass("unmarkedFrameDSO",
  representation(dist.breaks = "numeric", tlength = "numeric",
                 survey = "character", unitsIn = "character"),
  contains = "unmarkedFrameDailMadsen"
)

# Open-pop multinomial mixture
setClass("unmarkedFrameMMO",
  representation(piFun = "character", samplingMethod = "character"),
  contains = "unmarkedFrameDailMadsen"
)

# Open-pop N-mixture
setClass("unmarkedFramePCO", contains = "unmarkedFrameDailMadsen")


# unmarkedEstimate class-------------------------------------------------------

# Stores parameter estimates for a submodel (e.g. state, det)
setClass("unmarkedEstimate",
  representation(
    name = "character",
		short.name = "character",
    estimates = "numeric",
    covMat = "matrix",
    fixed = "numeric",
    covMatBS = "optionalMatrix",
    invlink = "character",
    invlinkGrad = "character",
    randomVarInfo= "list"),
  validity = function(object){
    errors <- character(0)
    if(nrow(object@covMat) != length(object@estimates)) {
    errors <- c(errors,
        "Size of covMat does not match length of estimates.")
    }
    if(length(errors) > 0)
        errors
    else
        TRUE
  }
)

# List of unmarkedEstimates
setClass("unmarkedEstimateList",
  representation(estimates = "list"),
  validity = function(object) {
    errors <- character(0)
    for(est in object@estimates) {
      if(!is(est, "unmarkedEstimate")) {
        errors <- c("At least one element of unmarkedEstimateList is not an unmarkedEstimate.")
        break
      }
    }
    if(length(errors) == 0) {
      return(TRUE)
    } else {
      return(errors)
    }
  }
)

# unmarkedFit classes----------------------------------------------------------

# Basic class
setClass("unmarkedFit",
  representation(fitType = "character",
    call = "call",
    formlist = "list",
    data = "unmarkedFrame",
    sitesRemoved = "numeric",  # vector of indices of removed sites
    estimates = "unmarkedEstimateList",
    AIC = "numeric",
    opt = "list",
    negLogLike = "numeric",
    nllFun = "function",
    bootstrapSamples = "optionalList",
    covMatBS = "optionalMatrix", # list of bootstrap sample fits
    TMB = "optionalList")  #TMB output object
)

### Single-season model types ###

# Distance sampling
setClass("unmarkedFitDS",
  representation(formula = "formula", keyfun = "character", unitsOut = "character", output = "character"), 
  contains = "unmarkedFit"
)

# Multinomial mixture model
setClass("unmarkedFitMPois", 
  representation(formula = "formula"),
  contains = "unmarkedFit"
)

# N-mixture TTD model
setClass("unmarkedFitNmixTTD",
  representation(K = "numeric"),
  contains = "unmarkedFit"
)

# Basic occupancy model
setClass("unmarkedFitOccu",
  representation(formula = "formula", knownOcc = "logical"),
  contains = "unmarkedFit"
)

# Community occupancy model
setClass("unmarkedFitOccuComm", contains="unmarkedFitOccu")

# Count-based occupancy model
setClass("unmarkedFitOccuCOP",
  representation(removed_obs = "matrix"),
  contains = "unmarkedFit"
)

# False positive occupancy
setClass("unmarkedFitOccuFP",
  representation(knownOcc = "logical", type = "numeric"),
  contains = "unmarkedFit"
)

# Penalized likelihood occupancy model
setClass("unmarkedFitOccuPEN",
  representation(formula = "formula", knownOcc = "logical", pen.type = "character", lambda = "numeric"),
  contains = "unmarkedFit"
)

setClass("unmarkedFitOccuPEN_CV",
  representation(formula = "formula", knownOcc = "logical", pen.type = "character", 
                 lambdaVec = "numeric", k = "numeric", foldAssignments = "numeric", 
                 lambdaScores = "numeric", chosenLambda = "numeric"),
  contains = "unmarkedFit"
)

# Rota multispecies occupancy model
setClass("unmarkedFitOccuMulti", contains = "unmarkedFit")

# Royle-Nichols occupancy model
setClass("unmarkedFitOccuRN",
  representation(formula = "formula", K = "numeric"),
  contains = "unmarkedFit"
)

# Twining multispecies occupancy-abundance model
setClass("unmarkedFitOccuRNMulti",
         representation(
            detformulas = "list",
            stateformulas = "list",
            modelOccupancy = "numeric",
            K = "numeric"),
         contains = "unmarkedFit")

# N-mixture model
setClass("unmarkedFitPCount",
  representation(formula = "formula", K = "numeric", mixture = "character"),
  contains = "unmarkedFit"
)


### Temporary emigration (TE) model types

setClass("unmarkedFitG3", contains = "unmarkedFit")

# TE occupancy
setClass("unmarkedFitGOccu", contains = "unmarkedFitG3")

# TE N-mixture
setClass("unmarkedFitGPC",
  representation(mixture = "character", K = "numeric"),
  contains = "unmarkedFitG3"
)

# TE multinomial N-mixture model
setClass("unmarkedFitGMM", contains = "unmarkedFitGPC")

# TE distance sampling
setClass("unmarkedFitGDS",
  representation(keyfun = "character", unitsOut = "character", output = "character"),
  contains = "unmarkedFitGMM"
)

# TE distance-removal model
setClass("unmarkedFitGDR", contains = "unmarkedFitGDS")


### Open-population models ###

# Dynamic occupancy model
setClass("unmarkedFitColExt",
  representation(phi = "matrix", projected = "array",
                 projected.mean = "matrix", smoothed = "array", smoothed.mean = "matrix",
                 projected.mean.bsse = "optionalMatrix", smoothed.mean.bsse = "optionalMatrix"),
  contains = "unmarkedFit"
)

# Multi-state occupancy model
setClass("unmarkedFitOccuMS",
  representation(parameterization = "character"),
  contains = "unmarkedFit"
)

# Time-to-detection occupancy model
setClass("unmarkedFitOccuTTD", contains = "unmarkedFit")

# This class is not used directly, just used as a base for for PCO, MMO, DSO
setClass("unmarkedFitDailMadsen",
  representation(K = "numeric", mixture = "character", dynamics = "character",
                 immigration = "logical", fix = "character"),
  contains = "unmarkedFit"
)

# Open-pop distance sampling
setClass("unmarkedFitDSO",
  representation(keyfun = "character", unitsOut = "character", output = "character"),
  contains = "unmarkedFitDailMadsen"
)

# Open-pop multinomial mixture
setClass("unmarkedFitMMO", contains = "unmarkedFitDailMadsen")

# Open-pop N-mixture
setClass("unmarkedFitPCO", contains = "unmarkedFitDailMadsen")


### Integrated models ###

setClassUnion("unmarkedFrameOrNULL", members=c("unmarkedFrame", "NULL"))
setClass("unmarkedFitIDS",
    representation(
        keyfun = "character",
        K = "numeric",
        dataPC = "unmarkedFrameOrNULL",
        dataOC = "unmarkedFrameOrNULL",
        maxDist = "list",
        surveyDurations = "list",
        unitsOut = "character"),
        contains = "unmarkedFit")


# Other output objects---------------------------------------------------------

# Linear combinations object
setClass("unmarkedLinComb",
  representation(parentEstimate = "unmarkedEstimate", estimate = "numeric",
                 covMat = "matrix", covMatBS = "optionalMatrix", coefficients = "matrix")
)

# Back-transformed linear comb
setClass("unmarkedBackTrans",
  representation(parentLinComb = "unmarkedLinComb", estimate = "numeric",
                 covMat = "matrix", covMatBS = "optionalMatrix")
)

# Union of linear comb and back transform object types
setClassUnion("linCombOrBackTrans", c("unmarkedLinComb", "unmarkedBackTrans"))

# Profile likelihood output
setClass("profile", representation(prof = "matrix"))

# ranef output (posterior of latent variables)
setClass("unmarkedRanef", representation(post = "array"))

# Posterior samples of latent variable
setClass("unmarkedPostSamples",
  representation(numSites="numeric", numPrimary="numeric",
                 nsims="numeric", samples="array")
)

# List of unmarkedFit objects
setClass("unmarkedFitList",
  representation(fits = "list"),
  validity = function(object) {
    fl <- object@fits
    umf1 <- getData(fl[[1]])
    y1 <- fl_getY(fl[[1]])
    dataTest <- sapply(fl, function(x) isTRUE(all.equal(umf1, getData(x))))
    yTest <- sapply(fl, function(x) isTRUE(all.equal(y1, fl_getY(x))))
    if(!all(dataTest)) {
      stop("Data are not the same among models. Make sure you use the same unmarkedFrame object for all models.")
    } else if(!all(yTest)) {
      stop("Data are not the same among models due to missing covariate values. Consider removing NAs before analysis.")
    }
    TRUE
  }
)

# Output of model selection (modSel function)
setClass("unmarkedModSel",
  representation(Full = "data.frame", Names = "matrix")
)

# Parametric bootstrap output
setClass("parboot",
  representation(call = "call", t0 = "numeric", t.star = "matrix")
)

# Cross validation output
setClass("unmarkedCrossVal",
  representation(stats = "data.frame", summary = "data.frame", method = "character",
                 folds = "numeric", holdoutPct = "numeric"),
  validity=function(object){
    errors <- character(0)
    hp <- object@holdoutPct
    if(hp<0|hp>1){
      errors <- c(errors,"holdoutPct must be between 0 and 1")
    }
  }
)

# List of cross validation objects
setClass("unmarkedCrossValList",
  representation(stats_list="list", method = "character", folds="numeric",
                 holdoutPct="numeric", sort="character")
)

# Power analysis output
setClass("unmarkedPower",
  representation(call="call", data="unmarkedFrame", M="numeric",
                 J="numeric", T="numeric", coefs="list", estimates="list",
                 alpha="numeric", nulls="list")
)

# List of power analysis outputs
setClass("unmarkedPowerList", representation(powerAnalyses="list"))
