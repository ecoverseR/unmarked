setMethod("update", "unmarkedFit", function(object, ..., evaluate=TRUE){
  inputs <- list(...)
  nm <- names(inputs)
  if(any(nm == "")) stop("All provided arguments to update should be named", call.=FALSE)
  # Re-create the original call with as many original args as possible
  cl <- rebuild_call(object)

  for (i in names(inputs)){
    cl[[i]] <- NA # Need to create blank entry first to fill in
    if(i == "data"){
      # Insert the quoted object name instead of the entire umf object
      # This prevents the @call slot from being filled with data
      cl[[i]] <- substitute(INP, list(INP = quote(inputs$data)))
    } else {
      cl[[i]] <- inputs[[i]]
    }
  }

  if(!evaluate) return(cl)

  eval(cl)
})

# The rebuild_call method starts with a fitted model's saved call, and inserts
# as much saved information from the object into the call as possible.
# This avoids situations where the original call referenced e.g. a 
# vector of formulas saved in the global environment, which would not be
# available if updating the model later in a different environment.
# In general we're just inserting everything from fit-specific slots
# in the unmarkedFit object that is also an input argument.

# This method is not actually used directly by any fitting function, instead
# it's a base piece of the call for most of them.
# Note that data is saved as quoted code instead of directly as an object
# This prevents the entire contents of the unmarkedFrame from being
# inserted directly into the call which clogs up the console
setMethod("rebuild_call", "unmarkedFit", function(object){ 
  cl <- object@call
  cl[["data"]] <- quote(object@data)
  cl
})

setMethod("rebuild_call", "unmarkedFitColExt", function(object){ 
  cl <- methods::callNextMethod(object) # calls base class method above
  cl[["psiformula"]] <- object@formlist$psi
  cl[["gammaformula"]] <- object@formlist$col
  cl[["epsilonformula"]] <- object@formlist$ext
  cl[["pformula"]] <- object@formlist$det
  cl
})

# Covers MMO and PCO
setMethod("rebuild_call", "unmarkedFitDailMadsen", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["lambdaformula"]] <- object@formlist$lambda
  cl[["gammaformula"]] <- object@formlist$gamma
  cl[["omegaformula"]] <- object@formlist$omega
  cl[["pformula"]] <- object@formlist$det
  cl[["iotaformula"]] <- object@formlist$iota
  cl[["mixture"]] <- object@mixture
  cl[["dynamics"]] <- object@dynamics
  cl[["immigration"]] <- object@immigration
  cl[["fix"]] <- object@fix
  if(methods::.hasSlot(object, "K")){
    cl[["K"]] <- object@K
  }
  cl
})

setMethod("rebuild_call", "unmarkedFitDS", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["formula"]] <- object@formula
  cl[["keyfun"]] <- object@keyfun
  cl[["unitsOut"]] <- object@unitsOut
  cl[["output"]] <- object@output
  cl
})

setMethod("rebuild_call", "unmarkedFitDSO", function(object){ 
  cl <- methods::callNextMethod(object) # Calls DailMadsen method
  cl[["keyfun"]] <- object@keyfun
  cl[["unitsOut"]] <- object@unitsOut
  cl[["output"]] <- object@output
  cl
})

# Also covers unmarkedFitGMM
setMethod("rebuild_call", "unmarkedFitGPC", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["lambdaformula"]] <- object@formlist$lambda
  cl[["phiformula"]] <- object@formlist$phi
  cl[["pformula"]] <- object@formlist$det
  cl[["mixture"]] <- object@mixture
  if(methods::.hasSlot(object, "K")){
    cl[["K"]] <- object@K
  }
  cl
})

setMethod("rebuild_call", "unmarkedFitGDS", function(object){ 
  cl <- methods::callNextMethod(object) # calls GMM method
  cl[["keyfun"]] <- object@keyfun
  cl[["unitsOut"]] <- object@unitsOut
  cl[["output"]] <- object@output
  cl
})

setMethod("rebuild_call", "unmarkedFitMPois", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["formula"]] <- object@formula
  cl
})

setMethod("rebuild_call", "unmarkedFitNmixTTD", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["stateformula"]] <- object@formlist$state
  cl[["detformula"]] <- object@formlist$det
  cl[["K"]] <- object@K
  cl
})

setMethod("rebuild_call", "unmarkedFitOccu", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["formula"]] <- object@formula
  if(is.null(object@TMB)){
    cl[["knownOcc"]] <- quote(object@knownOcc)
  }
  cl
})

setMethod("rebuild_call", "unmarkedFitOccuFP", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["detformula"]] <- object@formlist$det
  cl[["FPformula"]] <- object@formlist$fp
  cl[["Bformula"]] <- object@formlist$b
  cl[["stateformula"]] <- object@formlist$state
  cl
})

setMethod("rebuild_call", "unmarkedFitOccuMS", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["detformulas"]] <- quote(object@formlist$det)
  cl[["psiformulas"]] <- quote(object@formlist$state)
  if(!all(is.na(object@formlist$transition))){
    cl[["phiformulas"]] <- quote(object@formlist$transition)
  }
  cl[["parameterization"]] <- object@parameterization
  cl
})

setMethod("rebuild_call", "unmarkedFitOccuMulti", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["stateformulas"]] <- quote(object@formlist$state)
  cl[["detformulas"]] <- quote(object@formlist$det)
  cl
})

setMethod("rebuild_call", "unmarkedFitPCount", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["formula"]] <- object@formula
  cl[["mixture"]] <- object@mixture
  if(methods::.hasSlot(object, "K")){
    cl[["K"]] <- object@K
  }
  cl
})

setMethod("rebuild_call", "unmarkedFitOccuPEN", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["formula"]] <- object@formula
  cl[["knownOcc"]] <- quote(object@knownOcc)
  cl[["pen.type"]] <- object@pen.type
  cl[["lambda"]] <- object@lambda
  cl
})

setMethod("rebuild_call", "unmarkedFitOccuRN", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["formula"]] <- object@formula
  # For backwards compatibility
  if(methods::.hasSlot(object, "K")){
    cl[["K"]] <- object@K
  }
  cl
})

setMethod("rebuild_call", "unmarkedFitOccuTTD", function(object){ 
  cl <- methods::callNextMethod(object)
  cl[["psiformula"]] <- object@formlist$psi
  cl[["gammaformula"]] <- object@formlist$col
  cl[["epsilonformula"]] <- object@formlist$ext
  cl[["detformula"]] <- object@formlist$det
  cl
})
