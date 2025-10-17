setMethod("simulate", "character", function(object, ...){
  stop("Supplying the name of an unmarked fitting function to simulate is no longer supported.\n",
       "Use the simulate method for an unmarkedFrame instead. See the simulate vignette for more.")
})

setMethod("powerAnalysis", "unmarkedFit", function(object, ...){  
  stop("Using an unmarkedFit object to run a power analysis is no longer supported.\n",
       "Use an unmarkedFrame instead. See the power analysis vignette for more.")
})

shinyPower <- function(object, ...){
  stop("shinyPower is no longer supported; we plan to add a new shiny app in the future.\n",
       "See the power analysis vignette for an alternative approach to power analysis.")
}
