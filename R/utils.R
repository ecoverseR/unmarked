
genFixedNLL <- function(nll, whichFixed, fixedValues)
{
    function(params) {
        params[whichFixed] <- fixedValues
        do.call(nll, list(params))
        }
}

# nll the original negative log likelihood function
# MLE the full vector of MLE values
profileCI <- function(nll, whichPar, MLE, interval, level)
{
    stopifnot(length(whichPar) == 1)
    MLEnll <- nll(MLE)
    nPar <- length(MLE)
	chsq <- qchisq(level, 1)/2
    f <- function(value) {
        fixedNLL <- genFixedNLL(nll, whichPar, value)
            mleRestricted <- optim(MLE, fixedNLL)$value
        mleRestricted - MLEnll - chsq
        }
    lower <- tryCatch(uniroot(f, c(interval[1],MLE[whichPar]))$root,
        error = function(e) {
            warning("Lower endpoint of profile confidence interval is on the boundary.",
        call. = FALSE)
        -Inf
        })
    upper <- tryCatch(upper <- uniroot(f, c(MLE[whichPar], interval[2]))$root,
        error = function(e) {
            warning("Upper endpoint of profile confidence interval is on the boundary.",
        call. = FALSE)
        Inf
        })

    return(c(lower,upper))
}

## link functions and their gradients
logistic <- function(x) {
  1/(1 + exp(-x))
}


logistic.grad <- function(x) {
  exp(-x)/(exp(-x)+1)^2
}

#Complimentary log log link
cloglog <- function(x){
  1-exp(-exp(x))
}

cloglog.grad <- function(x){
  exp(-exp(x))
}

## use logarithms to vectorize row-wise products
## this speeds things up a LOT (vs. apply(x,1,prod))
rowProds <- function(x, na.rm = FALSE)
{
  exp(rowSums(log(x), na.rm = na.rm))
}

# This function converts an appropriatedly formated comma-separated
# values file (.csv) to a format usable by \emph{unmarked}'s fitting
# functions (see \emph{Details}).
csvToUMF <-
function(filename, long=FALSE, type, species = NULL, ...)
{

  dfin <- read.csv(filename, stringsAsFactors=TRUE)

  if(long == TRUE) return(formatLong(dfin, species, type = type, ...))
  else return(formatWide(dfin, type = type, ...))
}

# utility function to create a variable that follows the dates as 1,2,3,...
# site id is first column
# julian date is second column

dateToObs <- function(dfin)
{
    sitecol <- dfin[[1]]
    datecol <- dfin[[2]]

    # order by site, then obs date
    dfin <- dfin[order(sitecol,datecol),]
    sitecol <- dfin[[1]]
    datecol <- dfin[[2]]

    dTab <- table(datecol,sitecol)
    sites <- unique(sitecol)
    nSite <- length(sites)
    nStop <- colSums(dTab)
    nStop <- nStop[nStop > 0]  # get rid of stops for sites with no stops

    obsNum <- numeric(length(sitecol))
    # for each site i, replace stops with 1:nStop[i]
    for(i in 1:nSite){
        stops <- which(sitecol == sites[i])
        obsNum[stops] <- 1:nStop[i]
    }

    dfout <- cbind(dfin,obsNum)
    dfout
}

# take long data set and return data list
# column names must be
# site names, first
# date, one column
# response, one column
# obs vars, one per column
formatLong <- function(dfin, species = NULL, type, ...) {
  if (type %in% c("umarkedFrameMPois", "unmarkedFrameGMM"))
    stop("Multinomial data sets are not supported.")
  if(missing(type)) stop("type must be supplied")

  ## copy dates to last column so that they are also a covdata var
  nc <- ncol(dfin)
  dfin[[nc+1]] <- dfin[[2]]
  names(dfin)[nc+1] <- "JulianDate"

  if(!is.null(species)) {
    dfin$y <- ifelse(dfin$species == species, dfin$y, 0)
    dfin$y[is.na(dfin$y)] <- 0
    dfin$species = NULL
  }
# TODO: dbl check that multiple cells per site*time are handled correctly.
#  # sum up counts within time/site
#  expr <- substitute(recast(dfin[,1:3], sv + dv ~ ..., id.var = 1:2,
#                            fun.aggregate = sum),
#                     list(sv = as.name(names(dfin)[1]),
#                          dv = as.name(names(dfin)[2])))
#  dfin2 <- eval(expr)
#  dfin1 <- dfin[!duplicated(dfin[,1:2]),]
#
#  dfin <- merge(dfin1,dfin2, by = 1:2)
#  dfin[,3] <- dfin[,length(dfin)]
#  dfin <- dfin[,-length(dfin)]
    names(dfin)[3] <- "y"
    dfin <- dateToObs(dfin)

    #Create wide version of y with matching matrix of indices
    scol <- names(dfin)[1]
    df_sub <- dfin[c(scol, "y", "obsNum")]
    df_sub$yind <- 1:nrow(df_sub)
    ywide <- reshape(df_sub, idvar=scol, timevar="obsNum", direction="wide")

    y <- unname(as.matrix(ywide[grep("^y\\.", names(ywide))]))
    yind <- as.vector(t(as.matrix(ywide[grep("yind", names(ywide))])))

    #Reorder input data frame by y-index (=no factor issues)
    #Also drop site/time/y/numObs cols
    obsvars.df <- dfin[yind, -c(1:3, ncol(dfin)), drop=FALSE]
    rownames(obsvars.df) <- NULL

    ## check for siteCovs
    obsNum <- ncol(y)
    M <- nrow(y)
    site.inds <- matrix(1:(M*obsNum), M, obsNum, byrow = TRUE)
    siteCovs <- sapply(obsvars.df, function(x) {
      obsmat <- matrix(x, M, obsNum, byrow = TRUE)
      l.u <- apply(obsmat, 1, function(y) {
        row.u <- unique(y)
        length(row.u[!is.na(row.u)])
      })
      # if there are 0 or 1 unique vals per row, we have a sitecov
      if (all(l.u <= 1)) {
        u <- apply(obsmat, 1, function(y) {
          row.u <- unique(y)
          ## only remove NAs if there are some non-NAs.
          if(!all(is.na(row.u)))
            row.u <- row.u[!is.na(row.u)]
          row.u
        })
        u
      }
    })
    siteCovs <- as.data.frame(siteCovs[!sapply(siteCovs, is.null)],
                              stringsAsFactors=TRUE)
    if(nrow(siteCovs) == 0) siteCovs <- NULL

    ## remove sitecovs from obsvars
    obsvars.df <- obsvars.df[, !(names(obsvars.df) %in% names(siteCovs)), drop = FALSE]

    if (type %in% c("unmarkedFrameDS", "unmarkedFrameGDS"))
      # obsCovs cannot be used with distsamp
      do.call(type, list(y = y, siteCovs = siteCovs, ...))
    else
      do.call(type, list(y = y, siteCovs = siteCovs, obsCovs = obsvars.df, ...))
}

# column names must be
# site (optional, but if present, labeled "site")
# response: y.1, y.2, ..., y.J
# site vars: namefoo, namebar, ...
# obs v: namefoo.1, namefoo.2, ..., namefoo.J, namebar.1, .., namebar.J,..

formatWide <- function(dfin, sep = ".", obsToY, type, ...) {
  if (type %in% c("umarkedFrameMPois", "unmarkedFrameGMM", "double", "removal"))
    stop("Multinomial data sets are not supported.")

        # escape separater if it is regexp special
    reg.specials <- c('.', '\\', ':', '|', '(', ')', '[', '{', '^', '$',
                      '*', '+', '?')
    if(sep %in% reg.specials) {
        sep.reg <- paste("\\",sep,sep="")
    } else {
        sep.reg <- sep
    }

    dfnm <- colnames(dfin)

    y <- grep(paste("^y",sep.reg,"[[:digit:]]", sep=""),dfnm)
    J <- length(y)
    y <- as.matrix(dfin[,y])
    M <- nrow(y)

    if(identical(tolower(colnames(dfin))[1],"site")) {
        dfin <- dfin[,-1]
        dfnm <- dfnm[-1]
    }

    ncols <- length(dfnm)
    obsCovsPresent <- FALSE
    siteCovsPresent <- FALSE
    i <- J + 1
    while(i <= ncols) {     # loop through columns
        if(!identical(grep(paste(sep.reg,"[[:digit:]]+$",sep=""),
                           dfnm[i]),integer(0))) { # check if is obsdata
            newvar.name <- sub(paste(sep.reg,"[[:digit:]]+$",sep=""),'',
                               dfnm[i])
            newvar <- dfin[,grep(paste(newvar.name,sep.reg,
                                       "[[:digit:]]+$",sep=""),dfnm)]
            if(obsCovsPresent) {
                if(ncol(newvar) != R) {
                    stop("Not all observation-level covariates have the same number of columns.")
                } else {
                    obsCovs[newvar.name] <- as.vector(t(newvar))
                }
            } else {
                obsCovsPresent <- TRUE
                R <- ncol(newvar)
                obsCovs <- data.frame(newvar = as.vector(t(newvar)))
            }
            colnames(obsCovs)[length(obsCovs)] <- newvar.name
            i <- i + R
        }
        else {
            if(siteCovsPresent){
                siteCovs <- cbind(siteCovs,dfin[,i])
            } else {
                siteCovsPresent <- TRUE
                siteCovs <- data.frame(newvar = dfin[,i])
            }
            colnames(siteCovs)[length(siteCovs)] <- dfnm[i]
            i <- i + 1
        }
    }

    if(!obsCovsPresent) obsCovs <- NULL
    if(!siteCovsPresent) siteCovs <- NULL

    ## if don't know obsToY yet, use RxJ matrix of ones or diag if R == J
    if(missing(obsToY)) {
        if(identical(R,J)) {
            obsToY <- diag(J)
        } else {
            obsToY <- matrix(0, R, J)
        }
    }

    do.call(type, list(y = y, siteCovs = siteCovs, obsCovs = obsCovs, ...))
}


# This convenience function converts multi-year data in long format to
# unmarkedMultFrame Object.  See Details for more information.

formatMult <- function(df.in)
{
  .Deprecated("formatMult", package=NULL, 
              msg = paste("formatMult will be deprecated in the next version."),
             old = as.character(sys.call(sys.parent()))[1L])
    years <- sort(unique(df.in[[1]]))
    nY <- length(years)
    df.obs <- list()
    nsamp <- numeric()
    maxsamp <- max(table(df.in[[1]], df.in[[2]])) # the maximum samples/yr
    for(t in 1:nY){
        df.t <- df.in[df.in[[1]] == years[t],] # subset for current year
        df.t <- df.t[,-1] # remove year column
        df.t <- dateToObs(df.t)
        nsamp <- max(df.t$obsNum)
        if(nsamp < maxsamp) {
            newrows <- df.t[1:(maxsamp - nsamp), ] # just a placeholder
            newrows[,"obsNum"] <- ((nsamp + 1) : maxsamp)
            newrows[,3 : (ncol(df.t) - 1)] <- NA
            df.t <- rbind(df.t, newrows)
        }
        df.obs <- rbind(df.obs,cbind(year = years[t],df.t))
    }

    names(df.obs)[4] <- "y"
    scol <- names(df.obs)[2]
    df_sub <- df.obs[c("year", scol, "y", "obsNum")]
    df_sub$yind <- 1:nrow(df_sub)

    ywide <- reshape(df_sub, idvar=c(scol,"year"), timevar="obsNum",
                     direction="wide")
    ywide <- reshape(ywide, idvar=scol, timevar="year", direction="wide")
    #Reshape goofs up the order
    ywide <- ywide[order(ywide[,1]),]

    y <- unname(as.matrix(ywide[grep("^y\\.", names(ywide))]))
    yind <- as.vector(t(as.matrix(ywide[grep("yind", names(ywide))])))

    #Reorder input data frame by y-index (=no factor issues)
    #Also drop site/time/y/numObs cols
    obsvars.df <- df.obs[yind, -c(1:2, 4, ncol(df.obs)), drop=FALSE]
    rownames(obsvars.df) <- NULL

    ## check for siteCovs
    obsNum <- ncol(y)
    M <- nrow(y)
    site.inds <- matrix(1:(M*obsNum), M, obsNum, byrow = TRUE)
    siteCovs <- sapply(obsvars.df, function(x) {
        obsmat <- matrix(x, M, obsNum, byrow = TRUE)
        l.u <- apply(obsmat, 1, function(y) {
            row.u <- unique(y)
            length(row.u[!is.na(row.u)])
        })
        ## if there are 0 or 1 unique vals per row, we have a sitecov
        if(all(l.u %in% 0:1)) {
            u <- apply(obsmat, 1, function(y) {
                row.u <- unique(y)
                ## only remove NAs if there are some non-NAs.
                if(!all(is.na(row.u)))
                    row.u <- row.u[!is.na(row.u)]
                row.u
            })
            u
        }
    })
    siteCovs <- as.data.frame(siteCovs[!sapply(siteCovs, is.null)],
                              stringsAsFactors=TRUE)
    if(nrow(siteCovs) == 0) siteCovs <- NULL

    ## only check non-sitecovs
    obsvars.df2 <- as.data.frame(obsvars.df[, !(names(obsvars.df) %in%
                                                names(siteCovs))])
    names(obsvars.df2) <- names(obsvars.df)[!(names(obsvars.df) %in%
                                              names(siteCovs))]

    yearlySiteCovs <- sapply(obsvars.df2, function(x) {
        obsmat <- matrix(x, M*nY, obsNum/nY, byrow = TRUE)
        l.u <- apply(obsmat, 1, function(y) {
            row.u <- unique(y)
            length(row.u[!is.na(row.u)])
        })
        ## if there are 0 or 1 unique vals per row, we have a sitecov
        if(all(l.u %in% 0:1)) {
            u <- apply(obsmat, 1, function(y) {
                row.u <- unique(y)
                ## only remove NAs if there are some non-NAs.
                if(!all(is.na(row.u)))
                    row.u <- row.u[!is.na(row.u)]
                row.u
            })
            u
        }
    })
    yearlySiteCovs <- as.data.frame(yearlySiteCovs[!sapply(yearlySiteCovs,is.null)],
                                    stringsAsFactors=TRUE)
    if(nrow(yearlySiteCovs) == 0) yearlySiteCovs <- NULL

    # Extract siteCovs and yearlySiteCovs from obsvars
    finalobsvars.df <- as.data.frame(obsvars.df[, !(names(obsvars.df) %in%
                                                      c(names(siteCovs),
                                                        names(yearlySiteCovs)))])
    names(finalobsvars.df) <- names(obsvars.df)[!(names(obsvars.df) %in%
                                                    c(names(siteCovs),
                                                      names(yearlySiteCovs)))]

    umf <- unmarkedMultFrame(y = y, siteCovs = siteCovs,
                             obsCovs = finalobsvars.df, yearlySiteCovs =
                             yearlySiteCovs,
                             numPrimary = nY)
    return(umf)
}


truncateToBinary <- function(y) {
    if(max(y, na.rm = TRUE) > 1) {
        y <- ifelse(y > 0, 1, 0)
        warning("Some observations were > 1.  These were truncated to 1.")
    }
    return(y)
}


lambda2psi <- function(lambda)
{
.Deprecated("lambda2psi", package=NULL, 
              msg = paste("lambda2psi will be deprecated in the next version."),
             old = as.character(sys.call(sys.parent()))[1L])
if(any(lambda < 0))
    stop("lambda must be >= 0")
as.numeric(1 - exp(-lambda))
}


# Convert individual-level distance data to the
# transect-level format required by distsamp()

formatDistData <- function (distData, distCol, transectNameCol, dist.breaks, occasionCol,effortMatrix)
{
  .Deprecated("formatDistData", package=NULL, 
              msg = paste("formatDistData will be deprecated in the next version."),
             old = as.character(sys.call(sys.parent()))[1L])
  if (!is.numeric(distData[, distCol]))
    stop("The distances must be numeric")
  transects <- distData[, transectNameCol]
  if (!is.factor(transects)) {
    transects <- as.factor(transects)
    warning("The transects were converted to a factor")
  }

  if (missing(occasionCol)) {
    T <- 1
    occasions <- factor(rep(1, nrow(distData)))
  }
  else {
    occasions <- distData[, occasionCol]
    if (!is.factor(occasions)) {
      occasions <- as.factor(occasions)
      warning("The occasions were converted to a factor")
    }
    T <- nlevels(occasions)
  }
  M <- nlevels(transects)
  J <- length(dist.breaks) - 1
  if (missing(effortMatrix)) {
    effortMatrix <- matrix(nrow=M,ncol=T,1)
  }

  if (!is.numeric(effortMatrix)){
    stop("effortMatrix is not numeric")
    effortMatrix <- matrix(nrow=M,ncol=T,1)
  }


  dist.classes <- levels(cut(distData[, distCol], dist.breaks,
                             include.lowest = TRUE))
  ya <- array(NA, c(M, J, T), dimnames = list(levels(transects),
                                              dist.classes, paste("rep", 1:T, sep = "")))
  transect.levels <- levels(transects)
  occasion.levels <- levels(occasions)
  for (i in 1:M) {
    for (t in 1:T) {
      sub <- distData[transects == transect.levels[i] &
                        occasions == occasion.levels[t], , drop = FALSE]
      ya[i, , t] <- table(cut(sub[, distCol], dist.breaks,
                              include.lowest = TRUE))
    }
  }
  y <- matrix(ya, nrow = M, ncol = J * T)
  # takes into account the effortMatrix to allow for the insertion of NAs instead of 0s for surveys which were not completed
  ee <- array(NA, c(M,length(occasion.levels)*(length(dist.breaks)-1)))
  for(i in 1:length(occasion.levels)){
    ee[,((ncol(ee)/length(occasion.levels)*(i-1)+1):(ncol(ee)/length(occasion.levels)*i))] <- matrix(effortMatrix[,i], ncol=J, nrow=M)
  }
  ee[ee==0] <- NA
  y <- y * ee
  dn <- dimnames(ya)
  rownames(y) <- dn[[1]]
  if (T == 1)
    colnames(y) <- dn[[2]]
  else colnames(y) <- paste(rep(dn[[2]], times = T), rep(1:T, each = J), sep = "")
  return(y)
}


#Sum of squared errors method
setMethod("SSE", "unmarkedFit", function(fit, ...){
    sse <- sum(residuals(fit)^2, na.rm=TRUE)
    return(c(SSE=sse))
})

setMethod("SSE", "unmarkedFitOccuMulti", function(fit, ...){
    r <- do.call(rbind, residuals(fit))
    return(c(SSE = sum(r^2, na.rm=T)))
})

# Generate zero-inflated Poisson

rzip <- function(n, lambda, psi) {
    x <- rpois(n, lambda)
    x[runif(n) < psi] <- 0
    x
}

#Inverts Hessian. Returns blank matrix with a warning on a failure.
invertHessian <- function(optimOut, nparam, SE){

  blankMat <- matrix(NA, nparam, nparam)
  if(!SE) return(blankMat)

  tryCatch(solve(optimOut$hessian),
    error=function(e){
      warning("Hessian is singular. Try providing starting values or using fewer covariates.", call.=FALSE)
      return(blankMat)
  })
}

#Get u and a from distance sampling data
getUA <- function(umf){

  M <- numSites(umf)
  if(inherits(umf, "unmarkedFrameGDR")){
    J <- ncol(umf@yDistance) / umf@numPrimary
  } else if(methods::.hasSlot(umf, "numPrimary")){
    J <- ncol(getY(umf)) / umf@numPrimary
  } else {
    J <- ncol(getY(umf))
  }
  db <- umf@dist.breaks
  w <- diff(db)

  u <- a <- matrix(NA, M, J)
    switch(umf@survey,
    line = {
        for(i in 1:M) {
            a[i,] <- umf@tlength[i] * w
            u[i,] <- a[i,] / sum(a[i,])
            }
        },
    point = {
        for(i in 1:M) {
            a[i, 1] <- pi*db[2]^2
            if(J > 1){
              for(j in 2:J){
                a[i, j] <- pi*db[j+1]^2 - sum(a[i, 1:(j-1)])
              }
            }
            u[i,] <- a[i,] / sum(a[i,])
            }
        })
  list(a=a, u=u)

}

# Distance probability
# Consolidate this with functions below?
get_dist_prob <- function(object){
  db <- object@data@dist.breaks
  M <- numSites(object@data)
  T <- 1
  if(methods::.hasSlot(object@data, "numPrimary")){
    T <- object@data@numPrimary
  }

  if(inherits("object", "unmarkedFitGDR")){
    J <- ncol(object@data@yDist) / T
    type <- "dist"
  } else {
    J <- ncol(object@data@y) / T
    type <- "det"
  }
  stopifnot(J == length(db) - 1)

  par1 <- matrix(0, M, T)
  if(object@keyfun != "uniform"){
    par1 <- predict(object, type=type, level=NULL, na.rm=FALSE)$Predicted
    par1 <- matrix(par1, M, T, byrow=TRUE)
  }

  scale <- 0.0
  if(object@keyfun == "hazard"){
    scale <- exp(coef(object, type = "scale"))
  }

  ua <- getUA(object@data)
  u <- ua$u
  a <- ua$a
  w <- diff(db)

  cp <- array(NA, c(M,J,T))

  for (i in 1:M){
    for (t in 1:T){
      if(is.na(par1[i,t])) next    
      cp[i,,t] <- getDistCP(object@keyfun, par1[i,t], scale, object@data@survey,
                            db, w, a[i,], u[i,])
    }
  }
  matrix(cp, nrow=M)
}


# Get area for converting density to abundance in distance sampling models
get_ds_area <- function(umf, unitsOut){
  db <- umf@dist.breaks
 
  # Calculate area based on survey type
  A <- switch(umf@survey,
    line = umf@tlength * max(db) * 2,
    point = pi * max(db)^2
  )

  # Convert input unit to kmsq
  A <- switch(umf@unitsIn,
    m = A / 1e6,
    km = A
  )

  # Convert kmsq to output unit
  A <- switch(unitsOut,
    m = A * 1e6,
    ha = A * 100,
    kmsq = A
  )

  A
}

pHalfnorm <- function(sigma, survey, db, w, a){
  J <- length(w)
  cp <- rep(NA, J)
  switch(survey,
    line = {
      f.0 <- 2 * dnorm(0, 0, sd=sigma)
      int <- 2 * (pnorm(db[-1], 0, sd=sigma) - pnorm(db[-(J+1)], 0, sd=sigma))
      cp[1:J] <- int / f.0 / w
    },
    point = {
      for(j in 1:J) {
        cp[j] <- integrate(grhn, db[j], db[j+1],
                           sigma=sigma, rel.tol=1e-4)$value *
                           2 * pi / a[j]
      }
    }
  )
  cp
}

pExp <- function(rate, survey, db, w, a){
  J <- length(w)
  cp <- rep(NA, J)
  switch(survey,
    line = {
      for(j in 1:J) {
        cp[j] <- integrate(gxexp, db[j], db[j+1],
                           rate=rate, rel.tol=1e-4)$value / w[j]
      }
    },
    point = {
      for(j in 1:J) {
        cp[j] <- integrate(grexp, db[j], db[j+1],
                            rate=rate, rel.tol=1e-4)$value *
                            2 * pi / a[j]
      }
    }
  )
  cp
}

pHazard <- function(shape, scale, survey, db, w, a){
  J <- length(w)
  cp <- rep(NA, J)
  switch(survey,
    line = {
      for(j in 1:J) {
        cp[j] <- integrate(gxhaz, db[j], db[j+1],
                           shape=shape, scale=scale,
                           rel.tol=1e-4)$value / w[j]
      }
    },
    point = {
      for(j in 1:J) {
        cp[j] <- integrate(grhaz, db[j], db[j+1],
                           shape = shape, scale=scale,
                           rel.tol=1e-4)$value * 2 * pi / a[j]
      }
    })
    cp
}

getDistCP <- function(keyfun, param1, param2, survey, db, w, a, u){
  switch(keyfun,
    halfnorm = {
      cp <- pHalfnorm(param1, survey, db, w, a)
    },
    exp = {
      cp <- pExp(param1, survey, db, w, a)
    },
    hazard = {
      cp <- pHazard(param1, param2, survey, db, w, a)
    },
    uniform = {
      cp <- rep(1, length(u))
    })
    cp * u
}


#Modified rmultinom for handling NAs
rmultinom2 <- function(n, size, prob){
  if(is.na(size)|any(is.na(prob))){
    return(matrix(NA, length(prob), length(n)))
  }
  stats::rmultinom(n=n, size=size, prob=prob)

}

#Functions to convert character columns to factors
df_to_factor <- function(df, name="Input"){
  if(is.null(df)) return(NULL)
  stopifnot(inherits(df, "data.frame"))
  char_cols <- sapply(df, is.character)
  if(any(char_cols)){
    warning(paste(name, "contains characters. Converting them to factors."), call.=FALSE)
  }
  to_change <- df[,char_cols, drop=FALSE]
  df[,char_cols] <- lapply(to_change, as.factor)
  df
}

umf_to_factor <- function(umf){
  stopifnot(inherits(umf, "unmarkedFrame"))
  umf@siteCovs <- df_to_factor(siteCovs(umf), "siteCovs")
  umf@obsCovs <- df_to_factor(obsCovs(umf), "obsCovs")
  if(methods::.hasSlot(umf, "yearlySiteCovs")){
    umf@yearlySiteCovs <- df_to_factor(yearlySiteCovs(umf), "yearlySiteCovs")
  }
  umf
}

dzip <- function(x, lambda, psi) {
  stopifnot(length(lambda)==length(x) | length(lambda)==1)
  stopifnot(length(psi)==1)
  stopifnot(all(x>=0))
  if(length(lambda)==1) lambda <- rep(lambda, length(x))
  den <- rep(NA, length(x))
  zer <- x==0
  gr0 <- x > 0
  den[zer] <- psi + (1-psi)*exp(-lambda[zer])
  den[gr0] <- (1-psi)*dpois(x[gr0], lambda[gr0])
  den
}

sapply2 <- function(X, FUN, ..., cl = NULL){
  if(requireNamespace("pbapply", quietly=TRUE)){
    return(pbapply::pbsapply(X=X, FUN=FUN, ..., cl = cl))
  } else if(!is.null(cl)){
    return(parallel::parSapply(cl=cl, X=X, FUN=FUN, ...))
  }
  sapply(X=X, FUN=FUN, ...)
}

lapply2 <- function(X, FUN, ..., cl = NULL){
  if(requireNamespace("pbapply", quietly=TRUE)){
    return(pbapply::pblapply(X=X, FUN=FUN, ..., cl = cl))
  } else if(!is.null(cl)){
    return(parallel::parLapply(cl=cl, X=X, fun=FUN, ...))
  }
  lapply(X=X, FUN=FUN, ...)
}

# Determine automatic K or check provided K for multinomial-type models
# (gdistsamp, gmultmix, distsampOpen, multmixOpen, gdistremoval)
check_K_multinomial <- function(K, K_adjust = 0, y, T = 1){

  safe_sum <- function(x){
    if(all(is.na(x))) return(NA) else return(sum(x, na.rm=TRUE))
  }

  if(T == 1){
    yt <- apply(y, 1, safe_sum)
  } else {
    M <- nrow(y)
    J <- ncol(y) / T
    ya <- array(y, c(M, J, T))
    ya <- aperm(ya, c(1,3,2))
    yt <- apply(ya, 1:2, safe_sum)
  }
  Kmin <- max(yt, na.rm = TRUE)
  if(missing(K)){
    Kout <- Kmin + K_adjust
    warning("K was not specified and was set to ", Kout, ".", call.=FALSE)
  } else {
    if(K <= Kmin){
      stop("specified K is too small. Try a value larger than the max count at any site", call.=FALSE)
    }
    Kout <- K
  }
  Kout
}
