\name{unmarkedModSel-class}
\title{unmarkedModSel class and methods}
\alias{unmarkedModSel-class}
\alias{show,unmarkedModSel-method}
\alias{coef,unmarkedModSel-method}
\alias{SE,unmarkedModSel-method}

\description{The unmarkedModSel class and associated methods}

\section{Slots}{
  \describe{
	\item{Full}{data.frame with formula, estimates, standard errors and model 
		selection information. Converge is optim convergence code. CondNum is
		model condition number. n is the number of sites. delta is delta AIC.
		cumltvWt is cumulative AIC weight.
        Rsq is Nagelkerke's (1991) R-squared index, which is only returned when 
        the nullmod argument is specified. }
	\item{Names}{matrix referencing column names of estimates (row 1) and 
	   standard errors (row 2).}
  }
}

\section{Methods}{
  \describe{
    \item{show}{Print the AIC model selection table}
    \item{coef}{Data frame of coefficients from all models in model selection table}
    \item{SE}{Data frame of coefficient SEs from all models in model selection table}
  }
}

\seealso{\code{\link{modSel}}}

\examples{
data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 
lengths <- linetran$Length * 1000

ltUMF <- with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
	siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
	tlength = lengths, survey = "line", unitsIn = "m")
	})

fm1 <- distsamp(~ 1 ~1, ltUMF)
fm2 <- distsamp(~ area ~1, ltUMF)
fm3 <- distsamp( ~ 1 ~area, ltUMF)

fl <- fitList(Null=fm1, A.=fm2, .A=fm3)
fl

ms <- modSel(fl, nullmod="Null")
ms

coef(ms)                            # Estimates only
SE(ms)                              # Standard errors only
(toExport <- as(ms, "data.frame"))  # Everything

}
