# R package unmarked

<!-- badges: start -->
[![R build status](https://github.com/ecoverseR/unmarked/workflows/R-CMD-check/badge.svg)](https://github.com/ecoverseR/unmarked/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/unmarked)](https://cran.r-project.org/package=unmarked)
<!-- badges: end -->

`unmarked` is an [R](https://www.r-project.org/) package for analyzing ecological data arising from several popular sampling techniques.  The sampling methods include point counts, occurrence sampling, distance sampling, removal, double observer, and many others. `unmarked` uses hierarchical models to incorporate covariates of the latent abundance (or occupancy) and imperfect detection processes.

## Installation

The latest stable version of unmarked can be downloaded from [CRAN](https://cran.r-project.org/package=unmarked):

```r
install.packages("unmarked")
```

The latest development version can be installed from Github:

```r
install.packages("remotes")
remotes::install_github("ecoverseR/unmarked")
```

## Support

Support is provided through the [unmarked Google group](http://groups.google.com/group/unmarked).
The package [website](https://ecoverseR.github.io/unmarked) has more information.
You can report bugs [here](https://github.com/ecoverseR/unmarked/issues), by posting to the Google group, or by emailing [the current maintainer](https://kenkellner.com).

## Get Started

See the following vignettes for an introduction to `unmarked` and some example analyses:

[Overview of `unmarked`](https://cran.r-project.org/web/packages/unmarked/vignettes/unmarked.html)

[Dynamic occupancy models](https://cran.r-project.org/web/packages/unmarked/vignettes/colext.html)

[Multispecies occupancy models](https://cran.r-project.org/web/packages/unmarked/vignettes/occuMulti.html)

[Distance sampling](https://cran.r-project.org/web/packages/unmarked/vignettes/distsamp.html)

[Contributing to `unmarked`](https://cran.r-project.org/web/packages/unmarked/vignettes/contributing_to_unmarked.html)
