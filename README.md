<!-- badges: start -->

[![R-CMD-check](https://github.com/kasperskytte/ampvis2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kasperskytte/ampvis2/actions/workflows/R-CMD-check.yaml)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-ampvis2/README.html)

<!-- badges: end -->

Tools for visualising amplicon data
===================================

ampvis2 is an R-package to conveniently visualise and analyse 16S rRNA amplicon data in different ways.

Installing ampvis2
------------------

First, install [R (3.5.x or later)](https://mirrors.dotsrc.org/cran/) and [RStudio](https://www.rstudio.com/products/rstudio/download/#download). Windows users should also install [RTools](https://mirrors.dotsrc.org/cran/bin/windows/Rtools/). Then open RStudio as administrator (!) and run the commands below to install ampvis2 from the console:

``` r
install.packages("remotes")
remotes::install_github("kasperskytte/ampvis2", Ncpus = 4)
```

You can also [install ampvis2 into a conda environment](http://bioconda.github.io/recipes/r-ampvis2/README.html) from the [bioconda channel](https://bioconda.github.io/index.html), or use the Docker container:

```
docker pull quay.io/biocontainers/r-ampvis2:<tag>
```

Get started
-----------

For a quick guide on how to use ampvis2 go to the [Get Started](https://kasperskytte.github.io/ampvis2/articles/ampvis2.html) page. Detailed documentation of all ampvis2 functions can be found at the [Functions](https://kasperskytte.github.io/ampvis2/reference/index.html) page.

Shiny app
-----

An interactive Shiny app with some of the basic functionality of ampvis2 can be found at: <https://kasperskytte.shinyapps.io/shinyampvis>
