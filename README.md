<!-- badges: start -->
[![R-CMD-check](https://github.com/MadsAlbertsen/ampvis2/workflows/R-CMD-check/badge.svg)](https://github.com/MadsAlbertsen/ampvis2/actions)
<!-- badges: end -->

Tools for visualising amplicon data
===================================

ampvis2 is an R-package to conveniently visualise and analyse 16S rRNA amplicon data in different ways.

Installing ampvis2
------------------

First, install [R (3.5.x or later)](https://mirrors.dotsrc.org/cran/) and [RStudio](https://www.rstudio.com/products/rstudio/download/#download). Windows users should also install [RTools](https://mirrors.dotsrc.org/cran/bin/windows/Rtools/). Then open RStudio as administrator (!) and run the commands below to install ampvis2 from the console:

``` r
install.packages("remotes")
remotes::install_github("MadsAlbertsen/ampvis2")
```

Tip: For faster installation you can utilise multicore processors by setting the `Ncpus` argument, fx `remotes::install_github("madsalbertsen/ampvis2", Ncpus = 6)`. Most CPU's today can run 8 processes simultaneously, so setting it to 6 is a good starting point unless you know you have a CPU with more (logical) cores than 8.

Get started
-----------

For a quick guide on how to use ampvis2 go to the [Get Started](https://madsalbertsen.github.io/ampvis2/articles/ampvis2.html) page. Detailed documentation of all ampvis2 functions can be found at the [Functions](https://madsalbertsen.github.io/ampvis2/reference/index.html) page.

RStudio Docker container
------------------

A Docker container based on the [rocker/rstudio](https://hub.docker.com/r/rocker/rstudio) image is also provided with ampvis2 preinstalled. This is ideal for complete reproducibility and portability. All you need to do is [install Docker](https://docs.docker.com/) and then run:

``` bash
docker run -d \
  -e "PASSWORD=supersafepassword" \
  -v "local/folder/to/mount":/home/rstudio \
  -p 8787:8787 \
  ghcr.io/madsalbertsen/ampvis2:main
```

Access RStudio server through a browser at `http://localhost:8787` with username `rstudio`. Ideally use a specific version tag, fx `v2.7.12`, instead of `main` to not just pull the latest image every time.

Blog posts about ampvis2
------------------------

Check out the blog posts at <http://albertsenlab.org/> about selected ampvis2 plotting functions. The posts include details as well as example code:

-   [ampvis2: The bread and butter of our amplicon analyses: amp\_heatmap!](http://albertsenlab.org/ampvis2-heatmap/)
-   [ampvis2: A guide to ordination and how to use amp\_ordinate in R](http://albertsenlab.org/ampvis2-ordination/)
-   [Analysing amplicon data, how easy does it get?](http://albertsenlab.org/shinyampvis/)

Shiny app
-----

An interactive Shiny app with some of the basic functionality of ampvis2 can be found at: <https://kasperskytte.shinyapps.io/shinyampvis>
