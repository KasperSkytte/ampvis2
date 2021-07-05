FROM rocker/rstudio:4.1.0

ARG MAKEFLAGS="-j 10 "

COPY renv.lock .

#install nice-to-have system dependencies for R, and netstat to scan ports
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qqy && \
  apt-get install -y --no-install-recommends --no-install-suggests \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libjpeg-dev \
    net-tools

#set default renv cache path in container,
#install renv, install all packages in the lock file to /usr/local/lib/R/site-library/ in container
RUN echo "RENV_PATHS_CACHE=/usr/local/lib/R/renv-cache/" >> /usr/local/lib/R/etc/Renviron.site && \
  R -e "install.packages('renv', repos='http://cran.rstudio.com/', Ncpus = 10)" && \
  R -e "renv::consent(provided = TRUE)" && \
  R -e "renv::restore( \
    library = '/usr/local/lib/R/site-library/', \
    rebuild = TRUE, \
     clean = TRUE, \
     lockfile = 'renv.lock', \
     prompt = FALSE)"
RUN R -e "renv::install('madsalbertsen/ampvis2@2.7.7')"

#silence RStudio warnings about not being able to write dictionary stuff to /root
VOLUME /root
