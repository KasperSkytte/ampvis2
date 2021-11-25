FROM rocker/rstudio

#multithreaded make
ENV MAKEFLAGS="-j"

#nice-to-have system dependencies for R
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qqy && \
  apt-get install -y --no-install-recommends --no-install-suggests \
    libcairo2-dev \
    libxt-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libicu-dev \
    make \
    git \
    pandoc \
    libxml2-dev \
    libgit2-dev

RUN Rscript -e 'install.packages("devtools", Ncpus = parallel::detectCores())'
COPY . /opt/ampvis2/

RUN Rscript -e 'devtools::install("/opt/ampvis2/", reload = FALSE, dependencies = TRUE, Ncpus = parallel::detectCores())'

#set default renv cache path in container
#change CRAN mirror from https://packagemanager.rstudio.com to AAU mirror
#install renv, and install all packages in the lock file to /usr/local/lib/R/site-library/ in container
#RUN echo "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" >> /usr/local/lib/R/etc/Renviron.site && \\


#enable users to install R packages
#RUN chown 1000:1000 -R /usr/local/lib/R/site-library /usr/local/lib/R/library

#silence RStudio warnings about not being able to write dictionary stuff to /root
VOLUME /root

WORKDIR /home/rstudio