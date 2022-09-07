FROM rocker/rstudio

#multithreaded make for container build time
ENV MAKEFLAGS="-j"

#default password for RStudio Server
ENV PASSWORD=supersafepassword

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

# enable multithreaded make when installing packages by default
RUN mkdir -p /home/rstudio/.R \
  && echo "MAKEFLAGS = -j" > /home/rstudio/.R/Makevars

#install devtools and ampvis2
COPY . /opt/ampvis2
RUN Rscript -e 'install.packages("devtools", Ncpus = parallel::detectCores())' \
  && Rscript -e 'devtools::install("/opt/ampvis2/", reload = FALSE, dependencies = TRUE, Ncpus = parallel::detectCores())'
