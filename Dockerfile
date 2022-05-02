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


#install devtools and ampvis2
COPY . /opt/ampvis2
RUN Rscript -e 'install.packages("devtools", Ncpus = parallel::detectCores())' \
  && Rscript -e 'devtools::install("/opt/ampvis2/", reload = FALSE, dependencies = TRUE, Ncpus = parallel::detectCores())'

#enable users to install R packages
RUN chown 1000:1000 -R /usr/local/lib/R/site-library /usr/local/lib/R/library

# enable multithreaded make when installing packages by default
RUN mkdir -p /home/rstudio/.R \
  && echo "MAKEFLAGS = -j" > /home/rstudio/.R/Makevars

#silence RStudio warnings about not being able to write dictionary stuff to /root
VOLUME /root

USER rstudio
WORKDIR /home/rstudio
