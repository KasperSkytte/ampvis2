#!/usr/bin/env bash
set -eu

#set variables
r_ver="4.1.0"
ampvis2_rel="2.7.6"
image_name="kasperskytte/rstudio_r${r_ver}_ampvis2:${ampvis2_rel}"
password="supersafepassword"
port="8787" #will generate a random one if unavailable
RENV_PATHS_CACHE_HOST="${HOME}/.local/share/renv/cache" #path to renv cache on host, renv default is ~/.local/share/renv/cache
RENV_PATHS_CACHE_CONTAINER="/usr/local/lib/R/renv-cache/" #path to renv cache within the container (dont have to change)
num_threads=$(($(nproc) - 2)) #all cores except 2
#done setting variables

wget -q https://raw.githubusercontent.com/MadsAlbertsen/ampvis2/2.7.5/renv.lock -O renv.lock

cat << Dockerfile > Dockerfile
FROM rocker/rstudio:${r_ver}

ARG MAKEFLAGS="-j ${num_threads} "

COPY renv.lock .

#install nice-to-have system dependencies for R, and netstat to scan ports
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qqy && \\
  apt-get install -y --no-install-recommends --no-install-suggests \\
    libxml2-dev \\
    libcairo2-dev \\
    libxt-dev \\
    libjpeg-dev \\
    net-tools

#set default renv cache path in container,
#install renv, install all packages in the lock file to /usr/local/lib/R/site-library/ in container
RUN echo "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" >> /usr/local/lib/R/etc/Renviron.site && \\
  R -e "install.packages('renv', repos='http://cran.rstudio.com/', Ncpus = ${num_threads})" && \\
  R -e "renv::consent(provided = TRUE)" && \\
  R -e "renv::restore( \\
    library = '/usr/local/lib/R/site-library/', \\
    rebuild = TRUE, \\
     clean = TRUE, \\
     lockfile = 'renv.lock', \\
     prompt = FALSE)"
RUN R -e "renv::install('madsalbertsen/ampvis2@${ampvis2_rel}')"
Dockerfile

#optionally pull the image instead of building
if [ "$1" == "pull" ]
then
  docker pull "${image_name}"
else
  #build the container image
  docker build -t "${image_name}" .
fi

checkPort() {
  randomPort() {
    echo $(( ( RANDOM % 60000 )  + 1025 ))
  }
  local port
  local check_port
  check_port=${1:-"$(randomPort)"}

  while [ -n "$check_port" ]
  do
    port="$check_port"
    check_port=$(docker run --rm --net=host ${image_name} netstat -atn | grep "$port" || :)
    if [ -n "$check_port" ]
    then
      check_port="$(randomPort)"
    fi 
  done
  echo "$port"
}

port=$(checkPort "$port")

mkdir -p ${RENV_PATHS_CACHE_HOST}

#launch the container with the host cache mounted in the container
docker run -d --rm \
  -e "PASSWORD=${password}" \
  -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
  -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
  -v "$HOME":/home/rstudio \
  -p "$port":8787 \
  ${image_name}

echo
echo "Launch RStudio through a browser at one of these adresses:"
echo "http://127.0.0.1:${port} (this machine only)"
for IP in $(hostname -I)
do
  echo "http://${IP}:${port}"
done
echo
echo "Username: rstudio"
echo "Password: ${password}"
