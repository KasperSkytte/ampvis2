# R version: 4, 4.1, 4.0
ARG VARIANT="4"
FROM rocker/r-ver:${VARIANT}

# Use the [Option] comment to specify true/false arguments that should appear in VS Code UX
#
# [Option] Install zsh
ARG INSTALL_ZSH="true"
# [Option] Upgrade OS packages to their latest versions
ARG UPGRADE_PACKAGES="false"

# ampvis2 release to install
ARG ampvis2_rel="2.7.11"

# Install needed packages and setup non-root user. Use a separate RUN statement to add your own dependencies.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID
COPY library-scripts/*.sh /tmp/library-scripts/
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && /bin/bash /tmp/library-scripts/common-debian.sh "${INSTALL_ZSH}" "${USERNAME}" "${USER_UID}" "${USER_GID}" "${UPGRADE_PACKAGES}" "true" "true" \
    && usermod -a -G staff ${USERNAME} \
    && apt-get -y install \
        python3-pip \
        libgit2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libxt-dev \
        libcairo2-dev \
        libjpeg-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/* /tmp/library-scripts \
    && python3 -m pip --no-cache-dir install radian \
    && echo "options(repos = c(CRAN = 'https://mirrors.dotsrc.org/cran/'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site \
    && install2.r --error --skipinstalled --ncpus -1 \
        devtools \
        languageserver \
        httpgd \
        renv \
    && R -e "renv::consent(provided = TRUE)" \
    && R -e "renv::restore( \
      library = '/usr/local/lib/R/site-library/', \
      rebuild = TRUE, \
      clean = TRUE, \
      lockfile = 'https://raw.githubusercontent.com/MadsAlbertsen/ampvis2/${ampvis2_rel}/renv.lock', \
      prompt = FALSE)" \
    && R -e "renv::install('madsalbertsen/ampvis2@${ampvis2_rel}')" \
    && rm -rf /tmp/downloaded_packages

# VSCode R Debugger dependency, install the latest release version from GitHub.
RUN Rscript -e 'remotes::install_github("ManuelHentschel/vscDebugger@*release", dependencies = FALSE)'

# R Session watcher settings.
# See more details: https://github.com/REditorSupport/vscode-R/wiki/R-Session-watcher
RUN echo 'source(file.path(Sys.getenv("HOME"), ".vscode-R", "init.R"))' >> ${R_HOME}/etc/Rprofile.site