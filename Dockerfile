FROM rocker/verse:4.1.1
Label maintainer="Chris Jones <cmjone25@ncsu.edu>"

RUN apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    gdal-bin \
    lbzip2 \
    libfftw3-dev \
    libgdal-dev \
    libgeos-dev \
    libgsl0-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libhdf4-alt-dev \
    libhdf5-dev \
    libjq-dev \
    libpq-dev \
    libproj-dev \
    libprotobuf-dev \
    libnetcdf-dev \
    libsqlite3-dev \
    libssl-dev \
    libudunits2-dev \
    lsb-release \
    netcdf-bin \
    postgis \
    protobuf-compiler \
    sqlite3 \
    tk-dev \
    unixodbc-dev \
    git-core \
    libcurl4-gnutls-dev \
    curl \
    libsodium-dev \
    libxml2-dev


## Instal R packages
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv')"

WORKDIR /pops-model_api
COPY renv.lock renv.lock
COPY env env
COPY main.R main.R
COPY schedule.R schedule.R
ENV RENV_PATHS_LIBRARY renv/library

RUN R -e "renv::restore()"

## open port 8079 to traffic
EXPOSE 8079

# when the container starts, start the main.R script
ENTRYPOINT ["Rscript", "main.R"]
