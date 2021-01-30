FROM rocker/verse:4.0.2
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
RUN install2.r --error --skipinstalled \
    aws.s3 \
    RColorBrewer \
    RandomFields \
    RNetCDF \
    classInt \
    deldir \
    gstat \
    hdf5r \
    lidR \
    mapdata \
    maptools \
    mapview \
    ncdf4 \
    proj4 \
    raster \
    rgdal \
    rgeos \
    rlas \
    sf \
    sp \
    spacetime \
    spatstat \
    spdep \
    geoR \
    geosphere \
    devtools \
    doParallel \
    foreach \
    geojsonio \
    httr \
    iterators \
    jsonlite \
    plumber \
    protolite \
    raster \
    remotes \
    rgdal \
    sf \
    sp \
    terra \
    usethis
    
RUN R -e "install.packages(c('terra'), repos = 'http://cran.us.r-project.org')"

## RUN R -e "install.packages(c('RColorBrewer', 'RandomFields', 'RNetCDF', 'classInt', 'deldir', 'gstat', 'hdf5r', 'lidR', 'mapdata', 'maptools', 'mapview', 'ncdf4, 'proj4', 'raster', 'rgdal', 'rgeos', 'rlas', 'sf', 'sp', 'spacetime', 'spatstat', 'spdep', 'geoR', 'geosphere', 'aws.s3', 'devtools', 'doParallel', 'foreach', 'geojson', 'geojsonio', 'httr', 'iterators', 'jsonlite', 'plumber', 'protolite', 'raster', 'remotes', 'rgdal', 'sf', 'sp', 'terra', 'usethis'))"

RUN ["installGithub.r", "ncsu-landscape-dynamics/rpops@b12e865c490db7cf955767ed9ced1a1f7ba6b2d8"]
COPY / /

## open port 8079 to traffic
EXPOSE 8079

# when the container starts, start the main.R script
ENTRYPOINT ["Rscript", "main.R"]
