FROM rocker/r-ver:4.0.2
Maintainer "Chris" cmjone25@ncsu.edu

RUN apt-get -y update \
  && apt-get install -y --no-install-recommends \
    git-core \
    libcurl4-gnutls-dev \
    curl \
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
    liblwgeom-dev \
    libpq-dev \
    libproj-dev \
    libprotobuf-dev \
    libnetcdf-dev \
    libsqlite3-dev \
    libssl-dev \
    libudunits2-dev \
    netcdf-bin \
    postgis \
    protobuf-compiler \
    sqlite3 \
    tk-dev \
    libsodium-dev \
    libxml2-dev \
    unixodbc-dev
    
RUN install2.r --error \
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
    aws.s3 \
    devtools \
    doParallel \
    foreach \
    geojson \
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

RUN ["installGithub.r", "ncsu-landscape-dynamics/rpops@b12e865c490db7cf955767ed9ced1a1f7ba6b2d8"]
WORKDIR /payload/
COPY ["./", "./"]
EXPOSE 8080
ENTRYPOINT ["R", "-e", "pr <- plumber::plumb(commandArgs()[4]); pr$run(host='0.0.0.0', port=8080)"]
CMD ["schedule.R"]
