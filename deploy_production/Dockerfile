FROM trestletech/plumber
LABEL maintainer="Chris"
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y gdal-bin \
	git-core \
	libcurl4-openssl-dev \
	libgdal-dev \
	libgeos-dev \
	libgeos++-dev \
	libjq-dev \
	libproj-dev \
	libprotobuf-dev \
	libprotoc-dev \
	libssl-dev \
	libudunits2-dev \
	libv8-dev \
	make \
	protobuf-compiler
RUN ["install2.r", "askpass", "assertthat", "backports", "callr", "class", "classInt", "cli", "codetools", "crayon", "curl", "DBI", "desc", "devtools", "digest", "doParallel", "e1071", "foreach", "foreign", "fs", "geojson", "geojsonio", "glue", "googleAuthR", "googleCloudStorageR", "hms", "httr", "iterators", "jqr", "jsonlite", "KernSmooth", "landscapemetrics", "lattice", "lazyeval", "magrittr", "maptools", "memoise", "openssl", "pillar", "pkgbuild", "pkgconfig", "pkgload", "prettyunits", "processx", "protolite", "ps", "R6", "raster", "Rcpp", "readr", "remotes", "rgdal", "rgeos", "rlang", "rprojroot", "rstudioapi", "sessioninfo", "sf", "sp", "testthat", "tibble", "units", "usethis", "V8", "withr", "yaml", "zip"]
RUN ["installGithub.r", "ncsu-landscape-dynamics/rpops@0104394a78d4d98075e6c153fccf2d7367a397fa"]
WORKDIR /payload/
COPY ["./", "./"]
EXPOSE 8080
ENTRYPOINT ["R", "-e", "pr <- plumber::plumb(commandArgs()[4]); pr$run(host='0.0.0.0', port=8080)"]
CMD ["schedule.R"]
