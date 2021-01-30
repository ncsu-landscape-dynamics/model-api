FROM trestletech/plumber

LABEL maintainer="Chris"
RUN echo 'deb http://deb.debian.org/debian bullseye main' > /etc/apt/sources.list
RUN apt-get update \
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
	libxml2-dev \
	make \
	protobuf-compiler
RUN ["install2.r", "askpass", "assertthat", "aws.s3", "aws.signature", "backports", "base64enc", "callr", "class", "classInt", "cli", "codetools", "crayon", "curl", "DBI", "desc", "devtools", "digest", "doParallel", "e1071", "flyio", "foreach", "foreign", "fs", "geojson", "glue", "googleAuthR", "googleCloudStorageR", "hms", "httr", "iterators", "jqr", "jsonlite", "KernSmooth", "landscapemetrics", "lattice", "lazyeval", "magrittr", "maptools", "memoise", "openssl", "pillar", "pkgbuild", "pkgconfig", "pkgload", "prettyunits", "processx", "protolite", "ps", "R6", "raster", "Rcpp", "readr", "remotes", "rgdal", "rgeos", "rlang", "rprojroot", "rstudioapi", "sessioninfo", "sf", "sp", "stringi", "stringr", "testthat", "tibble", "units", "usethis", "V8", "vctrs", "withr", "xml2", "yaml", "zeallot", "zip"]
RUN ["installGithub.r", "ncsu-landscape-dynamics/rpops@b12e865c490db7cf955767ed9ced1a1f7ba6b2d8"]
WORKDIR /payload/
COPY ["./", "./"]
EXPOSE 8080
ENTRYPOINT ["R", "-e", "pr <- plumber::plumb(commandArgs()[4]); pr$run(host='0.0.0.0', port=8080)"]
CMD ["schedule.R"]
