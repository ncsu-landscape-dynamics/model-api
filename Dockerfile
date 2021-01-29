FROM rocker/r-ver:4.0.2
LABEL maintainer="Chris"
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y gdal-bin \
	git-core \
	libgdal-dev \
	libgeos-dev \
	libgeos++-dev \
	libprotobuf-dev \
	libprotoc-dev \
	protobuf-compiler
RUN ["install2.r", "aws.s3", "devtools", "doParallel", "foreach", "geojson", "geojsonio", "httr", "iterators", "jsonlite", "protolite", "raster", "remotes", "rgdal", "sf", "sp", "terra", "usethis"]
RUN ["installGithub.r", "ncsu-landscape-dynamics/rpops@b12e865c490db7cf955767ed9ced1a1f7ba6b2d8"]
WORKDIR /payload/
COPY ["deploy_model/", "deploy_model/"]
CMD ["R"]
