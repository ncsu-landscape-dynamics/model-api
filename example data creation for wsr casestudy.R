library(terra)
library(sp)
library(raster)
library(sf)
library(stars)
folderfun::setff("Out", "C:/Users/Chris/Desktop/wsr_online_cs/")

res <- 1.524 # resolution of a pixel or cell
ns_ext <- 5 # north south extent of the plots
we_ext <- 28 # east west extent of the plots

total_hosts <- terra::rast(nrow = ns_ext, ncol = we_ext,
                           xmin = -13677191.35, xmax = -13677191.35 + we_ext * res, ymin = 5554304.69, ymax = 5554304.69 + ns_ext * res)
crs(total_hosts) <- "epsg:3857"
values(total_hosts) <- rep(1000, ns_ext * we_ext)
plot(total_hosts)
th <- terra::as.matrix(total_hosts, wide = TRUE)



terra::writeRaster(total_hosts, ffOut("total_hosts.tif"))

hosts <- terra::rast(nrow = ns_ext, ncol = we_ext,
                           xmin = -13677191.35, xmax = -13677191.35 + we_ext * res, ymin = 5554304.69, ymax = 5554304.69 + ns_ext * res)
crs(hosts) <- "epsg:3857"
values(hosts) <- rep(1000, ns_ext * we_ext)
plot(hosts)
names(hosts) <- "outputs"
hosts_p <- terra::as.polygons(hosts, digits = 4, dissolve = TRUE,
                                   values = TRUE, na.rm = TRUE)
storage.mode(hosts_p$outputs) <- "integer"

setAs("SpatVector", "sf",
      function(from) {
        sf::st_as_sf(as.data.frame(from, geom=TRUE),
                     wkt="geometry", crs=crs(from))
      }
)

st_as_sf.SpatVector <- function(x, ...) {
  sf::st_as_sf(as.data.frame(x, geom=TRUE), wkt="geometry", crs=crs(x))
}

hosts_p <- terra::project(hosts_p, "epsg:4326")
hosts_p <- as(hosts_p, "sf")
hosts_p <- sf::st_as_sf(hosts_p)

hosts_p <-
  geojsonio::geojson_list(hosts_p, precision = 4, geometry = "polygon")
# class(hosts_p) <- "list"

geojsonio::geojson_write(hosts_p, file = "host.geojson")






terra::writeRaster(hosts, ffOut("hosts.tif"))

infected <- terra::rast(nrow = ns_ext, ncol = we_ext,
                           xmin = -13677191.35, xmax = -13677191.35 + we_ext * res, ymin = 5554304.69, ymax = 5554304.69 + ns_ext * res)
crs(infected) <- "epsg:3857"
values(infected) <- rep(0, ns_ext * we_ext)
infected[3, 4] <- 11
plot(infected)

names(infected) <- "outputs"
infected_p <- terra::as.polygons(infected, digits = 4, dissolve = TRUE,
                              values = TRUE, na.rm = TRUE)
storage.mode(infected_p$outputs) <- "integer"

setAs("SpatVector", "sf",
      function(from) {
        sf::st_as_sf(as.data.frame(from, geom=TRUE),
                     wkt="geometry", crs=crs(from))
      }
)

st_as_sf.SpatVector <- function(x, ...) {
  sf::st_as_sf(as.data.frame(x, geom=TRUE), wkt="geometry", crs=crs(x))
}

infected_p <- terra::project(infected_p, "epsg:4326")
infected_p <- as(infected_p, "sf")
infected_p <- sf::st_as_sf(infected_p)

infected_p <-
  geojsonio::geojson_list(infected_p, precision = 4, geometry = "polygon")
# class(infected_p) <- "list"

geojsonio::geojson_write(infected_p, file = "infected.geojson")


th <- terra::as.matrix(infected, wide = TRUE)
th <- matrix(infected, ncol = terra::ncol(infected), nrow = terra::nrow(infected))

terra::writeRaster(infected, ffOut("infected.tif"))

setAs("SpatVector", "sf",
      function(from) {
        sf::st_as_sf(as.data.frame(from, geom=TRUE), wkt="geometry", crs=crs(from))
      }
)

st_as_sf.SpatVector <- function(x, ...) {
  sf::st_as_sf(as.data.frame(x, geom=TRUE), wkt="geometry", crs=crs(x))
}

total_hosts_p <- terra::as.polygons(total_hosts, dissolve = TRUE)
total_hosts_p <- terra::project(total_hosts_p, "epsg:4326")
total_hosts_p <- as(total_hosts_p, "sf")
total_hosts_ps <- sf::st_as_sf(total_hosts_p)
st_write(total_hosts_ps, dsn = ffOut("total_hosts.geojson"), driver = "GeoJSON")

hosts_p <- terra::as.polygons(hosts, dissolve = TRUE)
hosts_p <- terra::project(hosts_p, "epsg:4326")
hosts_p <- as(hosts_p, "sf")
hosts_ps <- sf::st_as_sf(hosts_p)
st_write(hosts_ps, dsn = ffOut("hosts.geojson"), driver = "GeoJSON")

infected_p <- terra::as.polygons(infected, dissolve = TRUE)
infected_p <- terra::project(infected_p, "epsg:4326")
infected_p <- as(infected_p, "sf")
infected_ps <- sf::st_as_sf(infected_p)
st_write(infected_ps, dsn = ffOut("infected.geojson"), driver = "GeoJSON")



## change airport data
folderfun::setff("In",
                 "C:/Users/Chris/Desktop/High Risk Locations for PoPS Dashboard/High Risk Locations for PoPS Dashboard")

folderfun::setff("Out", "C:/Users/Chris/Desktop/Risk Locations/")

airports <- st_read(ffIn("airports_conus.shp"))
airports <- st_transform(airports, crs = 4326)
airports_g <-
  geojsonio::geojson_list(airports, precision = 4, geometry = "polygon")
geojsonio::geojson_write(airports_g, file = ffOut("airports.geojson"))


intermodal_frieght <- st_read(ffIn("intermodal_freight_conus.shp"))
intermodal_frieght <- st_transform(intermodal_frieght, crs = 4326)
intermodal_frieght_g <-
  geojsonio::geojson_list(intermodal_frieght, precision = 4, geometry = "polygon")
geojsonio::geojson_write(intermodal_frieght_g, file = ffOut("intermodal_frieght.geojson"))


intermodal_passenger <- st_read(ffIn("intermodal_passenger_conus.shp"))
intermodal_passenger <- st_transform(intermodal_passenger, crs = 4326)
intermodal_passenger_g <-
  geojsonio::geojson_list(intermodal_passenger, precision = 4, geometry = "polygon")
geojsonio::geojson_write(intermodal_passenger_g, file = ffOut("intermodal_passenger.geojson"))


seaports <- st_read(ffIn("seaports_conus.shp"))
seaports <- st_transform(seaports, crs = 4326)
seaports_g <-
  geojsonio::geojson_list(seaports, precision = 4, geometry = "polygon")
geojsonio::geojson_write(seaports_g, file = ffOut("seaports.geojson"))


trucking <- st_read(ffIn("truck_parking_conus.shp"))
trucking <- st_transform(trucking, crs = 4326)
trucking_g <-
  geojsonio::geojson_list(trucking, precision = 4, geometry = "polygon")
geojsonio::geojson_write(trucking_g, file = ffOut("trucking.geojson"))


folderfun::setff("In", "H:/Shared drives/Data/Vector/USA/")
rails <- st_read(ffIn("railroads_US.gpkg"))
st_write(rails, ffIn("railroads_US.shp"))
