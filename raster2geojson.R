library(raster)
library(geojsonio)

slf2014 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/initial_infections_2014_single_count_pm.tif")
slf2015 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/initial_infections_2015_single_count_pm.tif")
slf2016 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/initial_infections_2016_single_count_pm.tif")
slf2017 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/initial_infections_2017_single_count_pm.tif")
slf2017_2 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/initial_infections_2017_single_count_pm2.tif")
slf2018 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/initial_infections_2018_single_count_pm.tif")

rcl <- c(1, Inf, 1, 0, 0.99, NA)
rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
infected <- sum(slf2018[slf2018 >0])
slf2018 <- raster::reclassify(slf2018, rclmat)
slf2018[is.na(slf2018)] <- 0
area <- sum(slf2018[slf2018 > 0]) * raster::xres(slf2018) * raster::yres(slf2018)



slf2014[slf2014 == 0] <- NA
slf2014 <- as.integer(slf2014)
names(slf2014) <- "outputs"
slf2014_p <- raster::rasterToPolygons(slf2014, n = 4, digits = 0, dissolve = T, na.rm = TRUE)
storage.mode(slf2014_p$outputs) <- "integer"
geojsonio::geojson_write(slf2014_p, convert_wgs84 = TRUE, geometry = "polygon", precision = 4, file = "C:/Users/Chris/Desktop/Geojson SLF data/slf_2014.geojson")

slf2015[slf2015 == 0] <- NA
slf2015 <- as.integer(slf2015)
names(slf2015) <- "outputs"
slf2015_p <- raster::rasterToPolygons(slf2015, n = 4, digits = 0, dissolve = T, na.rm = TRUE)
storage.mode(slf2015_p$outputs) <- "integer"
geojsonio::geojson_write(slf2015_p, convert_wgs84 = TRUE, geometry = "polygon", precision = 4, file = "C:/Users/Chris/Desktop/Geojson SLF data/slf_2015.geojson")


slf2016[slf2016 == 0] <- NA
slf2016 <- as.integer(slf2016)
names(slf2016) <- "outputs"
slf2016_p <- raster::rasterToPolygons(slf2016, n = 4, digits = 0, dissolve = T, na.rm = TRUE)
storage.mode(slf2016_p$outputs) <- "integer"
geojsonio::geojson_write(slf2016_p, convert_wgs84 = TRUE, geometry = "polygon", precision = 4, file = "C:/Users/Chris/Desktop/Geojson SLF data/slf_2016.geojson")

slf2017[slf2017 == 0] <- NA
slf2017 <- as.integer(slf2017)
names(slf2017) <- "outputs"
slf2017_p <- raster::rasterToPolygons(slf2017, n = 4, digits = 0, dissolve = T, na.rm = TRUE)
storage.mode(slf2017_p$outputs) <- "integer"
geojsonio::geojson_write(slf2017_p, convert_wgs84 = TRUE, geometry = "polygon", precision = 4, file = "C:/Users/Chris/Desktop/Geojson SLF data/slf_2017.geojson")

slf2018[slf2018 == 0] <- NA
slf2018 <- as.integer(slf2018)
names(slf2018) <- "outputs"
slf2018_p <- raster::rasterToPolygons(slf2018, n = 4, digits = 0, dissolve = T, na.rm = TRUE)
storage.mode(slf2018_p$outputs) <- "integer"
geojsonio::geojson_write(slf2018_p, convert_wgs84 = TRUE, geometry = "polygon", precision = 4, file = "C:/Users/Chris/Desktop/Geojson SLF data/slf_2018.geojson")

## host for SLF case study
tree_of_heaven <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/tree_of_heaven_0.50.tif")
tree_of_heaven[tree_of_heaven == 0] <- NA
tree_of_heaven <- as.integer(tree_of_heaven)
rcl <- c(1, 10, 10, 11, 20, 20, 21, 30, 30, 31, 40, 40, 41, 50, 50, 51, 60, 60, 61, 70, 70, 71, 80, 80, 81, 90, 90, 91, 100, 100)
rcl <- c(1, 20, 20, 21, 40, 40, 41, 60, 60, 61, 80, 80, 81, 100, 100)
rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
tree_of_heaven <- raster::reclassify(tree_of_heaven, rclmat)
names(tree_of_heaven) <- "outputs"
tree_of_heaven_p <- raster::rasterToPolygons(tree_of_heaven, n = 4, digits = 0, dissolve = T, na.rm = TRUE)
storage.mode(tree_of_heaven_p$outputs) <- "integer"
geojsonio::geojson_write(tree_of_heaven_p, convert_wgs84 = TRUE, geometry = "polygon", precision = 3, file = "C:/Users/Chris/Desktop/Geojson SLF data/tree_of_heaven.geojson")

## host for SOD EU1 case study
tanoak <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/LEMMA Data Curry County/lide_300m_median_2018.tif")
tanoak[tanoak == 0] <- NA
tanoak <- as.integer(tanoak)
rcl <- c(1, 10, 10, 11, 20, 20, 21, 30, 30, 31, 40, 40, 41, 50, 50, 51, 60, 60, 61, 70, 70, 71, 80, 80, 81, 90, 90, 91, 100, 100)
rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
tanoak <- raster::reclassify(tanoak, rclmat)
names(tanoak) <- "outputs"
tanoak_p <- raster::rasterToPolygons(tanoak, n = 4, digits = 0, dissolve = T, na.rm = TRUE)
storage.mode(tanoak_p$outputs) <- "integer"
geojsonio::geojson_write(tanoak_p, convert_wgs84 = TRUE, geometry = "polygon", precision = 3, file = "C:/Users/Chris/Desktop/Geojson SLF data/tanoak_eu1.geojson")


na2018 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/LEMMA Data Curry County/100m/cum_inf_2019.tif")
rcl <- c(1, Inf, 1, 0, 0.99, NA)
rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
infected <- sum(na2018[na2018 >0])
na2018 <- raster::reclassify(na2018, rclmat)
na2018[is.na(na2018)] <- 0
area <- sum(na2018[na2018 > 0]) * raster::xres(na2018) * raster::yres(na2018)


eu2018 <- raster("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/LEMMA Data Curry County/100m/cum_inf_2019eu.tif")
rcl <- c(1, Inf, 1, 0, 0.99, NA)
rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
infected <- sum(eu2018[eu2018 >0])
eu2018 <- raster::reclassify(eu2018, rclmat)
eu2018[is.na(eu2018)] <- 0
area <- sum(eu2018[eu2018 > 0]) * raster::xres(eu2018) * raster::yres(eu2018)

