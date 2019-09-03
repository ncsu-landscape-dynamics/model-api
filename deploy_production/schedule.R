# devtools::install_github("ncsu-landscape-dynamics/rpops", ref="feature/updated_dispersal_kernel", force = TRUE)
Sys.setenv("GCS_AUTH_FILE" = "deploy/auth.json")
library(googleAuthR)         ## authentication
library(googleCloudStorageR)  ## google cloud storage
library(readr)                ##
# library(flyio)
# gcs auto authenticated via environment file
# pointed to via sys.env GCS_AUTH_FILE
library(PoPS)
library(httr)
library(geojsonio)
library(protolite)
library(geojson)
library(rgdal)
library(raster)
library(jsonlite)
library(sp)
library(devtools)
library(doParallel)
library(foreach)
library(parallel)

#' Plot out data from the iris dataset
#' @param run_id required to set the proper run id
#' @param case_study_id required to get the proper case study information
#' @get /status
modelapi <- function(run_id, case_study_id) {
  
  options(digits = 6)
  run_id <- as.numeric(run_id)
  json_run <- httr::GET(paste("https://popsmodel.org/api/run/", run_id, "/?format=json", sep = ""))
  run <- httr::content(json_run)
  run$status <- "READING DATA"
  run <- run[c(2:4,6:length(run))]
  httr::PUT(url = paste("https://popsmodel.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  case_study_id <- as.numeric(case_study_id)
  googleCloudStorageR::gcs_load(file = "casestudy2.Rdata", bucket = "testers-pops")
  
  end_time <- run$final_year
  short_distance_scale <- as.numeric(run$distance_scale)
  reproductive_rate <- as.numeric(run$reproductive_rate)
  efficacy = run$efficacy
  
  if (length(run$management_polygons) <= 1) {
    treatments_file <- ""
    treatment_years <- c(0)
    management <- FALSE
  } else if (length(run$management_polygons) > 1) {
    treatments_file <- run$management_polygons
    treatments <- geojsonio::as.json(treatments_file)
    treatments <- geojson::as.geojson(treatments)
    treatments <- geojsonio::geojson_sp(treatments)
    treatments <- spTransform(treatments, CRS = crs(host))
    treatment_map <- raster::rasterize(treatments, host, fun = "last")
    treatment_map[treatment_map > 1] <- 1 ## rasterize gives rasters in each polygon the polygon id value need to set those to 1
    treatment_map[is.na(treatment_map)] <- 0
    treatment_map <- treatment_map * (efficacy / 100)
    treatment_maps <- list(raster::as.matrix(treatment_map))
    treatment_years <- c(start_time)
    management <- TRUE
  }
  
  if (run$weather == "GOOD") {
    temperature <- high_temperature
  } else if (run$weather == "BAD") {
    temperature <- low_temperature
  } else {
    temperature <- temperature
  }
  
  core_count <- 10
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  run$status <- "RUNNING MODEL"
  httr::PUT(url = paste("https://popsmodel.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  years <- seq(start_time, end_time, 1)
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)

  infected_stack <- foreach::foreach(i = 1:10, .combine = c, .packages = c("raster", "PoPS"), .export = ls(globalenv())) %dopar% {
    random_seed <- round(stats::runif(1, 1, 1000000))
    data <- PoPS::pops_model(random_seed = random_seed,
                             lethal_temperature = lethal_temperature, use_lethal_temperature = use_lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                             reproductive_rate = reproductive_rate,
                             weather = weather, mortality_on = mortality_on,
                             short_distance_scale = short_distance_scale, infected = infected,
                             susceptible = susceptible, mortality_tracker = mortality_tracker, mortality = mortality,
                             total_plants = total_plants,
                             treatment_maps = treatment_maps, treatment_years = treatment_years,
                             temperature = temperature,
                             weather_coefficient = weather_coefficient,
                             ew_res = ew_res, ns_res = ns_res,
                             time_step = time_step, mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                             season_month_start = season_month_start, season_month_end = season_month_end,
                             start_time = start_time, end_time = end_time,
                             dispersal_kern = dispersal_kern, percent_short_distance_dispersal = percent_short_distance_dispersal,
                             long_distance_scale = long_distance_scale, treatment_method = "ratio",
                             wind_dir = wind_dir, kappa = kappa)
    comp_years <- raster::stack(lapply(1:length(data$infected_before_treatment), function(i) host))
    number_infected<- data.frame(t(years))
    for (q in 1:raster::nlayers(comp_years)) {
      comp_years[[q]] <- data$infected[[q]]
    }
    number_infected <- cellStats(comp_years, 'sum')
    comp_years <- raster::reclassify(comp_years, rclmat)
    comp_years[is.na(comp_years)] <- 0
    infected_stack <- comp_years
    data <- list(infected_stack, number_infected)
  }
  stopCluster(cl)
  prediction <- infected_stack[[1]]
  prediction[prediction > 0] <- 0
  infected_area <- data.frame(t(years))
  ii <- 0
  number_infected <- data.frame(t(years))
  
  for (i in seq(1,length(infected_stack), 2)) {
    prediction <- prediction + infected_stack[[i]]
    ii <- ii + 1
    qq <- i + 1
    number_infected[ii,] <- infected_stack[[qq]]
    infected_area[ii,] <- cellStats(infected_stack[[i]], 'sum') * xres(host) * yres(host) # converts to meters squared
  }
  
  probability <- (prediction/(length(infected_stack)/2)) * 100
  infected_areas <- round(sapply(infected_area, function(x) c( "Mean"= mean(x,na.rm=TRUE),
                                                               "Stand dev" = sd(x)
  )), digits = 0)
  number_infecteds <- round(sapply(number_infected, function(x) c( "Mean"= mean(x,na.rm=TRUE),
                                                                   "Stand dev" = sd(x)
  )), digits = 0)
  
  run$status <- "WRITING DATA"
  httr::PUT(url = paste("https://popsmodel.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  
  core_count <- 10
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  stuff <- foreach::foreach(q = 1:length(years), .packages =c("raster", "geojsonio", "httr"), .export = ls(globalenv())) %dopar% {
    number_infected <- number_infecteds[1, q]
    area_infected <- infected_areas[1, q]
    year <- years[q]
    spread_map <- probability[[q]]
    spread_map <- as.integer(spread_map)
    spread_map[spread_map <= 0] <- NA
    names(spread_map) <- "outputs"
    spread_map <- projectRaster(spread_map, crs = CRS("+proj=longlat +datum=WGS84"), method = "ngb")
    spread_map <- raster::rasterToPolygons(spread_map, n = 4, digits = 4, dissolve = T, na.rm = TRUE)
    storage.mode(spread_map$outputs) <- "integer"
    spread_map <- geojsonio::geojson_list(spread_map, convert_wgs84 = TRUE, geometry = "polygon")
    class(spread_map) <- "list"
    outs <- list()
    outs$number_infected <- number_infected
    outs$infected_area <- area_infected
    outs$years <- year
    outs$spread_map <- spread_map
    outs$run <- run_id
    post_code <- httr::POST(url = "https://popsmodel.org/api/output/", body = outs, encode = "json")
    if (post_code$status_code == 201) {
      run$status <- "SUCCESS"
    } else {
      run$status <- "FAILED"
    }
    stuff <- run$status
  }
  run$status <- stuff[[1]]
  
  stopCluster(cl)
  if (run$status == "SUCCESS") {
    httr::PUT(url = paste("https://popsmodel.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  } else {
    httr::PUT(url = paste("https://popsmodel.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  }
  
  status <- run$status
}