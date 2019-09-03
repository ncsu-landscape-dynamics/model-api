Sys.setenv("GCS_AUTH_FILE" = "deploy_production/auth.json")
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
library(flyio)

case_study_id = "2"
case_study_id <- as.numeric(case_study_id)
json_case_study <- httr::GET(paste("https://popsmodel.org/api/case_study/", case_study_id ,"/?format=json", sep = ""))

case_study <- content(json_case_study)


temp <- case_study$weather$temp_on
precip <- case_study$weather$precipitation_on
use_lethal_temperature <- case_study$weather$lethal_temp_on
season_month_start <- case_study$weather$seasonality$first_month 
season_month_end <- case_study$weather$seasonality$last_month 
time_step <- case_study$time_step
start_time <- case_study$end_year + 1
end_time <- case_study$future_years
dispersal_kern <- case_study$pest_set[[1]]$dispersal_type
if (dispersal_kern == "CAUCHY") {
  dispersal_kern <- "cauchy"
}

lethal_temperature <- as.numeric(case_study$weather$lethaltemperature$value)
lethal_temperature <- -14
lethal_temperature_month <- case_study$weather$lethaltemperature$month
if (case_study$weather$wind_on == FALSE) {
  wind_dir <- "NONE"
  kappa <- 0
} else {
  wind_dir <- case_study$weather$wind
  kappa <- case_study$weather$wind
}

percent_short_distance_dispersal <- 1.0
long_distance_scale <- 0.0

treatments_file <- ""
treatment_years <- c(0)
management <- FALSE


mortality_on <- case_study$host_set[[1]]$mortality_on
if (case_study$host_set[[1]]$mortality_on == FALSE) {
  mortality_rate <- 0
  mortality_time_lag <- 0
} else {
  mortality_rate <- case_study$host_set[[1]]$mortality # need to set with values from api
  mortality_time_lag <- case_study$host_set[[1]]$mortality # need to set with values from api
}

if (time_step == "week") {
  number_of_time_steps <- (end_time-start_time+1)*52
} else if (time_step == "month") {
  number_of_time_steps <- (end_time-start_time+1)*12
} else if (time_step == "day") {
  number_of_time_steps <- (end_time-start_time+1)*365
}

number_of_years <- end_time-start_time+1

flyio_set_datasource("gcs")
flyio_set_bucket("testers-pops")

infected <- import_raster(file = "initial_infections_2018_single_count_pm_prop.tif", FUN = raster)
infected[is.na(infected)] <- 0
host <- import_raster(file = "tree_of_heaven_0.50.tif", FUN = raster)
host[is.na(host)] <- 0
total_plants <- import_raster(file = "total_hosts.tif", FUN = raster)
total_plants[is.na(total_plants)] <- 0

if (!(raster::extent(infected) == raster::extent(host) && raster::extent(infected) == raster::extent(total_plants))) {
  return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
}

if (!(raster::xres(infected) == raster::xres(host) && raster::xres(infected) == raster::xres(total_plants) && raster::yres(infected) == raster::yres(host) && raster::yres(infected) == raster::yres(total_plants))) {
  return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
}

if (!(raster::compareCRS(host,infected) && raster::compareCRS(host, total_plants))) {
  return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
}

susceptible <- host - infected
susceptible[is.na(susceptible)] <- 0
susceptible[susceptible < 0] <- 0

# if (use_lethal_temperature == TRUE  && !file.exists(temperature_file)) {
#   return("Temperature file does not exist")
# }

# if (use_lethal_temperature == TRUE  && !(raster::extension(temperature_file) %in% c(".grd", ".tif", ".img"))) {
#   return("Temperature file is not one of '.grd', '.tif', '.img'")
# }

if (use_lethal_temperature == TRUE) {
  avg_temperature_stack <- import_raster(file = "avg_spread_crit_temp_slf_2018_2022_pm.tif", FUN = stack)
  avg_temperature_stack[is.na(avg_temperature_stack)] <- 0
  low_temperature_stack <- import_raster(file = "low_spread_crit_temp_slf_2018_2022_pm.tif", FUN = stack)
  low_temperature_stack[is.na(low_temperature_stack)] <- 0
  high_temperature_stack <- import_raster(file = "high_spread_crit_temp_slf_2018_2022_pm.tif", FUN = stack)
  high_temperature_stack[is.na(high_temperature_stack)] <- 0
  
  if (!(raster::extent(infected) == raster::extent(avg_temperature_stack))) {
    return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  }
  
  if (!(raster::xres(infected) == raster::xres(avg_temperature_stack) && raster::yres(infected) == raster::yres(avg_temperature_stack))) {
    return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  }
  
  if (!(raster::compareCRS(infected, avg_temperature_stack))) {
    return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  }
  
  avg_temperature <- list(raster::as.matrix(avg_temperature_stack[[1]]))
  low_temperature <- list(raster::as.matrix(low_temperature_stack[[1]]))
  high_temperature <- list(raster::as.matrix(high_temperature_stack[[1]]))
  for(i in 2:number_of_years) {
    avg_temperature[[i]] <- raster::as.matrix(avg_temperature_stack[[i]])
    low_temperature[[i]] <- raster::as.matrix(low_temperature_stack[[i]])
    high_temperature[[i]] <- raster::as.matrix(high_temperature_stack[[i]])
  }
} else {
  temperature <- host
  raster::values(temperature) <- 1
  temperature <- list(raster::as.matrix(temperature))
}

if (temp == TRUE  && !file.exists(temperature_coefficient_file)) {
  return("Temperature coefficient file does not exist")
}

if (temp == TRUE  && !(raster::extension(temperature_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
  return("Temperature coefficient file is not one of '.grd', '.tif', '.img'")
}

if (precip == TRUE  && !file.exists(precipitation_coefficient_file)) {
  return("Precipitation coefficient file does not exist")
}

if (precip == TRUE  && !(raster::extension(precipitation_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
  return("Precipitation coefficient file is not one of '.grd', '.tif', '.img'")
}

weather <- FALSE
if (temp == TRUE) {
  avg_temperature_coefficient <- import_raster(file = "avg_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
  low_temperature_coefficient <- import_raster(file = "low_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
  high_temperature_coefficient <- import_raster(file = "high_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
  
  if (!(raster::extent(infected) == raster::extent(avg_temperature_coefficient))) {
    return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  }
  
  if (!(raster::xres(infected) == raster::xres(avg_temperature_coefficient) && raster::yres(infected) == raster::yres(avg_temperature_coefficient))) {
    return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  }
  
  if (!(raster::compareCRS(infected, avg_temperature_coefficient))) {
    return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  }
  
  weather <- TRUE
  avg_weather_coefficient_stack <- avg_temperature_coefficient
  low_weather_coefficient_stack <- low_temperature_coefficient
  high_weather_coefficient_stack <- high_temperature_coefficient
  if (precip ==TRUE){
    precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
    
    if (!(raster::extent(infected) == raster::extent(precipitation_coefficient))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(precipitation_coefficient) && raster::yres(infected) == raster::yres(precipitation_coefficient))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, precipitation_coefficient))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
  }
} else if(precip == TRUE){
  precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
  
  if (!(raster::extent(infected) == raster::extent(precipitation_coefficient))) {
    return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  }
  
  if (!(raster::xres(infected) == raster::xres(precipitation_coefficient) && raster::yres(infected) == raster::yres(precipitation_coefficient))) {
    return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  }
  
  if (!(raster::compareCRS(infected, precipitation_coefficient))) {
    return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  }
  
  weather <- TRUE
  weather_coefficient_stack <- precipitation_coefficient
}

if (weather == TRUE){
  avg_weather_coefficient_stack[is.na(avg_weather_coefficient_stack)] <- 0
  avg_weather_coefficient <- list(raster::as.matrix(avg_weather_coefficient_stack[[1]]))
  low_weather_coefficient_stack[is.na(low_weather_coefficient_stack)] <- 0
  low_weather_coefficient <- list(raster::as.matrix(low_weather_coefficient_stack[[1]]))
  high_weather_coefficient_stack[is.na(high_weather_coefficient_stack)] <- 0
  high_weather_coefficient <- list(raster::as.matrix(high_weather_coefficient_stack[[1]]))
  
  for(i in 2:number_of_time_steps) {
    avg_weather_coefficient[[i]] <- raster::as.matrix(avg_weather_coefficient_stack[[i]])
    low_weather_coefficient[[i]] <- raster::as.matrix(low_weather_coefficient_stack[[i]])
    high_weather_coefficient[[i]] <- raster::as.matrix(high_weather_coefficient_stack[[i]])
  }
} else {
  weather_coefficient <- host
  raster::values(weather_coefficient) <- 1
  weather_coefficient <- list(raster::as.matrix(weather_coefficient))
}

if (management == TRUE  && !file.exists(treatments_file)) {
  return("Treatments file does not exist")
}

if (management == TRUE  && !(raster::extension(treatments_file) %in% c(".grd", ".tif", ".img"))) {
  return("Treatments file is not one of '.grd', '.tif', '.img'")
}

if (management == TRUE) {
  
  treatment_stack <- raster::stack(treatments_file)
  treatment_stack[is.na(treatment_stack)] <- 0
  
  if (!(raster::extent(infected) == raster::extent(treatment_stack))) {
    return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  }
  
  if (!(raster::xres(infected) == raster::xres(treatment_stack) && raster::yres(infected) == raster::yres(treatment_stack))) {
    return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  }
  
  if (!(raster::compareCRS(infected, treatment_stack))) {
    return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  }
  
  treatment_maps <- list(raster::as.matrix(treatment_stack[[1]]))
  if (raster::nlayers(treatment_stack) >= 2) {
    for(i in 2:raster::nlayers(treatment_stack)) {
      treatment_maps[[i]] <- raster::as.matrix(treatment_stack[[i]])
    }
  }
  treatment_years = treatment_years
} else {
  treatment_map <- host
  raster::values(treatment_map) <- 0
  treatment_maps = list(raster::as.matrix(treatment_map))
}

ew_res <- raster::xres(susceptible)
ns_res <- raster::yres(susceptible)

mortality_tracker <- infected
raster::values(mortality_tracker) <- 0

infected <- raster::as.matrix(infected)
susceptible <- raster::as.matrix(susceptible)
total_plants <- raster::as.matrix(total_plants)
mortality_tracker <- raster::as.matrix(mortality_tracker)
mortality <- mortality_tracker

random_seed <- round(stats::runif(1, 1, 1000000))
reproductive_rate <- 1.9
short_distance_scale <- 28


temperature <- avg_temperature
weather_coefficient <- avg_weather_coefficient

data <- pops_model(random_seed = random_seed,
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
                   long_distance_scale = long_distance_scale,
                   wind_dir = wind_dir, kappa = kappa)


rm(data)
rm(random_seed)
rm(reproductive_rate)
rm(short_distance_scale)
rm(json_case_study)
rm(case_study)
rm(case_study_id)
rm(high_weather_coefficient_stack)
rm(high_temperature_stack)
rm(low_weather_coefficient_stack)
rm(low_temperature_stack)
rm(avg_weather_coefficient_stack)
rm(avg_temperature_stack)
rm(avg)
rm(high_temperature_coefficient)
rm(avg_temperature_coefficient)
rm(low_temperature_coefficient)
rm(treatment_map)
rm(low_weather_coefficient)
rm(high_weather_coefficient)
rm(avg_weather_coefficient)
gcs_save_image(file = "casestudy2.Rdata", bucket = "testers-pops")
gcs_load(file = "casestudy2.Rdata", bucket = "testers-pops")
# gcs_load(file = ".Rdata", bucket = "testers-pops")
