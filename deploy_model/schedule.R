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
library(aws.s3)

#' Plot out data from the iris dataset
#' 
#' @param case_study_id required to get the proper case study information
#' @param session_id required to set the proper session id
#' @param run_collection_id required to set the proper run collection id
#' @param run_id required to set the proper run id
#' @get /status
modelapi <- function(case_study_id, session_id, run_collection_id, run_id) {
  
  options(digits = 6)
  run_id <- as.numeric(run_id)
  # json_run <- httr::GET(paste("https://pops-model.org/api/run/", run_id, "/?format=json", sep = ""))
  json_run <- httr::GET(paste("http://127.0.0.1/api/run_detail/", run_id, "/?format=json", sep = ""))
  run <- httr::content(json_run)
  run$status <- "READING DATA"
  # httr::PUT(url = paste("https://pops-model.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  httr::PUT(url = paste("http://127.0.0.1/api/run_detail/", run_id, "/", sep = ""), body = run, encode = "json")
  case_study_id <- as.numeric(case_study_id)
  session_id <- as.numeric(session_id)
  run_collection_id <- as.numeric(run_collection_id)
  # json_run_collection <- httr::GET(paste("https://pops-model.org/api/run_collection/", run_collection_id, "/?format=json", sep = ""))
  json_run_collection <- httr::GET(paste("http://127.0.0.1/api/run_collection/", run_collection_id, "/?format=json", sep = ""))
  run_collection <- httr::content(json_run_collection)
  # json_session <- httr::GET(paste("https://pops-model.org/api/session/", session_id, "/?format=json", sep = ""))
  json_session <- httr::GET(paste("http://127.0.0.1/api/session/", session_id, "/?format=json", sep = ""))
  session <- httr::content(json_session)
  
  ### this needs to change to the S3 bucket
  # googleCloudStorageR::gcs_load(file = paste("casestudy", case_study_id, ".Rdata", sep = ""), bucket = "test_pops_staging")
  
  config$end_date <- session$final_date
  ### potentially ignore this
  natural_distance_scale <- as.numeric(session$distance_scale)
  reproductive_rate <- as.numeric(session$reproductive_rate)
  ## need to pull this from the polygons now the data
  efficacy = run_collection$efficacy
  treatment_month <- session$management_month
  # susceptible_start <- susceptible
  # infected_start <- raster(infected, host)
  if (is.null(run$management_polygons) || class(run$management_polygons) != "list") {
    treatments_file <- ""
    treatment_years <- c(0)
    management <- FALSE
  } else if (length(run$management_polygons) >= 1) {
    treatments_file <- run$management_polygons
    treatments <- geojsonio::as.json(treatments_file)
    treatments <- geojson::as.geojson(treatments)
    treatments <- geojsonio::geojson_sp(treatments)
    treatments <- spTransform(treatments, CRS = crs(host))
    treatment_map <- raster::rasterize(treatments, host, fun = "last", getCover = TRUE)
    # treatment_map[treatment_map > 1] <- 1 ## rasterize gives rasters in each polygon the polygon id value need to set those to 1
    treatment_map[is.na(treatment_map)] <- 0
    treatment_map <- treatment_map * (efficacy / 100)
    treatment_maps <- list(raster::as.matrix(treatment_map))
    treatment_years <- c(run$steering_year)
    management <- TRUE
  }
  
  if (session$weather == "GOOD") {
    temperature <- high_temperature
  } else if (session$weather == "BAD") {
    temperature <- low_temperature
  } else {
    temperature <- temperature
  }
  
  
  if (!is.null(run$steering_year)) {
    
    if (run$steering_year > start_time) {
      infected_r <- flyio::import_raster(file = paste("infected_", case_study_id, "_", run_collection$second_most_recent_run, ".tif", sep = ""), data_source = "gcs", bucket = "test_pops_staging")
      susceptible_r <- flyio::import_raster(file = paste("susceptible_", case_study_id, "_", run_collection$second_most_recent_run,".tif", sep = ""), data_source = "gcs", bucket = "test_pops_staging")
      infected <- raster::as.matrix(infected_r)
      susceptible <- raster::as.matrix(susceptible_r)
      if (is.null(run$management_polygons) || class(run$management_polygons) != "list") {
        
      } else {
        run$management_cost <- round(sum(treatment_map[treatment_map > 0 & (infected_r > 0 | susceptible_r > 0)]) * xres(treatment_map) * yres(treatment_map) * as.numeric(run_collection$cost_per_meter_squared), digits = 2)
      }
    } else {
      if (is.null(run$management_polygons) || class(run$management_polygons) != "list") {
        
      } else {
        run$management_cost <- round(sum(treatment_map[treatment_map > 0 & (host > 0)]) * xres(treatment_map) * yres(treatment_map) * as.numeric(run_collection$cost_per_meter_squared), digits = 2)
      }
    }
    
    years_difference <- run$steering_year - start_time
    start_time <- run$steering_year
    years_move <- years_difference + 1
    if (use_lethal_temperature) {
      temperature <- temperature[years_move:length(temperature)]
    }
    time_step_move <- years_difference * steps_in_year + 1
    weather_coefficient <- weather_coefficient[time_step_move:length(weather_coefficient)]
  }
  
  core_count <- 10
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  run$status <- "RUNNING MODEL"
  httr::PUT(url = paste("https://pops-model.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  years <- seq(start_time, end_time, 1)
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  
  infected_stack <- foreach::foreach(i = 1:10, .combine = c, .packages = c("raster", "PoPS"), .export = ls(globalenv())) %dopar% {
    random_seed <- round(stats::runif(1, 1, 1000000))
    
    data <- PoPS::pops_model(random_seed = config$random_seed, 
                             use_lethal_temperature = config$use_lethal_temperature, 
                             lethal_temperature = config$lethal_temperature, 
                             lethal_temperature_month = config$lethal_temperature_month,
                             infected = config$infected,
                             exposed = config$exposed,
                             susceptible = config$susceptible,
                             total_populations  = config$total_populations,
                             mortality_on = config$mortality_on,
                             mortality_tracker = config$mortality_tracker,
                             mortality = config$mortality,
                             quarantine_areas = config$quarantine_areas,
                             treatment_maps = config$treatment_maps,
                             treatment_dates = config$treatment_dates,
                             pesticide_duration = config$pesticide_duration,
                             resistant = config$resistant,
                             use_movements = config$use_movements,
                             movements = config$movements,
                             movements_dates = config$movements_dates,
                             weather = config$weather,
                             temperature = config$temperature,
                             weather_coefficient = config$weather_coefficient,
                             ew_res = config$ew_res,
                             ns_res = config$ns_res,
                             num_rows = config$num_rows,
                             num_cols = config$num_cols,
                             time_step = config$time_step,
                             reproductive_rate = config$reproductive_rate,
                             spatial_indices = config$spatial_indices,
                             mortality_rate = config$mortality_rate,
                             mortality_time_lag = config$mortality_time_lag,
                             season_month_start = config$season_month_start,
                             season_month_end = config$season_month_end,
                             start_date = config$start_date,
                             end_date = config$end_date,
                             treatment_method = config$treatment_method,
                             natural_kernel_type = config$natural_kernel_type,
                             anthropogenic_kernel_type =
                               config$anthropogenic_kernel_type,
                             use_anthropogenic_kernel =
                               config$use_anthropogenic_kernel,
                             percent_natural_dispersal =
                               config$percent_natural_dispersal,
                             natural_distance_scale =
                               config$natural_distance_scale,
                             anthropogenic_distance_scale =
                               config$anthropogenic_distance_scale,
                             natural_dir = config$natural_dir,
                             natural_kappa = config$natural_kappa,
                             anthropogenic_dir = config$anthropogenic_dir,
                             anthropogenic_kappa = config$anthropogenic_kappa,
                             output_frequency = config$output_frequency,
                             output_frequency_n = config$output_frequency_n,
                             quarantine_frequency = config$quarantine_frequency,
                             quarantine_frequency_n = config$quarantine_frequency_n,
                             use_quarantine = config$use_quarantine,
                             spreadrate_frequency = config$spreadrate_frequency,
                             spreadrate_frequency_n = config$spreadrate_frequency_n,
                             use_spreadrates = config$use_spreadrates,
                             model_type_ = config$model_type,
                             latency_period = config$latency_period,
                             generate_stochasticity =
                               config$generate_stochasticity,
                             establishment_stochasticity =
                               config$establishment_stochasticity,
                             movement_stochasticity = config$movement_stochasticity,
                             deterministic = config$deterministic,
                             establishment_probability =
                               config$establishment_probability,
                             dispersal_percentage = config$dispersal_percentage)
    
    comp_years <- raster::stack(lapply(1:length(data$infected_before_treatment), function(i) host))
    susceptible_runs <- raster::stack(lapply(1:length(data$infected_before_treatment), function(i) host))
    
    for (q in 1:raster::nlayers(comp_years)) {
      comp_years[[q]] <- data$infected[[q]]
      susceptible_runs[[q]] <- data$susceptible[[q]]
    }
    
    number_infected <- data$number_infected
    spread_rate <- data$rates
    infected_area <- data$area_infected
    single_run <- comp_years
    comp_years <- raster::reclassify(comp_years, rclmat)
    comp_years[is.na(comp_years)] <- 0
    infected_stack <- comp_years
    data <- list(single_run, infected_stack, number_infected, susceptible_runs, infected_area, spread_rate)
    
  }
  
  stopCluster(cl)
  single_runs <- infected_stack[seq(1,length(infected_stack),6)]
  probability_runs <- infected_stack[seq(2,length(infected_stack),6)]
  number_infected_runs <- infected_stack[seq(3,length(infected_stack),6)]
  susceptible_runs <- infected_stack[seq(4,length(infected_stack),6)]
  area_infected_runs <- infected_stack[seq(5,length(infected_stack),6)]
  spread_rate_runs <- infected_stack[seq(6,length(infected_stack),6)]
  
  prediction <- probability_runs[[1]]
  prediction[prediction > 0] <- 0
  infected_area <- data.frame(t(years))
  infected_number <- data.frame(t(years))
  west_rates <- data.frame(t(years))
  east_rates <- data.frame(t(years))
  south_rates <- data.frame(t(years))
  north_rates <- data.frame(t(years))
  max_values <- data.frame(t(years))
  
  for (i in 1:length(probability_runs)) {
    prediction <- prediction + probability_runs[[i]]
    infected_number[i,] <- number_infected_runs[[i]]
    infected_area[i,] <- area_infected_runs[[i]]
    rates <- do.call(rbind, spread_rate_runs[[i]])
    west_rates[i,] <- rates[,4]
    east_rates[i,] <- rates[,3]
    south_rates[i,] <- rates[,2]
    north_rates[i,] <- rates[,1]
    max_values[i,] <- raster::maxValue(single_runs[[i]])
  }
  
  
  probability <- (prediction/(length(probability_runs))) * 100
  
  infected_areas <- round(sapply(infected_area, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  number_infecteds <- round(sapply(infected_number, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  west_rate <- round(sapply(west_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  east_rate <- round(sapply(east_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  south_rate <- round(sapply(south_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  north_rate <- round(sapply(north_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  which_median <- function(x) raster::which.min(abs(x - median(x)))
  
  median_run_index <- which_median(infected_number[[1]])
  
  single_run <- single_runs[[median_run_index]]
  susceptible_run <- susceptible_runs[[median_run_index]]
  
  run$status <- "WRITING DATA"
  httr::PUT(url = paste("https://pops-model.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  
  if (run_id == session$default_run) {
    max_value <- round(sapply(max_values, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
    max_value_out <- max_value[1,ncol(max_value)]
    session$max_value <- max_value_out
    httr::PUT(url = paste("https://pops-model.org/api/session/", session_id, "/", sep = ""), body = session, encode = "json")
    
  }
  
  single_run_out <- single_run[[1]]
  susceptible_run_out <- susceptible_run[[1]]
  
  flyio::flyio_auth(auth_list = c("GCS_AUTH_FILE"))
  flyio::export_raster(x = single_run_out, file = paste("infected_", case_study_id , "_", run_id ,".tif", sep = ""), data_source = "gcs", bucket = "test_pops_staging")
  flyio::export_raster(x = susceptible_run_out, file = paste("susceptible_", case_study_id , "_", run_id ,".tif", sep = ""), data_source = "gcs", bucket = "test_pops_staging", overwrite = TRUE)
  
  # core_count <- 10
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  stuff <- foreach::foreach(q = 1:length(years), .packages =c("raster", "geojsonio", "httr"), .export = ls(globalenv())) %dopar% {
    if (q == 1) {
      number_infected <- infected_number[median_run_index, q]
      area_infected <- infected_area[median_run_index, q]
      if (is.nan(west_rates[median_run_index,q])) {
        west_rate_r <- 0
      } else {
        west_rate_r <- west_rates[median_run_index,q]
      }
      if (is.nan(east_rates[median_run_index,q])) {
        east_rate_r <- 0
      } else {
        east_rate_r <- east_rates[median_run_index,q]
      }
      if (is.nan(south_rates[median_run_index,q])) {
        south_rate_r <- 0
      } else {
        south_rate_r <- south_rates[median_run_index,q]
      }
      if (is.nan(north_rates[median_run_index,q])) {
        north_rate_r <- 0
      } else {
        north_rate_r <- north_rates[median_run_index,q]      
      }
    } else if (q > 1) {
      number_infected <- number_infecteds[1, q]
      area_infected <- infected_areas[1, q]
      if (is.nan(west_rate[1,q])) {
        west_rate_r <- 0
      } else {
        west_rate_r <- west_rate[1,q]      }
      if (is.nan(east_rate[1,q])) {
        east_rate_r <- 0
      } else {
        east_rate_r <- east_rate[1,q]
      }
      if (is.nan(south_rate[1,q])) {
        south_rate_r <- 0
      } else {
        south_rate_r <- south_rate[1,q]      }
      if (is.nan(north_rate[1,q])) {
        north_rate_r <- 0
      } else {
        north_rate_r <- north_rate[1,q]     
      }
    }
    
    year <- years[q]
    
    single_map <- single_run[[q]]
    single_map <- as.integer(single_map)
    single_map[single_map <= 0] <- NA
    names(single_map) <- "outputs"
    if (cellStats(single_map, stat = 'sum') == 0) {
      single_map[infected_start > 0] <- 0
    }
    # single_map <- projectRaster(single_map, crs = CRS("+proj=longlat +datum=WGS84"), method = "ngb")
    single_map <- raster::rasterToPolygons(single_map, n = 4, digits = 4, dissolve = T, na.rm = TRUE)
    storage.mode(single_map$outputs) <- "integer"
    single_map <- geojsonio::geojson_list(single_map, precision = 4, convert_wgs84 = TRUE, geometry = "polygon")
    class(single_map) <- "list"
    
    spread_map <- probability[[q]]
    spread_map <- as.integer(spread_map)
    spread_map[spread_map <= 0] <- NA
    names(spread_map) <- "outputs"
    if (cellStats(spread_map, stat = 'sum') == 0) {
      spread_map[infected_start > 0] <- 0
    }
    # spread_map <- projectRaster(spread_map, crs = CRS("+proj=longlat +datum=WGS84"), method = "ngb")
    spread_map <- raster::rasterToPolygons(spread_map, n = 4, dissolve = T, na.rm = TRUE)
    storage.mode(spread_map$outputs) <- "integer"
    spread_map <- geojsonio::geojson_list(spread_map, precision = 4, convert_wgs84 = TRUE, geometry = "polygon")
    class(spread_map) <- "list"
    
    outs <- list()
    outs$run <- run_id
    outs$number_infected <- number_infected
    outs$infected_area <- area_infected
    outs$year <- year
    outs$single_spread_map <- single_map
    outs$probability_map <- spread_map
    outs$susceptible_map <- "null"
    outs$escape_probability <- 0
    spreadrate <- list()
    spreadrate$west_rate <- west_rate_r
    spreadrate$east_rate <- east_rate_r
    spreadrate$south_rate <- south_rate_r
    spreadrate$north_rate <- north_rate_r
    
    timetoboundary <- list()
    timetoboundary$west_time <- 0
    timetoboundary$east_time <- 0
    timetoboundary$south_time <- 0
    timetoboundary$north_time <- 0
    
    distancetoboundary <- list()
    distancetoboundary$west_distance <- 0
    distancetoboundary$east_distance <- 0
    distancetoboundary$south_distance <- 0
    distancetoboundary$north_distance <- 0
    outs$timetoboundary <- timetoboundary
    outs$spreadrate <- spreadrate
    outs$distancetoboundary <- distancetoboundary
    post_code <- httr::POST(url = "https://pops-model.org/api/output/", body = outs, encode = "json")
    if (post_code$status_code == 201) {
      run$status <- "SUCCESS"
    } else {
      run$status <- "FAILED"
    }
    stuff <- run$status
  }
  
  stopCluster(cl)
  run$status <- stuff[[1]]
  
  if (run$status == "SUCCESS") {
    httr::PUT(url = paste("https://pops-model.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  } else {
    httr::PUT(url = paste("https://pops-model.org/api/run/", run_id, "/", sep = ""), body = run, encode = "json")
  }
  
  if (run_collection$default == TRUE || run$steering_year == end_time) {
    if (run$status == "SUCCESS") {
      run_collection$status <- "SUCCESS"
      httr::PUT(url = paste("https://pops-model.org/api/run_collection/", run_collection_id, "/", sep = ""), body = run_collection, encode = "json")
    } else {
      run_collection$status <- "FAILED"
      httr::PUT(url = paste("https://pops-model.org/api/run_collection/", run_collection_id, "/", sep = ""), body = run_collection, encode = "json")
    }
  }
  
  status <- run$status
  
}