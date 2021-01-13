# devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "feature/spread_rate")
library(readr)
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


## 2 is SLF, 3 is SOD EU1, 5 is SOD NA1

case_studey_setup <- function(case_study_id) {
  
  case_study_id <- as.numeric(case_study_id)
  # json_case_study <- httr::GET(paste("https://popsmodel.org/api/case_study/", case_study_id ,"/?format=json", sep = ""))
  json_case_study <- httr::GET(paste("http://127.0.0.1/api/case_study/", case_study_id ,"/?format=json", sep = ""))
  case_study <- httr::content(json_case_study)
  
  config <- c()
  random_seed <- round(stats::runif(1, 1, 1000000))
  config$random_seed <- random_seed
  
  for (i in seq_len(length(case_study$pest_set))) {
    
  }
  
  # set up host data and initial infected conditions data
  config$infected_file <- case_study$pest_set[[1]]$infestation$user_file
  config$infected_file <- 
    stringr::str_split(config$infected_file, pattern = ".com/")[[1]][2]
  config$host_file <- 
    case_study$pest_set[[1]]$pesthostinteraction_set[[1]]$clippedhostlocation$raster_map
  config$host_file <- 
    stringr::str_split(config$host_file, pattern = ".com/")[[1]][2]
  config$total_populations_file <- case_study$allpopulationsdata$user_file
  config$total_populations_file <- 
    stringr::str_split(config$total_populations_file, pattern = ".com/")[[1]][2]
  config$parameter_means <- unlist(case_study$pest_set[[1]]$parameters$means)
  config$parameter_cov_matrix <- 
    matrix(unlist(case_study$pest_set[[1]]$parameters$covariance_matrix), nrow = 6, ncol = 6)
  # set up temperature and precipitation data
  config$temp <- case_study$weather$temp_on
  config$temperature_coefficient_file <- temperature_coefficient_file
  config$precip <- case_study$weather$precipitation_on
  config$precipitation_coefficient_file <- precipitation_coefficient_file
  # determine model type and latency period if SEI
  config$model_type <- case_study$pest_set[[1]]$model_type
  if (config$model_type == 'SI') {
    config$latency_period <- 0
  } else {
    config$latency_period <- case_study$pest_set[[1]]$latencyperiod
  }
  
  # set up model timing variables
  config$time_step <- case_study$time_step_unit
  config$time_step_n <- case_study$time_step_n
  if (case_study$weather$seasonality_on) {
    config$season_month_start <- case_study$weather$seasonality$first_month 
    config$season_month_end <- case_study$weather$seasonality$last_month 
  } else {
    config$season_month_start <- 1
    config$season_month_end <- 12
  }
  config$start_date <- case_study$first_forecast_date
  config$end_date <- case_study$last_forecast_date
  ## from old setup
  start_time <- case_study$end_year + 1
  end_time <- case_study$future_years
  # Set up lethal temperature
  config$use_lethal_temperature <- case_study$weather$lethal_temp_on
  
  if (config$use_lethal_temperature) {
    config$lethal_temperature <- as.numeric(case_study$weather$lethaltemperature$value)
    config$lethal_temperature_month <- case_study$weather$lethaltemperature$month
  } else {
    config$lethal_temperature <- -20
    config$lethal_temperature_month <- 1
  }
  config$temperature_file <- case_study$weather$temperature_file
  
  # set up mortality
  config$mortality_on <- 
    case_study$pest_set[[1]]$pesthostinteraction_set[[1]]$mortality_on
  if (!config$mortality_on) {
    config$mortality_rate <- 0
    config$mortality_time_lag <- 0
  } else {
    config$mortality_rate <- 
      as.numeric(case_study$pest_set[[1]]$pesthostinteraction_set[[1]]$mortality$rate) 
    config$mortality_time_lag <- 
      case_study$pest_set[[1]]$pesthostinteraction_set[[1]]$mortality$time_lag 
  }
  
  # set up management variables
  config$management <- case_study$pest_set[[1]]$use_treatment
  if (config$management) {
    config$treatment_dates <- c(0)
    config$treatments_file <- case_study$pest_set[[1]]$priortreatment
    config$treatment_method <- "ratio"
  } else {
    config$treatment_dates <- c(0)
    config$treatments_file <- ""
    config$treatment_method <- "ratio"
  }
  config$treatment_dates <- c(0)
  config$treatments_file <- ""
  config$treatment_method <- "ratio"
  ## still need to add in check for treatment and treatment files (probably 
  ## something that isn't needed)
  
  # Setup dispersal parameters
  config$natural_kernel_type <- case_study$pest_set[[1]]$natural_dispersal_type
  config$anthropogenic_kernel_type <- 
    case_study$pest_set[[1]]$anthropogenic_dispersal_type
  if (config$natural_kernel_type == "CAUCHY") {
    config$natural_kernel_type <- "cauchy"
  } else if (config$natural_kernel_type == "EXPONENTIAL") {
    config$natural_kernel_type <- "exponential"
  }
  
  if (config$anthropogenic_kernel_type == "CAUCHY") {
    config$anthropogenic_kernel_type <- "cauchy"
  } else if (config$anthropogenic_kernel_type == "EXPONENTIAL") {
    config$anthropogenic_kernel_type <- "exponential"
  }
  
  if (case_study$weather$wind_on == FALSE) {
    config$natural_dir <- "NONE"
    config$natural_kappa <- 0
  } else {
    config$natural_dir <- case_study$weather$wind$wind_direction
    config$natural_kappa <- case_study$weather$wind$kappa
  }
  if (length(case_study$pest_set[[1]]$anthropogenicdistance_set) == 0) {
    config$anthropogenic_dir <- "NONE"
    config$anthropogenic_kappa <- 0
  } else {
    config$anthropogenic_dir <- case_study$weather$wind
    config$anthropogenic_kappa <- case_study$weather$wind
  }
  # set up pesticide duration and efficacy if used
  config$pesticide_duration <- pesticide_duration
  config$pesticide_efficacy <- pesticide_efficacy
  # set output frequency
  config$output_frequency <- case_study$output_frequency_unit
  config$output_frequency_n <- case_study$output_frequency_n
  # read in movement variables
  
  config$use_movements <- case_study$use_movements
  if (config$use_movements) {
    config$movements_file <- movements_file
  } else {
    config$movements_file <- ""
  }
  
  # determine if model should have the populations start off exposed
  config$start_exposed <- case_study$start_exposed
  # determine whether or not to use stochasticity
  config$generate_stochasticity <- TRUE
  config$establishment_stochasticity <- TRUE
  config$movement_stochasticity <- TRUE
  config$deterministic <- FALSE
  config$establishment_probability <- 0.5
  config$dispersal_percentage <- 0.99
  # set up quarantine and spreadrates
  
  config$use_quarantine <- case_study$pest_set[[1]]$use_quarantine
  if (config$use_quarantine) {
    config$quarantine_areas_file <- case_study$pest_set[[1]]$quarantinelink_set
  } else {
    config$quarantine_areas_file <- ""
  }
  config$use_spreadrates <- case_study$use_spread_rate
  # setup parallel processing
  config$number_of_iterations <- 10
  config$number_of_cores <- 10
  # add function name for use in configuration function to skip
  # function specific configurations namely for validation and
  # calibration.
  config$function_name <- "multirun"
  config$failure <- NULL
  
  config <- configuration(config)
  
  
  ## For SLF
  # infected <- import_raster(file = "initial_infections_2018_single_count_pm_prop.tif", FUN = raster)
  # infected[is.na(infected)] <- 0
  # host <- import_raster(file = "tree_of_heaven_0.50.tif", FUN = raster)
  # host[is.na(host)] <- 0
  # total_plants <- import_raster(file = "total_hosts.tif", FUN = raster)
  # total_plants[is.na(total_plants)] <- 0
  
  ## For SOD EU1
  infected <- import_raster(file = "cum_inf_2019eu.tif", FUN = raster)
  infected[is.na(infected)] <- 0
  host <- import_raster(file = "lide_100m_median_2018.tif", FUN = raster)
  host[is.na(host)] <- 0
  total_plants <- import_raster(file = "lemma_max100m.tif", FUN = raster)
  total_plants[is.na(total_plants)] <- 0
  
  ## For SOD NA1
  # infected <- import_raster(file = "cum_inf_2019.tif", FUN = raster)
  # infected[is.na(infected)] <- 0
  # host <- import_raster(file = "lide_100m_median_2018.tif", FUN = raster)
  # host[is.na(host)] <- 0
  # total_plants <- import_raster(file = "lemma_max100m.tif", FUN = raster)
  # total_plants[is.na(total_plants)] <- 0
  
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
    ## For SLF
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
  
  # if (temp == TRUE  && !file.exists(temperature_coefficient_file)) {
  #   return("Temperature coefficient file does not exist")
  # }
  # 
  # if (temp == TRUE  && !(raster::extension(temperature_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
  #   return("Temperature coefficient file is not one of '.grd', '.tif', '.img'")
  # }
  # 
  # if (precip == TRUE  && !file.exists(precipitation_coefficient_file)) {
  #   return("Precipitation coefficient file does not exist")
  # }
  # 
  # if (precip == TRUE  && !(raster::extension(precipitation_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
  #   return("Precipitation coefficient file is not one of '.grd', '.tif', '.img'")
  # }
  
  weather <- FALSE
  if (temp == TRUE) {
    ## For SLF
    # avg_temperature_coefficient <- import_raster(file = "avg_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
    # low_temperature_coefficient <- import_raster(file = "low_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
    # high_temperature_coefficient <- import_raster(file = "high_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
    
    ## For SOD
    avg_temperature_coefficient <- import_raster(file = "average_weather_2019_2025.tif", FUN = stack)
    # low_temperature_coefficient <- import_raster(file = "low_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
    # high_temperature_coefficient <- import_raster(file = "high_spread_temp_coefficient_slf_2018_2022_pm.tif", FUN = stack)
    
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
    # low_weather_coefficient_stack <- low_temperature_coefficient
    # high_weather_coefficient_stack <- high_temperature_coefficient
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
    # low_weather_coefficient_stack[is.na(low_weather_coefficient_stack)] <- 0
    # low_weather_coefficient <- list(raster::as.matrix(low_weather_coefficient_stack[[1]]))
    # high_weather_coefficient_stack[is.na(high_weather_coefficient_stack)] <- 0
    # high_weather_coefficient <- list(raster::as.matrix(high_weather_coefficient_stack[[1]]))
    
    for(i in 2:number_of_time_steps) {
      avg_weather_coefficient[[i]] <- raster::as.matrix(avg_weather_coefficient_stack[[i]])
      # low_weather_coefficient[[i]] <- raster::as.matrix(low_weather_coefficient_stack[[i]])
      # high_weather_coefficient[[i]] <- raster::as.matrix(high_weather_coefficient_stack[[i]])
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
  
  if(percent_natural_dispersal == 1.0) {
    use_anthropogenic_kernel = FALSE
  } else if (percent_natural_dispersal < 1.0  && percent_natural_dispersal >= 0.0) {
    use_anthropogenic_kernel = TRUE
  } else {
    return("Percent natural dispersal must be between 0.0 and 1.0")
  }
  
  config$ew_res <- raster::xres(susceptible)
  config$ns_res <- raster::yres(susceptible)
  config$num_cols <- raster::ncol(susceptible)
  config$num_rows <- raster::nrow(susceptible)
  
  mortality_tracker <- infected
  raster::values(mortality_tracker) <- 0
  
  infected <- raster::as.matrix(infected)
  susceptible <- raster::as.matrix(susceptible)
  total_plants <- raster::as.matrix(total_plants)
  mortality_tracker <- raster::as.matrix(mortality_tracker)
  mortality <- mortality_tracker
  
  
  reproductive_rate <- 1.6
  natural_distance_scale <- 242
  
  treatment_month <- 12
  treatment_method <- "ratio"
  
  # temperature <- avg_temperature
  # weather_coefficient <- avg_weather_coefficient
  end_time <- 2021
  
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
  
  
  # weather_coefficient <- weather_coefficient[1:156]
  
  rm(data)
  rm(i)
  rm(random_seed)
  rm(reproductive_rate)
  rm(natural_distance_scale)
  rm(json_case_study)
  rm(case_study)
  # rm(case_study_id)
  rm(high_weather_coefficient_stack)
  rm(high_temperature_stack)
  rm(low_weather_coefficient_stack)
  rm(low_temperature_stack)
  rm(avg_weather_coefficient_stack)
  rm(avg_temperature_stack)
  rm(high_temperature_coefficient)
  rm(avg_temperature_coefficient)
  rm(low_temperature_coefficient)
  rm(treatment_map)
  rm(low_weather_coefficient)
  rm(high_weather_coefficient)
  rm(avg_weather_coefficient)
  
  
  ## save to dashboard project bucket
  case_study$r_data
  gcs_save_image(file = paste("casestudy", case_study_id, ".Rdata", sep = ""), bucket = "test_pops_staging")
  googleCloudStorageR::gcs_load(file = paste("casestudy", case_study_id, ".Rdata", sep = ""), bucket = "test_pops_staging")
}


