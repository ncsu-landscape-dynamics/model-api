# devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "terra")
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

case_studey_setup <- function(case_study_id, bucket = "") {
  
  case_study_id <- as.numeric(case_study_id)
  # change api_url to access either dev, staging, or production Database API.
  # api_url <- "https://popsmodel.org/api/case_study/"
  api_url <- "http://127.0.0.1:8000/api/"
  json_case_study <- httr::GET(paste(api_url, "case_study/", case_study_id ,"/?format=json", sep = ""))
  case_study <- httr::content(json_case_study)
  
  config <- c()
  random_seed <- round(stats::runif(1, 1, 1000000))
  config$random_seed <- random_seed
  
  ## marked for later for handling multiple pests
  # for (i in seq_len(length(case_study$pest_set))) {
  #   
  # }
  
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
  if (config$parameter_means[3] < 1.0) {
    config$use_anthropogenic_kernel <- TRUE
  } else {
    config$use_anthropogenic_kernel <- FALSE
  }
  
  # set up temperature and precipitation data
  config$temp <- case_study$pest_set[[1]]$weather$temp_on
  if (config$temp) {
    config$temperature_coefficient_file <- 
      case_study$pest_set[[1]]$weather$temperature$temperature_data
    config$temperature_coefficient_file <-
      stringr::str_split(config$temperature_coefficient_file, pattern = ".com/")[[1]][2]
  } else {
    config$temperature_coefficient_file <- ""
  }

  config$precip <- case_study$pest_set[[1]]$weather$precipitation_on
  if (config$precip) {
    config$precipitation_coefficient_file <- 
      case_study$pest_set[[1]]$weather$precipitation$precipitation_data
    config$precipitation_coefficient_file <- 
      stringr::str_split(config$precipitation_coefficient_file, pattern = ".com/")[[1]][2]
  } else {
    config$precipitation_coefficient_file <- ""
  }

  
  # determine model type and latency period if SEI
  config$model_type <- case_study$pest_set[[1]]$model_type
  if (config$model_type == 'SI') {
    config$latency_period <- 0
  } else {
    config$latency_period <- 
      mean(case_study$pest_set[[1]]$latencyperiod$minimum, case_study$pest_set[[1]]$latencyperiod$maximum)
  }
  
  # set up model timing variables
  config$time_step <- case_study$time_step_unit
  config$time_step_n <- case_study$time_step_n
  if (case_study$pest_set[[1]]$weather$seasonality_on) {
    config$season_month_start <- 
      case_study$pest_set[[1]]$weather$seasonality$first_month 
    config$season_month_end <- 
      case_study$pest_set[[1]]$weather$seasonality$last_month 
  } else {
    config$season_month_start <- 1
    config$season_month_end <- 12
  }
  config$start_date <- case_study$first_forecast_date
  config$end_date <- case_study$last_forecast_date
  
  # Set up lethal temperature
  config$use_lethal_temperature <- case_study$pest_set[[1]]$weather$lethal_temp_on
  
  if (config$use_lethal_temperature) {
    config$lethal_temperature <- 
      as.numeric(case_study$pest_set[[1]]$weather$lethaltemperature$value)
    config$lethal_temperature_month <- 
      case_study$pest_set[[1]]$weather$lethaltemperature$month
    config$temperature_file <- 
      case_study$pest_set[[1]]$weather$lethaltemperature$lethal_temperature_data
    config$temperature_file <- 
      stringr::str_split(config$temperature_file, pattern = ".com/")[[1]][2]
  } else {
    config$lethal_temperature <- -20
    config$lethal_temperature_month <- 1
    config$temperature_file <- ""
  }
  
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
  ## still need to add in check for treatment and treatment files (probably 
  ## something that isn't needed)
  config$management <- case_study$pest_set[[1]]$use_treatment
  if (config$management) {
    for (i in seq_len(length(case_study$pest_set))) {
      
    }
    config$treatment_dates <- c(0)
    config$treatments_file <- case_study$pest_set[[1]]$priortreatment
    config$treatment_method <- "ratio"
    # config$pesticide_duration <- pesticide_duration
    # config$pesticide_efficacy <- pesticide_efficacy
  } else {
    config$treatment_dates <- c(0)
    config$treatments_file <- ""
    config$treatment_method <- "ratio"
    config$pesticide_duration <- c(0)
    config$pesticide_efficacy <- c(0)
  }

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
  
  if (case_study$pest_set[[1]]$weather$wind_on == FALSE) {
    config$natural_dir <- "NONE"
    config$natural_kappa <- 0
  } else {
    config$natural_dir <- case_study$pest_set[[1]]$weather$wind$wind_direction
    config$natural_kappa <- case_study$pest_set[[1]]$weather$wind$kappa
  }
  if (is.null(case_study$pest_set[[1]]$anthropogenicdirection)) {
    config$anthropogenic_dir <- "NONE"
    config$anthropogenic_kappa <- 0
  } else {
    # need to test whether these are correct
    config$anthropogenic_dir <- 
      case_study$pest_set[[1]]$anthropogenicdirection$direction
    config$anthropogenic_kappa <- 
      case_study$pest_set[[1]]$anthropogenicdirection$kappa
  }

  # set output frequency
  config$output_frequency <- case_study$output_frequency_unit
  config$output_frequency_n <- case_study$output_frequency_n
  # read in movement variables
  
  config$use_movements <- case_study$use_movements
  if (config$use_movements) {
    config$movements_file <- 
      case_study$pest_set[[1]]$pesthostinteraction_set[[1]]$clippedhostmovement
  } else {
    config$movements_file <- ""
  }
  
  # determine if model should have the populations start off exposed
  config$start_exposed <- case_study$start_exposed
  
  # set case study stochasticity off by default and allow user to adjust in
  # session settings
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
  config$function_name <- "casestudy_creation"
  config$failure <- NULL
  config$use_s3 <- TRUE
  config$bucket <- bucket
  
  config <- configuration(config)
  config$crs <- terra::crs(config$host)
  config$xmax <- terra::xmax(config$host)
  config$xmin <- terra::xmin(config$host)
  config$ymax <- terra::ymax(config$host)
  config$ymin <- terra::ymin(config$host)
  
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
                           reproductive_rate = config$reproductive_rate[1],
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
                             config$percent_natural_dispersal[1],
                           natural_distance_scale =
                             config$natural_distance_scale[1],
                           anthropogenic_distance_scale =
                             config$anthropogenic_distance_scale[1],
                           natural_dir = config$natural_dir,
                           natural_kappa = config$natural_kappa[1],
                           anthropogenic_dir = config$anthropogenic_dir,
                           anthropogenic_kappa = config$anthropogenic_kappa[1],
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
  rm(random_seed)
  rm(json_case_study)
  rm(case_study)
  rm(run_collection_id)
  rm(bucket)
  rm(run_id)
  rm(session_id)
  save.image("case_study.RData")
  
  ## save to dashboard project bucket
  casestudy_cs <- httr::upload_file("case_study.RData")
  httr::PUT(url = paste(api_url, "case_study_r_data/", case_study_id, "/", sep = ""), body = list(r_data = casestudy_cs))
}


