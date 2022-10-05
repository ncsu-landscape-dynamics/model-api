library(PoPS)
library(httr)
# library(geojsonio)
# library(protolite)
# library(geojson)
library(rgdal)
library(raster)
# library(jsonlite)
library(sp)
library(devtools)
library(doParallel)
library(foreach)
library(parallel)
library(aws.s3)
library(sf)
library(terra)
library(plumber)
library(stars)
library(geojsonsf)

#' Return the status of a model call
#'
#' @param case_study_id required to get the proper case study information
#' @param session_id required to set the proper session id
#' @param run_collection_id required to set the proper run collection id
#' @param run_id required to set the proper run id
#' @get /status
modelapi <- function(case_study_id, session_id, run_collection_id, run_id) {

  # change api_url to access either dev, staging, or production Database API.
  # api_url <- "http://127.0.0.1:8000/api/"
  api_url <- "https://pops-model.org/api/"
  readRenviron("env")

  setAs("SpatVector", "sf",
        function(from) {
          sf::st_as_sf(as.data.frame(from, geom = TRUE),
                       wkt = "geometry", crs = crs(from))
        }
  )

  st_as_sf.SpatVector <- function(x, ...) {
    sf::st_as_sf(as.data.frame(x, geom = TRUE), wkt = "geometry",
                 crs = crs(x))
  }

  run_id <- as.numeric(run_id)
  json_run <-
    httr::GET(paste(api_url, "run_write/", run_id, "/?format=json", sep = ""))
  run <- httr::content(json_run)
  run2 <- run[c(1:10, 12)]
  httr::PUT(url = paste(api_url, "run_write/", run_id, "/", sep = ""),
            body = run2, encode = "json")

  case_study_id <- as.numeric(case_study_id)
  session_id <- as.numeric(session_id)
  run_collection_id <- as.numeric(run_collection_id)
  json_run_collection <-
    httr::GET(paste(api_url, "run_collection_write/",
                    run_collection_id, "/?format=json", sep = ""))
  run_collection <- httr::content(json_run_collection)
  json_session <-
    httr::GET(paste(api_url, "session_write/",
                    session_id, "/?format=json", sep = ""))
  session <- httr::content(json_session)
  json_case_study <-
    httr::GET(paste(api_url, "case_study/", case_study_id,
                    "/?format=json", sep = ""))
  case_study <- httr::content(json_case_study)

  if (run_collection$second_most_recent_run == "null") {
    run_collection$second_most_recent_run <- NULL
  }

  ## Read in Rdata file
  if (is.null(run$steering_year) ||
      lubridate::year(case_study$first_forecast_date) >= run$steering_year ||
      is.null(run_collection$second_most_recent_run)) {
    run_file <- case_study$r_data
    run_file <- stringr::str_split(run_file, pattern = ".com/")[[1]][2]
  } else {
    prev_run_j <-
      httr::GET(paste(api_url, "run_r_data/",
                      run_collection$second_most_recent_run,
                      "/?format=json", sep = ""))
    prev_run <- httr::content(prev_run_j)
    run_file <- prev_run$r_data
    run_file <- stringr::str_split(run_file, pattern = ".com/")[[1]][2]
  }

  s3load(object = run_file, bucket = "pops-production")

  config$end_date <- session$final_date
  ### potentially add in a section for controlling natural dispersal and
  ### reproductive rate
  # natural_distance_scale <- as.numeric(session$distance_scale)
  # reproductive_rate <- as.numeric(session$reproductive_rate)

  host <-
    terra::rast(nrow = config$num_rows, ncol = config$num_cols,
                xmin = config$xmin, xmax = config$xmax,
                ymin = config$ymin, ymax = config$ymax, crs = config$crs)

  ## need to pull this from the polygons now the data
  if (is.null(run$management_polygons) ||
      class(run$management_polygons) != "list") {
    if (is.null(run$steering_year) ||
        run$steering_year <= lubridate::year(config$start_date)) {
      config$treatment_dates <- config$start_date
    } else {
      config$treatment_dates <- paste(run$steering_year, "-01-01", sep = "")
    }
    config$pesticide_duration <- c(0)
    config$management <- FALSE

  } else if (length(run$management_polygons) >= 1) {
    config$management <- TRUE
    treatments_file <- run$management_polygons
    treatments <- geojsonio::as.json(treatments_file)
    treatments <- sf::st_read(treatments)
    treatments <- sf::st_transform(treatments, crs = config$crs)
    treatments_table <-
      data.frame(treatments[, c("efficacy", "duration", "date")])
    treatments_table <- treatments_table[, c("efficacy", "duration", "date")]
    unique_treatments <- unique(treatments_table)
    pesticide_efficacy <- as.numeric(c())
    pesticide_duration <- as.integer(c())
    treatment_maps <- c()
    treatment_dates <- as.character(c())
    management_costs <- 0
    treatments2 <- vect(treatments)
    treatments2$efficacy <- treatments$efficacy
    treatments2$date <- treatments$date
    treatments2$duration <- treatments$duration
    for (i in seq_len(nrow(unique_treatments))) {
      current_treatments <-
        treatments2[treatments2$date == unique_treatments$date[i] &
                     treatments2$duration == unique_treatments$duration[i] &
                     treatments2$efficacy == unique_treatments$efficacy[i], ]
      treatment_map <-
        terra::rasterize(current_treatments, host, touches = TRUE,
                                        fun = "last", cover = TRUE)
      treatment_map[is.na(treatment_map)] <- 0
      treatment_map[treatment_map > 1] <- 1
      treatment_map <- treatment_map * as.numeric(unique_treatments$efficacy[i])
      treatment_map <- terra::as.matrix(treatment_map, wide = TRUE)
      management_cost <-
        round(sum(treatment_map[treatment_map > 0 &
                                  (config$infected > 0 |
                                     config$susceptible > 0)]) *
                config$ew_res * config$ns_res *
                as.numeric(run_collection$cost_per_meter_squared), digits = 2)
      pesticide_dur <- as.integer(unique_treatments$duration[i])
      if (config$time_step == "month") {
        pesticide_dur <- round(pesticide_dur / 30)
      } else if (config$time_step == "week") {
        pesticide_dur <- round(pesticide_dur / 7)
      } else {
        pesticide_dur <- pesticide_dur
      }
      pesticide_duration[[i]] <- pesticide_dur
      pesticide_efficacy[[i]] <-  as.numeric(unique_treatments$efficacy[i])
      treatment_dates[[i]] <- as.character(unique_treatments$date[i])
      treatment_maps[[i]] <- treatment_map
      management_costs <- management_costs + management_cost
    }
    config$treatment_maps <- treatment_maps
    config$treatment_dates <- treatment_dates
    config$pesticide_duration <- pesticide_duration
    config$pesticide_efficacy <- pesticide_efficacy
  }

  ## set up if statement for this and weather
  # if (config$temp) {
  #   if (session$weather == "GOOD") {
  #     temperature <- high_temperature
  #   } else if (session$weather == "BAD") {
  #     temperature <- low_temperature
  #   } else {
  #     temperature <- temperature
  #   }
  # }

  if (!is.null(run$steering_year)) {
    years_difference <-
      run$steering_year - lubridate::year(config$start_date)
    config$start_date <- paste(run$steering_year, "-01-01", sep = "")
    years_move <- years_difference + 1
    if (config$use_lethal_temperature) {
      config$temperature <- config$temperature[years_move:length(config$temperature)]
    }
    steps_in_year <- floor(config$number_of_time_steps / config$number_of_years)
    time_step_move <- years_difference * steps_in_year + 1
    if (config$weather) {
      config$weather_coefficient <-
        config$weather_coefficient[
          time_step_move:length(config$weather_coefficient)]
    }

  }

  run2$status <- "RUNNING MODEL"
  httr::PUT(url = paste(api_url, "run_write/", run_id, "/", sep = ""),
            body = run2, encode = "json")
  years <-
    seq(lubridate::year(config$start_date), lubridate::year(config$end_date), 1)
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol = 3, byrow = TRUE)

  core_count <- 10
  cl <- makeCluster(core_count)
  registerDoParallel(cl)

  infected_stack <- foreach::foreach(i = seq_len(10),
                                     .combine = c,
                                     .packages = c("PoPS", "stats")) %dopar% {
    config$random_seed <- round(stats::runif(1, 1, 1000000))
    data <- PoPS::pops_model(random_seed = config$random_seed,
                             use_lethal_temperature =
                               config$use_lethal_temperature,
                             lethal_temperature = config$lethal_temperature,
                             lethal_temperature_month =
                               config$lethal_temperature_month,
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
                             reproductive_rate = config$reproductive_rate[i],
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
                               config$percent_natural_dispersal[i],
                             natural_distance_scale =
                               config$natural_distance_scale[i],
                             anthropogenic_distance_scale =
                               config$anthropogenic_distance_scale[i],
                             natural_dir = config$natural_dir,
                             natural_kappa = config$natural_kappa[i],
                             anthropogenic_dir = config$anthropogenic_dir,
                             anthropogenic_kappa =
                               config$anthropogenic_kappa[i],
                             output_frequency = config$output_frequency,
                             output_frequency_n = config$output_frequency_n,
                             quarantine_frequency = config$quarantine_frequency,
                             quarantine_frequency_n =
                               config$quarantine_frequency_n,
                             use_quarantine = config$use_quarantine,
                             spreadrate_frequency = config$spreadrate_frequency,
                             spreadrate_frequency_n =
                               config$spreadrate_frequency_n,
                             use_spreadrates = config$use_spreadrates,
                             model_type_ = config$model_type,
                             latency_period = config$latency_period,
                             generate_stochasticity =
                               config$generate_stochasticity,
                             establishment_stochasticity =
                               config$establishment_stochasticity,
                             movement_stochasticity =
                               config$movement_stochasticity,
                             deterministic = config$deterministic,
                             establishment_probability =
                               config$establishment_probability,
                             dispersal_percentage = config$dispersal_percentage)

    number_infected <- data$number_infected
    spread_rate <- data$rates
    infected_area <- data$area_infected
    quarantine_escape <- data$quarantine_escape
    quarantine_distance <- data$quarantine_escape_distance
    quarantine_direction <- data$quarantine_escape_directions

    single_run <- data$infected
    exposed_runs <- data$exposed
    susceptible_runs <- data$susceptible
    infected_stack <- data$infected

    runs <-
      list(single_run, infected_stack, number_infected, susceptible_runs,
           infected_area, spread_rate, exposed_runs,
           quarantine_escape, quarantine_distance, quarantine_direction)

  }

  stopCluster(cl)
  single_runs <- infected_stack[seq(1, length(infected_stack), 10)]
  probability_runs <- infected_stack[seq(2, length(infected_stack), 10)]
  number_infected_runs <- infected_stack[seq(3, length(infected_stack), 10)]
  susceptible_runs <- infected_stack[seq(4, length(infected_stack), 10)]
  area_infected_runs <- infected_stack[seq(5, length(infected_stack), 10)]
  spread_rate_runs <- infected_stack[seq(6, length(infected_stack), 10)]
  exposed_runs <- infected_stack[seq(7, length(infected_stack), 10)]
  quarantine_escape_runs <- infected_stack[seq(8, length(infected_stack), 10)]
  quarantine_escape_distance_runs <-
    infected_stack[seq(9, length(infected_stack), 10)]
  quarantine_escape_directions_runs <-
    infected_stack[seq(10, length(infected_stack), 10)]

  prediction <- probability_runs[[1]]
  for(w in seq_len(length(prediction))) {
    prediction[[w]] <- 0
  }
  escape_probability <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))
  infected_area <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  infected_number <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  west_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  east_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  south_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  north_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  max_values <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  quarantine_escapes <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))
  quarantine_escape_distances <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))
  quarantine_escape_directions <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))

  for (p in seq_len(length(probability_runs))) {
    for (w in seq_len(length(prediction))) {
      prob <- probability_runs[[p]][[w]]
      max_values[p, w] <- max(prob)
      prob[prob <= 1] <- 0
      prob[prob > 1] <- 1
      prediction[[w]] <- prediction[[w]] + prob
    }
    infected_number[p, ] <- number_infected_runs[[p]]
    infected_area[p, ] <- area_infected_runs[[p]]
    rates <- do.call(rbind, spread_rate_runs[[p]])
    if (!is.null(rates)) {
      west_rates[p, ] <- rates[, 4]
      east_rates[p, ] <- rates[, 3]
      south_rates[p, ] <- rates[, 2]
      north_rates[p, ] <- rates[, 1]
    }else {
      west_rates[p, ] <- 0
      east_rates[p, ] <- 0
      south_rates[p, ] <- 0
      north_rates[p, ] <- 0
    }

    if (config$use_quarantine &
        length(quarantine_escape_runs[[p]]) ==
        length(probability_runs[[p]])) {
      escape_probability <- escape_probability + quarantine_escape_runs[[p]]
      quarantine_escapes[p, ] <- quarantine_escape_runs[[p]]
      quarantine_escape_distances <- quarantine_escape_distance_runs[[p]]
      quarantine_escape_directions <- quarantine_escape_directions_runs[[p]]
    }
  }

  probability <- prediction
  for(w in seq_len(length(prediction))) {
    probability[[w]] <- (prediction[[w]] / (length(probability_runs))) * 100
  }

  infected_areas <-
    round(sapply(infected_area, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                              "Stand dev" = sd(x))), digits = 0)
  number_infecteds <-
    round(sapply(infected_number, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                                "Stand dev" = sd(x))),
          digits = 0)
  west_rate <-
    round(sapply(west_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                           "Stand dev" = sd(x))), digits = 0)
  east_rate <-
    round(sapply(east_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                           "Stand dev" = sd(x))), digits = 0)
  south_rate <-
    round(sapply(south_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                            "Stand dev" = sd(x))), digits = 0)
  north_rate <-
    round(sapply(north_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                            "Stand dev" = sd(x))), digits = 0)

  west_rate[is.na(west_rate)] <- 0
  east_rate[is.na(east_rate)] <- 0
  south_rate[is.na(south_rate)] <- 0
  north_rate[is.na(north_rate)] <- 0

  which_median <- function(x) raster::which.min(abs(x - median(x)))

  median_run_index <- which_median(infected_number[[ncol(infected_number)]])
  min_run_index <- which.min(infected_number[[ncol(infected_number)]])
  max_run_index <- which.max(infected_number[[ncol(infected_number)]])

  single_run <- single_runs[[median_run_index]]
  susceptible_run <- susceptible_runs[[median_run_index]]
  exposed_run <- exposed_runs[[median_run_index]]

  min_run <- single_runs[[min_run_index]]
  max_run <- single_runs[[max_run_index]]

  for (q in seq_len(length(single_runs[[1]]))) {
    for (j in seq_len(length(single_runs))) {
      if (j == 1) {
        raster_stacks <- list(single_runs[[j]][[q]])
      } else {
        raster_stacks[[j]] <- single_runs[[j]][[q]]
      }
    }
    raster_stacks2 <- do.call(cbind, raster_stacks)
    raster_stacks2 <- array(raster_stacks2, dim=c(dim(raster_stacks[[1]]), length(raster_stacks)))
    simulation_mean <- round(apply(raster_stacks2, c(1, 2), mean, na.rm = TRUE), digits = 0)
    # simulation_min <- apply(raster_stacks2, c(1, 2), min, na.rm = TRUE)
    # simulation_max <- apply(raster_stacks2, c(1, 2), max, na.rm = TRUE)
    simulation_sd <- apply(raster_stacks2, c(1, 2), sd, na.rm = TRUE)

    simulation_stack <-
      terra::rast(nrow = config$num_rows, ncol = config$num_cols,
                  xmin = config$xmin, xmax = config$xmax,
                  ymin = config$ymin, ymax = config$ymax, crs = config$crs)
    terra::values(simulation_stack) <- simulation_mean
    names(simulation_stack) <- 'mean'
    simulation_stack$probability <- probability[[q]]
    simulation_stack$max <- max_run[[q]]
    simulation_stack$min <- min_run[[q]]
    simulation_stack$standard_deviation <- simulation_sd
    simulation_stack$median <- single_run[[q]]

    if (q == 1) {
      simulation_stacks <- list(simulation_stack)
    } else {
      simulation_stacks[[q]] <- simulation_stack
    }
  }

  run2$status <- "WRITING DATA"
  httr::PUT(url = paste(api_url, "run_write/", run_id, "/", sep = ""),
            body = run2, encode = "json")

  if (run_id == session$default_run) {
    max_value <-
      round(sapply(max_values, function(x)
        c("Mean" = mean(x, na.rm = TRUE), "Stand dev" = sd(x))), digits = 0)
    max_value_out <- max_value[1, ncol(max_value)]
    session$max_value <- max_value_out
    httr::PUT(url = paste(api_url, "session_write/", session_id, "/", sep = ""),
              body = session, encode = "json")
  }

  single_run_out <- single_run[[1]]
  susceptible_run_out <- susceptible_run[[1]]
  exposed_run_out <- exposed_run[[1]]

  # stuff <-
    for(q in 1:length(years)) {
    # foreach::foreach(q = seq_len(length(years)),
    #                  .packages = c("terra", "geojsonio", "httr")) %do% {
    if (q == 1) {
      number_infected <- infected_number[median_run_index, q]
      area_infected <- infected_area[median_run_index, q]
      if (is.nan(west_rates[median_run_index, q])) {
        west_rate_r <- 0
      } else {
        west_rate_r <- round(west_rates[median_run_index, q], digits = 0)
      }
      if (is.nan(east_rates[median_run_index, q])) {
        east_rate_r <- 0
      } else {
        east_rate_r <- round(east_rates[median_run_index, q], digits = 0)
      }
      if (is.nan(south_rates[median_run_index, q])) {
        south_rate_r <- 0
      } else {
        south_rate_r <- round(south_rates[median_run_index, q], digits = 0)
      }
      if (is.nan(north_rates[median_run_index, q])) {
        north_rate_r <- 0
      } else {
        north_rate_r <- round(north_rates[median_run_index, q], digits = 0)
      }
    } else if (q > 1) {
      number_infected <- number_infecteds[1, q]
      area_infected <- infected_areas[1, q]
      if (is.nan(west_rate[1, q])) {
        west_rate_r <- 0
      } else {
        west_rate_r <- west_rate[1, q]
      }
      if (is.nan(east_rate[1, q])) {
        east_rate_r <- 0
      } else {
        east_rate_r <- east_rate[1, q]
      }
      if (is.nan(south_rate[1, q])) {
        south_rate_r <- 0
      } else {
        south_rate_r <- south_rate[1, q]
      }
      if (is.nan(north_rate[1, q])) {
        north_rate_r <- 0
      } else {
        north_rate_r <- north_rate[1, q]
      }
    }

    year <- years[q]

    single_map <- simulation_stacks[[q]]
    infec <- single_map$mean
    expos <- single_map$mean
    values(expos) <- 0
    values(infec) <- 0
    values(infec) <- config$infected
    for (m in seq_len(length(config$exposed))) {
      values(expos) <- as.matrix(expos, wide = TRUE) + config$exposed[[m]]
    }
    values(infec) <- values(infec) + values(expos)
    single_map[single_map <= 0] <- NA
    if (terra::global(single_map, fun = "sum", na.rm = TRUE) == 0 ||
        is.na(terra::global(single_map, fun = "sum", na.rm = TRUE))) {
      single_map[infec > 0] <- 0
    }
    reso <- config$ew_res

    # if (terra::ncell(single_map) <= 100000) {
    #   single_map_s <- st_as_stars(single_map)
    # } else if (terra::ncell(single_map) > 100000) {
    #   single_map_b <- terra::aggregate(single_map, fact = 2, fun = "mean")
    #   single_map_s <- st_as_stars(single_map_b)
    #   reso <- reso * 2
    # }
    single_map_s <- st_as_stars(single_map)
    # st_crs(single_map_s) <- 3857
    single_map_s <- st_transform(single_map_s, 4326)
    single_map_p <- st_as_sf(single_map_s)
    # single_map_p <- terra::as.polygons(single_map, dissolve = FALSE,
                                       # values = TRUE, trunc = TRUE)
    storage.mode(single_map_p$median) <- "integer"
    storage.mode(single_map_p$max) <- "integer"
    storage.mode(single_map_p$probability) <- "integer"
    storage.mode(single_map_p$min) <- "integer"
    storage.mode(single_map_p$mean) <- "integer"
    storage.mode(single_map_p$standard_deviation) <- "integer"
    if (length(single_map_p$min[is.na(single_map_p$min)])) {
      single_map_p$min[is.na(single_map_p$min)] <- 0
    }
    if (length(single_map_p$mean[is.na(single_map_p$mean)])) {
      single_map_p$mean[is.na(single_map_p$mean)] <- 0
    }
    if (length(single_map_p$max[is.na(single_map_p$max)])) {
      single_map_p$max[is.na(single_map_p$max)] <- 0
    }
    if (length(single_map_p$median[is.na(single_map_p$median)])) {
      single_map_p$median[is.na(single_map_p$median)] <- 0
    }
    if (length(single_map_p$probability[is.na(single_map_p$probability)])) {
      single_map_p$probability[is.na(single_map_p$probability)] <- 0
    }
    if (length(single_map_p$standard_deviation[is.na(single_map_p$standard_deviation)])) {
      single_map_p$standard_deviation[is.na(single_map_p$standard_deviation)] <- 0
    }
    # single_map_p <- terra::project(single_map_p, "epsg:4326")
    # single_map_p <- as(single_map_p, "sf")
    # single_map_p <- sf::st_as_sf(single_map_p)

    # single_map_p <- as(single_map_p, "Spatial")

    if (reso <= 10) {
      single_map_p <-
        # geojsonio::geojson_list(single_map_p, precision = 6, geometry = "polygon")
        sf_geojson(single_map_p, digits = 6, simplify = FALSE)
    } else if (reso > 10 && reso <= 100) {
      single_map_p <-
        # geojsonio::geojson_list(single_map_p, precision = 5, geometry = "polygon")
        sf_geojson(single_map_p, digits = 5, simplify = FALSE)
    } else if (reso > 100 && reso <= 500) {
      single_map_p <-
        # geojsonio::geojson_list(single_map_p, precision = 4, geometry = "polygon")
        sf_geojson(single_map_p, digits = 4, simplify = FALSE)
    } else if (reso > 500) {
      single_map_p <-
        # geojsonio::geojson_list(single_map_p, precision = 3, geometry = "polygon")
        sf_geojson(single_map_p, digits = 3, simplify = FALSE)
      # single_map_p3  <- geojsonio::geojson_list(single_map_p, precision = 3, geometry = "polygon")
    }

    # class(single_map_p) <- "list"

    outs <- list()
    outs$run <- run_id
    outs$number_infected <- number_infected
    outs$infected_area <- round(area_infected, digits = 2)
    outs$year <- year
    outs$median_spread_map <- NULL
    outs$probability_map <- NULL
    outs$max_spread_map <- NULL
    outs$min_spread_map <- NULL
    outs$standard_deviation_map <- NULL
    outs$mean_spread_map <- NULL
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
    post_code <-
      httr::POST(url = paste(api_url, "output/", sep = ""),
                 body = outs, encode = "json")

    if (post_code$status_code == 201) {
      run2$status <- "SUCCESS"
      output_id <- httr::content(post_code)
      output_id <- output_id$pk
      s <- httr::PUT(url = paste(api_url, "output_spread_map/", output_id, "/", sep = ""),
                     body = list(median_spread_map = single_map_p))
    } else {
      run2$status <- "FAILED"
    }

    print(q)
    stuff <- run2$status
  }

  # run2$status <- stuff[[1]]
  run2$status <- stuff

  if (run2$status == "SUCCESS") {
    httr::PUT(url = paste(api_url, "run_write/", run_id, "/", sep = ""),
              body = run2, encode = "json")
  } else {
    httr::PUT(url = paste(api_url, "run_write/", run_id, "/", sep = ""),
              body = run2, encode = "json")
  }

  if (run_collection$default == TRUE ||
      is.null(run$steering_year) ||
      run$steering_year == lubridate::year(config$end_date)) {
    if (run2$status == "SUCCESS") {
      run_collection$status <- "SUCCESS"
      httr::PUT(url = paste(api_url, "run_collection_write/",
                            run_collection_id, "/", sep = ""),
                body = run_collection, encode = "json")
    } else {
      run_collection$status <- "FAILED"
      httr::PUT(url = paste(api_url, "run_collection_write/",
                            run_collection_id, "/", sep = ""),
                body = run_collection, encode = "multipart")
    }
  }

  config$infected <- single_run_out
  config$susceptible <- susceptible_run_out
  config$exposed <- exposed_run_out

  ## write out config to run api
  save(config, file = "current_run.RData")
  run_cs <- httr::upload_file("current_run.RData")
  httr::PUT(url = paste(api_url, "run_r_data/", run_id, "/", sep = ""),
            body = list(r_data = run_cs))

  status <- run2$status
  }
