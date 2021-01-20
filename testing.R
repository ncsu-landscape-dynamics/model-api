
json_output<- httr::GET(paste("https://pops-model.org/api/output/", 158, "/?format=json", sep = ""))
output <- httr::content(json_output)
post_code <- httr::POST(url = "https://pops-model.org/api/output/", body = output, encode = "json")
code <- httr::content(post_code)
output$run <- 1488

## run with management for parsing
case_study_id = "1"
session_id = "1"
run_collection_id = "5"
run_id = "2"

library(aws.s3)
bucket = 'pops-production'

t <- config$infected_file
save_file = aws.s3::save_object(object = t, bucket = bucket, file = t, check_region = FALSE)
result = raster::stack(t)

save_object("s", object = t, bucket = 'pops-production')
s <- raster("H:/Shared drives/Data/Raster/Regional/SLF_100m/tree_of_heaven_100m.tif")
