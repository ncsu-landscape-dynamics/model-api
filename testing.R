case_study_id = "3"
session_id = "84"
run_collection_id = "522"
run_id = "773"

case_study_id = "3"
session_id = "81"
run_collection_id = "522"
run_id = "753"


set.seed(33)

random_seeds <- sample(1:10000000, 10)
random_seeds

httr::DELETE("https://pops-model.org/api/output/941")


case_study_id = "2"
session_id = "46"
run_collection_id = "350"
run_id = "526"


case_study_id = "2"
run_id = "1488"
run_collection_id = "925"
session_id = "130"


json_output<- httr::GET(paste("https://pops-model.org/api/output/", 158, "/?format=json", sep = ""))
output <- httr::content(json_output)
post_code <- httr::POST(url = "https://pops-model.org/api/output/", body = output, encode = "json")
code <- httr::content(post_code)
output$run <- 1488
