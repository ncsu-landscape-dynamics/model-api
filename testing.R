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
