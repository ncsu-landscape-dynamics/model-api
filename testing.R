case_study_id = "3"
session_id = "81"
run_collection_id = "545"
run_id = "796"

case_study_id = "5"
session_id = "92"
run_collection_id = "561"
run_id = "817"


set.seed(33)

random_seeds <- sample(1:10000000, 10)
random_seeds

httr::DELETE("https://pops-model.org/api/output/941")


case_study_id = "2"
session_id = "46"
run_collection_id = "350"
run_id = "526"
