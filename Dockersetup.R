# devtools::install_github("ChrisJones687/geojsonio")
library(containerit)

dockerfile <- dockerfile("deploy_model/schedule.R", copy = "script_dir", soft = TRUE)
write(dockerfile, file = "Dockerfile")
