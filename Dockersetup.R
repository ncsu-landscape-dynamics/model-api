library(containerit)

dockerfile <- dockerfile("deploy_staging/schedule.R", copy = "script_dir", soft = TRUE)
write(dockerfile, file = "Dockerfile")
