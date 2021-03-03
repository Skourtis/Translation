#Settting up directory
##Settting up directory
options(repos = getOption("repos")["CRAN"])
install.packages("PTXQC")
install.packages("pacman")
pacman::p_load(piggyback, renv, here, tidyverse, targets,
               visNetwork)
install.packages('piggyback')
install.packages('renv')
install.packages('here')
install.packages('tidyverse')
install.packages('targets')
install.packages('visNetwork')
#testthat::use_test()

## Created a first release directly on Github
#pb_new_release("Skourtis/Project_Template")
piggyback::pb_track(c("Datasets/Raw/*.txt",
                      "Datasets/Raw/*.dat",
                      "Datasets/Raw/*.zip",
                      "Datasets/Raw/*.csv",
                      "Datasets/Raw/*.RData"))

piggyback::pb_track() %>%
    pb_upload(repo = "Skourtis/Project_Template")

##end
renv::snapshot()


