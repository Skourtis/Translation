#Downloading from Git template
#Cachem install error bc 'make' was not installed, Rtools installed
#Manually created .Renv and paste the path from Rtools webpage.
renv::restore()
install.packages('devtools')

options(repos = getOption("repos")["CRAN"])
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)
pacman::p_load(piggyback, renv, here, tidyverse)
pb_download()
