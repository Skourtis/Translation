library(targets)
library(pacman)
setwd("/nfs/users2/ssdelci/skourtis/Translation/")
library(here)
library(future)

tar_make_future(
  names = NULL,
  reporter = "verbose",
  workers = future::availableCores(),
  callr_function = callr::r,
  callr_arguments = list()
)
