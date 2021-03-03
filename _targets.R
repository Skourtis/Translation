library(targets)
library(tarchetypes)
source(here::here("Codes","functions_minimal.R"))
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("biglm", "tidyverse"))
#in the console run targets::tar_visnetwork()
#and targets::tar_make()
list(
    tar_target(
        raw_data_file,
        here::here("Datasets","Raw","raw_data.csv"),
        format = "file"
    ),
    tar_target(
        raw_data,
        read_csv(raw_data_file, col_types = cols())
    ),
    tar_target(
        data,
        raw_data %>%
            mutate(Ozone = replace_na(Ozone, mean(Ozone, na.rm = TRUE)))
    ),
    tar_target(hist, create_plot(data)),
    tar_render(report, here::here("Output","report.Rmd")),
    tar_target(fit, biglm(Ozone ~ Wind + Temp, data)),
    tar_target(dataset, data.frame(x = letters))
 )
