source_file_name <- here::here("Codes",
                               "functions.R")
source(source_file_name)
library(testthat)


test_that("addition results is correct",{
    expect_true(is.data.frame(addition(5,6) %>% as.data.frame()))
    testthat::expect_error(addition("ab","cd"))
})

