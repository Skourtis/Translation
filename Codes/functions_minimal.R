pacman::p_load("biglm",
               "readxl",
               "rmarkdown",
               "tarchetypes",
               "targets",
               "tidyverse")

#' @title Plot ozone from the preprocessed air quality data.
#' @description Plot a histogram of ozone concentration.
#' @return A ggplot histogram showing ozone content.
#' @param data Data frame, preprocessed air quality dataset.
#' @examples
#' library(ggplot2)
#' library(tidyverse)
#' data <- airquality %>%
#'   mutate(Ozone = replace_na(Ozone, mean(Ozone, na.rm = TRUE)))
#' create_plot(data)
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone), bins = 12) +
    theme_gray(24)+
        ggtitle("Hist_1")
}

####


loading_CCLE_prot <- function(x){
    CCLE_proteins <- openxlsx::read.xlsx(x, sheet =2) %>%
        .[,!str_detect(colnames(.),"Peptides|Column")]
    colnames(CCLE_proteins) <- str_remove_all(colnames(CCLE_proteins),"_TenPx..")
    rownames(CCLE_proteins) <- NULL
    CCLE_proteins <- CCLE_proteins %>%
        column_to_rownames(var = "Uniprot_Acc")
    CCLE_proteins <- CCLE_proteins[,-c(1:5)]
    colnames(CCLE_proteins) <- str_match(colnames(CCLE_proteins),"^([:graph:]*?)_")[,2]
    CCLE_proteins <- CCLE_proteins[(CCLE_proteins %>% is.na() %>% rowSums()) < ncol(CCLE_proteins)*0.75, # removing lines with more than 75% NA
                                   which(!duplicated(colnames(CCLE_proteins)))] ##Removing cell_lines with the same name
    
}
Z_transf <- function(x){
    x <- as.matrix(x)
    (x - mean(x, na.rm = T))/(sd(x, na.rm = T))
}