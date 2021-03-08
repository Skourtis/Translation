#' @title Proteomics_imputation
#' @description Import and clean proteomic CCLE Data
#' @return An imputted matrix in the same format as input
#' @param Takes in a wide matrix with samples in columns, proteins in row names. (Requires column and row names)
#' @param N_condition, N_rep
#' @examples set.seed(123) x <- matrix(sample(17:35,9000,replace = T),nrow = 100, ncol = 9, dimnames = list(paste0("Uniprot_",1:100),paste0("sample_",1:9))) 
#' @examples x[sample(1:900,300)] <- NA

Proteomics_imputation <- function(Data, N_condition, N_rep){
    # set.seed(123)
    # x <- matrix(sample(17:35,9000,replace = T),nrow = 100, ncol = 9, dimnames = list(paste0("Uniprot_",1:100),paste0("sample_",1:9))) 
    # x[sample(1:900,300)] <- NA
    # Data <- x
    # N_condition = 3
    # N_rep <- 3
    warning("Replicates should be next to each other in the Matrix")
    colnames_data <- data.frame(Raw_names = colnames(Data), #Creates_new column names that cen be used to find replicates
                                Processed = paste0(rep(LETTERS[1:N_condition],each = N_rep),"_",rep(1:N_condition,N_rep)))
    
   
    
    deciles <- quantile(x %>% as.vector() %>% na.omit(), #### finding the bottom 2 deciles of the distribution 
                        prob = seq(0, 1, length = 11), type = 5) %>%
        .[1:2]
    Data <- as.data.frame(x) %>%
        rownames_to_column("Uniprot") %>%
        pivot_longer(-Uniprot, names_to = "Raw_names", values_to = "Abundance") %>%
        left_join(colnames_data) %>% 
        mutate(Processed = str_remove(Processed, "_.")) %>%
        group_by(Uniprot,Processed) %>%
        mutate(N_NA = sum(is.na(Abundance)),
               Expression = case_when(is.na(Abundance) ~ "Low", # decides on whether an Abundance is Low or High
                                      Abundance< deciles[3] ~ "Low",
                                      Abundance>=deciles[3] ~ "Mid/High")) %>%
        mutate(Group_Abundance = if_else(length(unique(Expression))<2,"Single","Mixed")) %>% #decides on whether replicates show the same behaviour
        mutate(Imputted_Abundance = case_when(
            is.na(Abundance) & N_NA >1 & Expression == "Low" & Group_Abundance == "Single" ~rnorm(1,mean = mean(deciles), sd = 1),
            TRUE ~ as.numeric(Abundance)) #If the Original Abundance is low, The replicates show the same behaviour and the expression is low , Imputes with ditribution
        ) %>% ungroup()
    Data %>% 
        dplyr::select(Uniprot, Raw_names,Imputted_Abundance) %>%
        pivot_wider(names_from = "Raw_names",values_from = "Imputted_Abundance") %>%
        column_to_rownames("Uniprot") %>% as.matrix()
    
}
