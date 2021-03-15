source_file_name <- here::here("Codes",
                               "functions_minimal.R")
source(source_file_name)
library(testthat)


test_that("Intra_omic Correlation",{ #checks that the columns and row names from the intra_omic correlation are identical
    Omic_test <- data.frame(HEL = runif(3,0,1),MDAMB231 = runif(3,0,1), T47D = runif(3,0,1), 
                            row.names = c("RANBP6","UNC119","HES4")) %>% 
        Intra_omic_Feature_correlation(Genes = head(rownames(.),2))
    
    expect_identical(colnames(Omic_test),rownames(Omic_test))
    expect_true(between(Omic_test,-1,1) %>% na.omit() %>% all) #checks that all numbers produced are between -1 1
    #testthat::expect_error(addition("ab","cd"))
})

test_that("Heatmap_Output",{ #Simply checks that the output from this test is of Class heatmap
    Data_df <- data.frame(RANBP6 = c(1,0.16,-1),UNC119 = c(0.16,1,-0.25), HES4 = c(-1,0.25,1), row.names = c("RANBP6","UNC119","HES4"))
   Main_top_annotation <- data.frame(Genes = c("UNC119", "HES4"), Pathways = c("b","b"))
   expect_true(Convert_df_to_heatmap(Data_df = Data_df, 
                         Main_top_annotation = Main_top_annotation,
               Subset_row = "RANBP6")%>% class(.) %>%.[[1]] == "Heatmap")
    
})

test_that("Matrix_Alignment",{ #Checks that the two matrices have the same dimensions
    Matrix_1 <- matrix(1:20,nrow = 4,ncol=5, dimnames = list(letters[1:4],LETTERS[1:5]))
    Matrix_2 <- matrix(1:30,nrow = 5,ncol=6, dimnames = list(letters[2:6],LETTERS[3:8]))
    colnames <- Aligning_two_matrices(Matrix_1,Matrix_2)%>% map(.x =.,~colnames(x = .x))
    rownames <- Aligning_two_matrices(Matrix_1,Matrix_2)%>% map(.x =.,~rownames(x = .x))
    expect_true(c(colnames[[1]]==colnames[[2]],rownames[[1]] == rownames[[2]]) %>% all)
})
