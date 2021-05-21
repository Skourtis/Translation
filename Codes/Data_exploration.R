### Data exploration
library(tidyverse)
library(tidymodels)
library(modeldata)
library(furrr)
library(progressr)
library(ComplexHeatmap)
data(meats)
targets::tar_load(CCLE_proteins)
targets::tar_load(PTR_CCLE)
targets::tar_load(uniprot_factors)
Proteins <-  inner_join(CCLE_proteins %>% t() %>%  
                            as.data.frame() %>% 
                            #dplyr::select(Uniprot) %>%  
                            set_names(.,glue::glue("{colnames(.)}_Uniprot")) %>% 
                            rownames_to_column("Cell_line") ,
                        PTR_CCLE %>% t() %>% as.data.frame %>%
                            dplyr::select(any_of(uniprot_factors %>% 
                                                     subset(type = "eukaryotic_translation") %>% 
                                                     pull("Gene names") %>% 
                                                     str_split( " ") %>% 
                                                     unlist %>% str_subset("EIF|ETF|EEF"))) %>% 
                            
                            rownames_to_column("Cell_line")) 
Proteins <-  select(Proteins,which((Proteins %>% 
                                        is.na() %>% 
                                        colSums())==0)) #%>% column_to_rownames("Cell_line") #%>% as.matrix()# <nrow(IDH1)/10))
list_of_Uniprot_models <- map(.x = colnames(Proteins) %>% 
                   str_subset("Uniprot"),
               ~dplyr::select(Proteins, contains(.x)|!contains("Uniprot"))) %>% set_names(colnames(Proteins) %>% 
                                                                                              str_subset("Uniprot"))

Linear_model_function <- function(Protein_list){

office_split <- initial_split(Protein_list)# , strata = season)
office_train <- training(office_split)
office_test <- testing(office_split)
# Importance = data.frame(Variable = NULL,
#                         Importance = NULL,
#                         Sign = NULL,
#                         Uniprot = NULL)
# Importance = data.frame(.metric = NULL,
#                         .estimator = NULL,
#                         .estimate = NULL,
#                         .config = NULL,
#                         Uniprot = NULL)
    

rec <- recipe(Protein_list) %>% 
    update_role(contains("Uniprot"), 
                new_role = "outcome") %>% 
    update_role(!contains("Uniprot"), 
                new_role = "predictor") %>% 
    update_role(contains("Cell_line"), 
                new_role = "ID") %>% 
    step_normalize(all_predictors()) #%>%
#step_nzv(all_predictors())


office_prep <- rec %>%
    prep(strings_as_factors = FALSE)

lasso_spec <- linear_reg(penalty = 0.1, mixture = 0) %>%
    set_engine("glmnet")
wf <- workflow() %>%
    add_recipe(rec)
lasso_fit <- wf %>%
    add_model(lasso_spec) %>%
    fit(data = office_train)
lasso_fit %>%
    pull_workflow_fit() %>%
    tidy()
set.seed(1234)
office_boot <- bootstraps(office_train)

tune_spec <- linear_reg(penalty = tune(), mixture = 0) %>%
    set_engine("glmnet")

lambda_grid <- grid_regular(penalty(), levels = 10)
set.seed(2020)
lasso_grid <- tune_grid(
    wf %>% add_model(tune_spec),
    resamples = office_boot,
    grid = lambda_grid
)
lasso_grid %>%
    collect_metrics()

lasso_grid %>%
    collect_metrics() %>%
    ggplot(aes(mixture, mean, color = .metric)) +
    geom_errorbar(aes(
        ymin = mean - std_err,
        ymax = mean + std_err
    ),
    alpha = 0.5
    ) +
    geom_line(size = 1.5) +
    facet_wrap(~.metric, scales = "free", nrow = 2) +
    scale_x_log10() +
    theme(legend.position = "none")
lowest_rmse <- lasso_grid %>%
    select_best("rmse")
final_lasso <- finalize_workflow(
    wf %>% add_model(tune_spec),
    lowest_rmse
)
#library(vip)
importance <- final_lasso %>%
    fit(office_train) %>%
    pull_workflow_fit() %>%
    vip::vi(lambda = lowest_rmse$penalty) %>%
    mutate(
        Importance = abs(Importance),
        Variable = fct_reorder(Variable, Importance)
    ) 
# %>% #subset(Importance>0.01) %>% 
#     ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
#     geom_col() +
#     scale_x_continuous(expand = c(0, 0)) +
#     labs(y = NULL)
final_fitted <- last_fit(
    final_lasso,
    office_split
) %>%
    collect_metrics()
list(importance = importance,
     final_fitted = final_fitted)

}
options(future.globals.maxSize= 891289600)
Linear_model_function_safe <- safely(Linear_model_function)
library(future)
plan(multisession, workers = 5)
Linear_model_f <-  furrr::future_map(list_of_Uniprot_models[1:20], 
                                                  Linear_model_function_safe,.options = furrr_options(seed = 1234),
                                     .progress = TRUE)
future::plan(sequential)

Highly_predicted <- Linear_model_f %>% imap_dbl(.x = . ,~.x[["result"]][['final_fitted']] %>% 
                                subset(.metric == "rsq") %>% 
                                pull(.estimate)) %>% 
    sort() %>% 
    tail(20)

Linear_model_f %>% map(.x = . ,~.x[["result"]][["importance"]]) %>% .[names(Highly_predicted)] %>% 
    imap(.x =.,~ subset(.x,str_detect(Variable,"EIF3")) %>% ggplot(aes(x = Importance,fill = Sign,y = Variable))+geom_col() +
             ggtitle(.y))
Linear_model_f %>% imap_dfr(.x = . ,~.x[["result"]][["importance"]] %>% 
                                mutate(Uniprot = .y,
                                       Importance = if_else(Sign == "POS",Importance,Importance*(-1)))) %>% 
    subset(str_detect(Uniprot,paste(names(Highly_predicted),collapse = "|" ))) %>% 
    pivot_wider(id_cols = Variable,
                names_from = Uniprot,
                values_from = Importance) %>% subset(str_detect(Variable,"EIF3"))%>% column_to_rownames("Variable")  %>% t() %>%  cor() %>% 
    
    Convert_df_to_heatmap(Data_df = ., Heatmap_name = "eTFs Proteome against themselves")
########TCGA Example from https://www.jkangpathology.com/post/tidymodel-and-glmnet/

### Data exploration
pacman::p_load(targets,tidyverse,ComplexHeatmap)

tar_load(CCLE_prot_Intra_omics_corr_matrix)
tar_load(uniprot_factors)
Convert_df_to_heatmap(Data_df = CCLE_prot_Intra_omics_corr_matrix, Heatmap_name = "eTFs Proteome against themselves",
                      Subset_col =uniprot_factors %>% subset(Type == "eukaryotic_translation") %>%
                          pull(`Gene names`) %>% str_split(" ")  %>% unlist()%>% str_subset("^(EIF3)"), 
                      Subset_row = uniprot_factors %>% subset(Type == "eukaryotic_translation") %>%
                          pull(`Gene names`) %>% str_split(" ")  %>% unlist()%>% str_subset("^(EIF3)"),
                      show_names = "both")
##### PTR From Kuster####
library("dgof")
testing <- mRNA_Prot_Cor_CCLE %>% 
    group_by(lineage) %>% 
    group_split() %>%
    set_names(.,purrr::map_chr(.x = .,~ .x %>% pull('lineage') %>% unique()))%$%
    purrr::map2(.x = map(.x = .,~.x %>% subset(pathway == "OTHER") %>% pull('Gene_mRNA_prot')),
                .y = map(.x = .,~.x %>% subset(pathway == "OXPHOS") %>% pull('Gene_mRNA_prot')),
                ~ks.test(.x,.y))
CCLE_Pathway_Cor_diff  <-  
    purrr::map(.x = unique(KEGG_genes$pathway),~mRNA_Prot_Cor_CCLE %>%
                   arrange(-Gene_mRNA_prot) %>% ungroup()%>%
                   mutate(pathway_KEGG = if_else(GeneName %in% 
                                                     (KEGG_genes %>% subset(pathway == .x) %>% pull(ID)),
                                                 .x,
                                                 "OTHER"),
                          mean_dif = Gene_mRNA_prot-mean(Gene_mRNA_prot )) %>%
                   subset(pathway_KEGG == .x  ) %>% pull(mean_dif,GeneName))%>% 
    set_names(unique(KEGG_genes$pathway))%>% .[purrr::map_lgl(.x =., ~length(.x)>10)]

Kuster_Pathway_Cor_diff <- purrr::map(.x = unique(KEGG_genes$pathway),~mRNA_Prot_Cor_Kuster %>%
                                          arrange(-Gene_mRNA_prot) %>% ungroup()%>%
                                          mutate(pathway_KEGG = if_else(GeneName %in% 
                                                                            (KEGG_genes %>% subset(pathway == .x) %>% pull(ID)),
                                                                        .x,
                                                                        "OTHER"),
                                                 mean_dif = Gene_mRNA_prot-mean(Gene_mRNA_prot )) %>%
                                          subset(pathway_KEGG == .x  ) %>% pull(mean_dif,GeneName))%>% 
    set_names(unique(KEGG_genes$pathway)) %>% .[purrr::map_lgl(.x =., ~length(.x)>10)]

healthy_cancer_corr_pathways<-purrr::map_dbl(.x = intersect(names(Kuster_Pathway_Cor_diff),
                          names(CCLE_Pathway_Cor_diff)),
                          ~ks.test(Kuster_Pathway_Cor_diff[[.x]],
                                   CCLE_Pathway_Cor_diff[[.x]]) %>% 
                            pluck("p.value") )%>% set_names(intersect(names(Kuster_Pathway_Cor_diff),
                                                                      names(CCLE_Pathway_Cor_diff)))

purrr::map(c("KEGG_genes","PTR_Kuster","KEGG_pathways"), tar_load)
targets::tar_load(KEGG_genes)
targets::tar_load(KEGG_pathways)
tar_load(PTR_Kuster)
#### testing which pathways are different overall #### needs to be doen on unnormalised data
PTR_different_pathways <- purrr::map(.x = unique(pathway_PTR_Kuster$pathway),
           ~ks.test(pathway_PTR_Kuster %>% subset(pathway == .x) %>% pull(value),
                    pathway_PTR_Kuster %>% subset(pathway != .x) %>% pull(value)) %>%
               pluck("p.value")) %>% set_names(unique(pathway_PTR_Kuster$pathway))
Tissue_pathway_PTR_test <- pathway_PTR_CCLE %>% 
    group_split(pathway)  %>%
    purrr::map(.x = .,Calc_Tissu_PTR ) %>% 
    set_names(unique(pathway_PTR_Kuster$pathway))
Calc_Tissu_PTR <- function(df_pathway){
    #df_pathway = pathway_PTR_Kuster %>% 
    #    group_split(pathway) %>% .[[1]]
    map(.x = unique(df_pathway$Tissue),
        ~ks.test(df_pathway %>% subset(Tissue == .x) %>% pull(value),
            df_pathway %>% subset(Tissue != .x) %>% pull(value)) %>%
            pluck("p.value")) %>% set_names(unique(df_pathway$Tissue))
}
####
Get_PC_complexes <- function(Identifiers,Class_Name, Data_Matrix){
    #Data_Matrix =  CCLE_proteins
    #Identifiers = Complex_complete_units[[1]] %>% unlist()
    
    pca_fit_original <- Data_Matrix%>% 
        subset(.,rownames(.) %in% unlist(Identifiers)) %>%
        subset(is.na(.) %>% rowSums<(ncol(.)/4))%>% as.matrix() %>%  
        impute::impute.knn() %>% .[["data"]] %>% t()%>% prcomp()
    
    PC_for_60<-pca_fit_original %>%
        broom::tidy(matrix = "eigenvalues") %>%
        subset(cumulative>0.6) %>% head(1) %>% pull(PC)
    pca_fit_original$x[,1:PC_for_60] %>% as.data.frame() %>% set_names(.,paste0(colnames(.),Class_Name))
}

BiocManager::install("OmnipathR")
library(OmnipathR)
