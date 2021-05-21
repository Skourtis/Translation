#targets::tar_load(eTFs_PTR_corr)
#targets::tar_load(CCLE_RNA_seq)
humanreactome_df <- read_csv2(here::here("Datasets","Raw","Reactome_pathways.csv"))[,-1]
eTFs_PTR_corr_wide <-  eTFs_PTR_corr %>% 
    pivot_wider(names_from = Uniprot_PTR, values_from = Correlation) %>%
    left_join(tar_read(uniprot_factors) %>% dplyr::select(Entry,`Entry name`) %>% distinct, by = c("Factor" =  "Entry"))
pca_fit <- eTFs_PTR_corr_wide %>% 
    dplyr::select(where(is.numeric))%>%
    as.matrix() %>%
    prcomp() 

# Families <- getBM(attributes=c("external_gene_name","superfamily"),
#                   filters = 'external_gene_name',
#                   values = unique(eIF_correlations_wide$ID) , mart =ensembl) %>%
#     set_names("ID","family") %>%
#     subset(!duplicated(ID)) 
#library(bioMart)
library(cowplot)
pca_fit %>%
    broom::augment(eTFs_PTR_corr_wide) %>% # add original dataset back in
    distinct()%>%
    #left_join(Families,by = "ID" ) %>%
    #mutate(family = if_else(str_detect(family,"SS"),family,"NA"))%>%
    ggplot(aes(.fittedPC1, .fittedPC2)) + 
    geom_point(size = 1.5, show.legend = FALSE) + theme(legend.position = "none") +
    ggrepel::geom_text_repel(aes(label = `Entry name`),show.legend = FALSE) + theme(legend.position = "none") +
    theme_half_open(12) + background_grid()
tar_load(PTR_CCLE)
PCA_dataset <- pca_fit$rotation[1:2]
PTR_CCLE_na_rm <- PTR_CCLE %>% subset(is.na(.) %>% rowSums<(ncol(.)/4)) %>% #is.na() %>% table() 
 t() %>% as.data.frame() %>%  rownames_to_column("Cell_line")
tar_read(PTR_Kuster)  %>% subset(is.na(.) %>% rowSums<(ncol(.)/4)) %>% #is.na() %>% table() 
  t() %>% as.data.frame() %>%  rownames_to_column("Cell_line")

tar_load(Protein_Kuster)

Combined_for_multivariate <- inner_join(  tar_read(PTR_Kuster)  %>%
                                            subset(is.na(.) %>% rowSums<(ncol(.)/4)) %>% #is.na() %>% table() 
                                            t() %>% as.data.frame() %>%
                                            dplyr::select(-matches("[a-z]", ignore.case = FALSE)) %>%
                                            rownames_to_column("Cell_line")
                                            ,
                                        CCLE_Complex_PCs %>% as.data.frame() %>%  rownames_to_column("Cell_line")) %>% 
    column_to_rownames("Cell_line") %>% as.matrix() %>% impute::impute.knn() %>% .[["data"]] %>% as.data.frame() %>% 
    set_names(.,str_remove_all(colnames(.),"-|\\.")) 
  

#Combined_for_multivariate[,c(7:10,7480:7482)] %>% pairs()
Combined_for_multivariate[,c(11:23,ncol(PTR_Kuster):ncol(PTR_Kuster)+5)] %>% pairs()
# ami_data <- read.table("http://static.lib.virginia.edu/statlab/materials/data/ami_data.DAT")
# names(ami_data) <- c("TOT","AMI","GEN","AMT","PR","DIAP","QRS")
# summary(ami_data)
# pairs(ami_data)

#mlm1 <- lm(cbind(TOT, AMI) ~ GEN + AMT + PR + DIAP + QRS, data = ami_data)
#summary(mlm1)
Ind <-  colnames(Combined_for_multivariate) %>% str_subset("^PC(1|2|3|4)")
Dep <-  colnames(Combined_for_multivariate) %>% str_subset("^PC[0-9]",negate= T) %>% unique()

fo <- sprintf("cbind(%s) ~ %s",toString(Dep),paste(Ind,collapse = "+"))
mlm2 <- do.call("lm",list(fo,quote(Combined_for_multivariate)))
testing<-list(PCs =  broom::tidy(mlm2), #%>% mutate(Pred_class = x),
     R2 = map_dbl(summary(mlm2),"adj.r.squared") %>% set_names(Dep) %>% stack() %>% set_names(c("R2","GeneName")))# %>% mutate(Pred_class = x))
testing$R2 %>% subset(R2>0) %>%  left_join(KEGG_genes %>% dplyr::select(-Uniprot) %>% 
            inner_join(KEGG_pathways %>%
                         subset(Path_type =="metabolic") %>% 
                         dplyr::select(Path_id,Path_description), by = c("pathway" = "Path_id")
            ) %>% distinct(), by = c("GeneName" = "ID")) %>%
  group_by(pathway) %>% 
  add_count(name="Genes_per_pathway") %>% 
  subset(#pathway %in% unique(.$pathway)[1:50] &
    Genes_per_pathway>15) %>% ungroup %>%
  ggplot(aes(x = R2, y = Path_description,fill = Path_description))+
  ggridges::geom_density_ridges(rel_min_height = 0.01,alpha = 0.5)+
  geom_vline( xintercept = 0) +
  ggtitle("Metabolism PTR Across Cancer tissues", 
          subtitle = "subsetted for more that 15 genes per pathway")
highly_predicted <- testing$R2 %>% subset(R2>0.9999) %>% pull(GeneName) %>% as.character()
testing$PCs %>% subset(response %in% highly_predicted) %>%
  pivot_wider(id_cols = "response", names_from = "term", values_from = "estimate") %>% 
  column_to_rownames("response") %>%
ComplexHeatmap::Heatmap()
#testing$R2$R2[testing$R2$R2<0]<- 0  
#mlm2 <- lm.fit(cbind(Combined_for_multivariate[,1:3]),Combined_for_multivariate[,Dep)
#R2 <- map_dbl(summary(mlm2),"adj.r.squared") %>% set_names(Dep)

#R2[R2>0.75]

hist_of_R2_safe <- safely(hist_of_R2)
library(future)
plan(multisession, workers = 4)
Reactome_pathways_of_interest <- humanreactome_df %>% pull("Reactome_Pathway") %>% unique() %>% 
    str_subset(pattern = "tRNA|Translation |roteasome|Transcription|plicing")

testing3_f <- future::future({map(.x  = Reactome_pathways_of_interest,
        ~hist_of_R2_safe(.x))},
        seed = T)
testing3 <- value(testing3_f) %>% set_names(Reactome_pathways_of_interest) #%>% 
testing34  <- testing3 %>%   map("result") %>%  compact()
heatmap_df <- testing34 %>%  map_df("R2") %>% pivot_wider(names_from =Pred_class, values_from = R2)
plan(sequential)

testing4 <- testing3 %>% reduce(cbind)  
heatmap_df %>% column_to_rownames("Predicted") %>% dplyr::select((matches("tRNA| Eukaryotic Translation|Proteasome |Metabolism of RNA|Gene expression |mRNA Splicing|Citric Acid |Cell Cycle |Apoptosis", ignore.case = FALSE))) %>%  cor(use = "pairwise.complete.obs") %>%round(2) %>% pheatmap::pheatmap(main = "Corr matrix between samples", show_rownames = F)
heatmap_df %>% column_to_rownames("Predicted") %>% dplyr::select((matches("tRNA| Eukaryotic Translation|Proteasome |Metabolism of RNA|Gene expression |mRNA Splicing|Citric Acid |Cell Cycle |Apoptosis", ignore.case = FALSE))) %>%  pheatmap::pheatmap(main = "R2", show_rownames = F,  clustering_distance_rows = "euclidean",  clustering_distance_cols =  "correlation")

testing4 <- testing4 %>% as.data.frame()%>% rownames_to_column("ID") %>% left_join(KEGG_genes[,-2] )  %>% 
    mutate(pathways = case_when(pathway == "hsa00020" ~ "TCA Cycle",
                                pathway == "hsa00520" ~ "Sugar and Amino acid",
                                pathway == "hsa00240" ~ "Pyrimidine metabolism",
                                pathway == "hsa00230" ~ "Purine metabolism",
                                TRUE ~"ELSE")) %>% distinct()
testing4 %>%    ggplot(aes(x=`Mito Ribosome`, colour = pathway)) +
    geom_histogram(alpha=0.8, position="stack", aes(y=..density..)) + theme(legend.position = "none") +
    facet_wrap(~ pathway)
testing345 <- testing34 %>%  map("R2") %>% .[str_detect(names(.),"tRNA| Eukaryotic Translation|Proteasome |Metabolism of RNA|Gene expression |mRNA Splicing|Citric Acid |Cell Cycle |Apoptosis")]
testing345 <- testing345 %>% imap(.x = .,~.x %>% mutate(in_class = (as.character(Predicted) %in% (subset(Reactome_pathways,Reactome_Pathway == .y) %>% pull(Gene_name)))))%>% purrr::reduce(rbind)  
    
ggplot(testing345, aes(x=R2, fill = in_class)) +
    geom_histogram(alpha=0.6, position="dodge" , aes(y=..density..)) + #theme(legend.position = "none") +
    facet_wrap(~ Pred_class)


pca_fit_mito <- CCLE_proteins%>% 
    subset(.,rownames(.) %in% (humanreactome_df %>%  subset(str_detect(Reactome_Pathway, "itochondrial")) %>%                                    pull(Uniprot))) %>%  # (uniprot_factors %>% subset(str_detect(Type,"plice" )) %>% pull(Entry))) %>%
    subset(is.na(.) %>% rowSums<(ncol(.)/4))%>% as.matrix() %>%  
    impute::impute.knn() %>% .[["data"]] %>% t()%>% prcomp()
pca_fit_cell_cycle <- CCLE_proteins%>% 
    subset(.,rownames(.) %in% (humanreactome_df %>%  subset(Reactome_Pathway == "Cell Cycle Checkpoints") %>% 
                                   pull(Uniprot))) %>%  # (uniprot_factors %>% subset(str_detect(Type,"plice" )) %>% pull(Entry))) %>%
    subset(is.na(.) %>% rowSums<(ncol(.)/4))%>% as.matrix() %>%  
    impute::impute.knn() %>% .[["data"]] %>% t()%>% prcomp()

pca_fit_splicing <- CCLE_proteins%>% 
    subset(.,rownames(.) %in% (humanreactome_df %>%  subset(Reactome_Pathway == "mRNA Splicing") %>% 
                                   pull(Uniprot))) %>%  # (uniprot_factors %>% subset(str_detect(Type,"plice" )) %>% pull(Entry))) %>%
    subset(is.na(.) %>% rowSums<(ncol(.)/4))%>% as.matrix() %>%  
    impute::impute.knn() %>% .[["data"]] %>% t()%>% prcomp()


Combined_for_multivariate <- inner_join(PTR_CCLE_na_rm_2,
                                        pca_fit_mito$x[,1:3] %>% as.data.frame() %>%  rownames_to_column("Cell_line")) %>%
    #inner_join(pca_fit_cell_cycle$x[,1:3] %>% as.data.frame() %>%  rownames_to_column("Cell_line"), by = "Cell_line") %>%
    #inner_join(pca_fit_splicing$x[,1:3] %>% as.data.frame() %>%  rownames_to_column("Cell_line"), by = "Cell_line") %>% 
    column_to_rownames("Cell_line") %>% as.matrix() %>% impute::impute.knn() %>% .[["data"]] %>% as.data.frame() %>% 
    set_names(.,str_replace(colnames(.),"-", "_"))


#Combined_for_multivariate[,c(7:10,7480:7482)] %>% pairs()
#Combined_for_multivariate[,c(11:23,7480:7482)] %>% pairs()
# ami_data <- read.table("http://static.lib.virginia.edu/statlab/materials/data/ami_data.DAT")
# names(ami_data) <- c("TOT","AMI","GEN","AMT","PR","DIAP","QRS")
# summary(ami_data)
# pairs(ami_data)

#mlm1 <- lm(cbind(TOT, AMI) ~ GEN + AMT + PR + DIAP + QRS, data = ami_data)
#summary(mlm1)
Ind <-  colnames(Combined_for_multivariate) %>% .[(ncol(Combined_for_multivariate)-2):ncol(Combined_for_multivariate)]
Dep <-  colnames(Combined_for_multivariate) %>% .[1:(ncol(Combined_for_multivariate)-3)]
fo <- sprintf("cbind(%s) ~ %s",toString(Dep),paste(Ind,collapse = "+"))
mlm2 <- do.call("lm",list(fo,quote(Combined_for_multivariate)))
testing456 <- map_dbl(summary(mlm2),"adj.r.squared") %>% set_names(Dep) %>% stack() %>% set_names(c("R2","Predicted")) %>% mutate(in_class = Predicted %in% (humanreactome_df %>% 
                                                                                                                                                                         subset(str_detect(Reactome_Pathway,"itochondrial")) %>% 
                                                                                                                                                                         pull(Gene_name)),
                                                                                                                                  in_euk_tra = Predicted %in% (humanreactome_df %>% 
                                                                                                                                                                 subset(str_detect(Reactome_Pathway,"tic Translation ")) %>% 
                                                                                                                                                                 pull(Gene_name)),
                                                                                                                                  in_TCA = Predicted %in% (humanreactome_df %>% 
                                                                                                                                                               subset(str_detect(Reactome_Pathway,"Citric acid cycle ")) %>% 
                                                                                                                                                               pull(Gene_name)),
                                                                                                                                  in_amino_acid = Predicted %in% (humanreactome_df %>% 
                                                                                                                                                               subset(str_detect(Reactome_Pathway,"Nucleotide salvage")) %>% 
                                                                                                                                                               pull(Gene_name)),
                                                                                                                                  MT_= str_detect(Predicted,"MT_"),
                                                                                                                                  in_MT_tr = Predicted %in% (humanreactome_df %>% 
                                                                                                                                                                  subset(str_detect("Reactome_Pathway,Mitochondrial tRNA aminoacylation")) %>% 
                                                                                                                                                               pull(Gene_name)))
testing456 %>% ggplot(aes(x = R2, fill= MT_)) + geom_histogram( position = 'dodge', aes(y = ..density..))
list(PCs =  broom::tidy(mlm2) %>% mutate(Pred_class = x),
     R2 = map_dbl(summary(mlm2),"adj.r.squared") %>% set_names(Dep) %>% stack() %>% set_names(c("R2","Predicted")) %>% mutate(Pred_class = x))
################################
#RNA-seq vs Protein vs PTR
Combined_for_multivariate_RNA_seq <- inner_join(CCLE_RNA_seq %>% t() %>% as.data.frame() %>%  rownames_to_column("Cell_line"),
                                        pca_fit_eEF$x[,1:3] %>% as.data.frame() %>%  rownames_to_column("Cell_line")) %>%
    column_to_rownames("Cell_line") %>% as.matrix() %>% impute::impute.knn() %>% .[["data"]] %>% as.data.frame() %>% 
    set_names(.,str_remove_all(colnames(.),"-"))

Combined_for_multivariate_Prot <- inner_join(CCLE_proteins %>% t() %>% as.data.frame() %>%  rownames_to_column("Cell_line"),
                                                pca_fit_eEF$x[,1:3] %>% as.data.frame() %>%  rownames_to_column("Cell_line")) %>%
    column_to_rownames("Cell_line") %>% as.matrix() %>% impute::impute.knn() %>% .[["data"]] %>% as.data.frame() %>% 
    set_names(.,str_remove_all(colnames(.),"-"))    

Combined_for_multivariate_PTR <- inner_join(PTR_CCLE_na_rm_2,
                                        pca_fit_eEF$x[,1:3] %>% as.data.frame() %>%  rownames_to_column("Cell_line")) %>%
      column_to_rownames("Cell_line") %>% as.matrix() %>% impute::impute.knn() %>% .[["data"]] %>% as.data.frame() %>% 
    set_names(.,str_remove_all(colnames(.),"-"))

multivariate_datasets <- function(df,y){
    Ind <-  colnames(df) %>% .[(ncol(df)-8):ncol(df)]
    Dep <-  colnames(df) %>% .[1:(ncol(df)-9)]
    fo <- sprintf("cbind(%s) ~ %s",toString(Dep),paste(Ind,collapse = "+"))
    mlm2 <- do.call("lm",list(fo,quote(df)))
    map_dbl(summary(mlm2),"adj.r.squared") %>% set_names(Dep) %>% stack() %>% set_names(c("R2","Predicted")) %>% mutate(omic = y)
}
omics_eEFs <- imap(list(RNA_seq = Combined_for_multivariate_RNA_seq,Prot = Combined_for_multivariate_Prot,PTR =  Combined_for_multivariate_PTR), multivariate_datasets)
omics_eEFs[["Prot"]] <- omics_eEFs[["Prot"]] %>% left_join(HUMAN_9606_idmapping %>% subset(Type == "Gene_Name"), by = c("Predicted" = "Uniprot")) %>% dplyr::select(R2,ID,omic) %>% set_names(c("R2","Predicted","omic"))
omics_eEFs_2 <- omics_eEFs %>% purrr::reduce(rbind) %>% distinct(Predicted,omic,.keep_all = T) %>% pivot_wider(names_from = omic,values_from = R2)
omics_eEFs_3 <- omics_eEFs_2 %>% subset((PTR>Prot) &(PTR>RNA_seq))  %>% mutate(Diff = PTR - Prot )


####
testing <- models %>% map("result") %>% compact() %>% imap_dfr(.x = ., ~.x[["final_fitted"]] %>% subset(.metric =="rsq") %>% mutate(Uniprot = .y))
testing$.estimate %>% hist()              
