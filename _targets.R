pacman::p_load(targets,openxlsx,tarchetypes,ComplexHeatmap, magrittr,future,readr)
source(here::here("Codes","functions_minimal.R"))
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("biglm", "tidyverse"))
#in the console run targets::tar_visnetwork()
#and targets::tar_make()
list(
    tar_target(
        sample_info_file, #Raw from https://depmap.org/portal/download/ downloaded 2/11/2020 
        here::here("Datasets","Raw","sample_info.csv"),
        format = "file"),
    tar_target(
        uniprot_factors, #Downloaded from uniprot by searching for specific
        purrr::map_df(.x = c("ribosome","proteasome","membrane","spliceosome","eukaryotic_translation"), 
                                 ~read_tsv(file = here::here("Datasets","Raw",paste0("uniprot_",.x,"_factor.tab"))) %>% 
                                     mutate(Type = .x))),
    tar_target(
        gene_clusters_eIFs_file,
        here::here("Datasets","Processed","eIFS_gene_clusters.csv"),
        format = "file"),
    tar_target(
        gene_clusters_eIFs,
        read_csv(gene_clusters_eIFs_file)),
    tar_target(
        sample_info,
        read.csv(sample_info_file, stringsAsFactors = FALSE)),
    tar_target(
        Metab_CCLE_file, ####Raw from https://portals.broadinstitute.org/ccle/data CCLE_metabolomics_20190502.csv 2/11/2020
        here::here("Datasets","Raw","CCLE_metabolites_landscape_of_cancer.xlsx"),
        format = "file"),
    tar_target(
        Metab_CCLE,
        openxlsx::read.xlsx(Metab_CCLE_file,
                            sheet = "1-clean data")[-760,] %>% 
            mutate(X1 = str_match(X1,"^([:graph:]*?)_")[,2]) %>%
            column_to_rownames("X1")%>% t()    ),
    tar_target(
        CCLE_prot_file, #Raw from https://gygi.med.harvard.edu/publications/ccle 2/11/2020
        here::here("Datasets","Raw","Table_S2_Protein_Quant_Normalized.xlsx"),
        format = "file"),
    tar_target(
        CCLE_proteins,
        loading_CCLE_prot(x= CCLE_prot_file)),
    tar_target(
        CCLE_RNA_seq_file_original,#Raw from https://depmap.org/portal/download/ 2/11/2020
        here::here("Datasets","Raw","CCLE_expression.csv"),
        format = "file"),
    tar_target(
        CCLE_RNA_seq_file,#Raw from https://portals.broadinstitute.org/ccle/data 3/18/2020
        here::here("Datasets","Raw","CCLE_RNAseq_rsem_genes_tpm_20180929.txt"),
        format = "file"),
    tar_target(
        CCLE_RNA_seq,#Raw from https://portals.broadinstitute.org/ccle/data 3/18/2020
        loading_CCLE_TPM(CCLE_RNA_seq_file)),
    tar_target(
        CCLE_RNA_seq_original,
        inner_join(sample_info[,1:2],read.csv(CCLE_RNA_seq_file_original), by = c("DepMap_ID" = "X"))[,-1] %>%
            column_to_rownames(var = "stripped_cell_line_name")  %>% t()  %>% 
            magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2]) %>% 
            .[!(rownames(.) %>% duplicated),] %>%
            as_tibble(.,rownames = NA) %>% rownames_to_column("Gene") %>% 
            #mutate(across(where(is.numeric), ~na_if(., 0)))%>% 
            column_to_rownames("Gene") %>% 
            subset(.,((. == 0) %>% rowSums())<(ncol(.)/100))%>%
            `+`(.,0.01)%>%
            `/`(.,matrixStats::rowMedians(as.matrix(.),na.rm = T)) %>% log10() %>% as.matrix) ,    
    tar_target(
        NCI_60_metabolites_file,#Raw from  https://www.nature.com/articles/s41467-019-09695-9#Sec28 2/11/2020
        here::here("Datasets","Raw","41467_2019_9695_MOESM2_ESM.xlsx"),
        format = "file"),
    tar_target(
        NCI_60_metabolites,
        openxlsx::read.xlsx(NCI_60_metabolites_file, sheet = 3, startRow = 4, na.strings = "NaN") %>%
            .[str_detect(.$`Annotation.ID`,"H|C"),-c(1:3,5)] %>% setNames(str_remove_all(colnames(.), "^[:graph:]*?_"))  %>%
            na.omit()%>% separate_rows(Annotation.ID, sep = " ; ") %>%
            subset(str_detect(Annotation.ID, "^H|C") & !duplicated(Annotation.ID)) %>%
            mutate(Annotation.ID =  str_replace_all(Annotation.ID,"HMDB", "HMDB00")) ),
    tar_target(
        NCI_60_proteins_file,#Raw from https://depmap.org/portal/download/ 2/11/2020
        here::here("Datasets","Raw","1-s2.0-S2589004219304407-mmc2.xlsx"),
        format = "file"),
    tar_target(
        NCI_60_proteins,
        openxlsx::read.xlsx(NCI_60_proteins_file, sheet =6) %>%
            .[,-c(2:8)] %>% remove_rownames() %>%
            column_to_rownames(var = "protein.accession.number") %>%
            setNames(str_remove(colnames(.), "^[:graph:]*?_")) %>% log() %>% Z_transf()),
    tar_target(
        NCI_60_RNA_file,#Raw from https://discover.nci.nih.gov/cellminer/loadDownload.do 3/11/2020 - RNAseq Composite expression
        here::here("Datasets","Raw","RNA__RNA_seq_composite_expression.xls"),
        format = "file"),
    tar_target(
        NCI_60_RNA,
        readxl::read_xls(NCI_60_RNA_file, skip = 10) %>%
            .[,-c(2:6)] %>% setNames(str_remove_all(colnames(.), "^[:graph:]*:|-| ")) %>%
            column_to_rownames("Genenamed") %>%
            .[!(matrixStats::rowSds(as.matrix(.)) == 0),]),
    tar_target(
        Achilles_file,#Raw from https://depmap.org/portal/download/ 2/11/2020
        here::here("Datasets","Raw","Achilles_gene_effect.csv"),
        format = "file"),
    tar_target(
        Achilles,
        inner_join(sample_info[,1:2],read.csv(Achilles_file), by = "DepMap_ID")[,-1] %>%
        column_to_rownames(var = "stripped_cell_line_name") %>% t() %>%
        magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2])),
    tar_target(
        RNAi_file,#Raw from https://depmap.org/portal/download/ 2/11/2020
        here::here("Datasets","Raw","D2_combined_gene_dep_scores.csv"),
        format = "file"),
    tar_target(
        RNAi,
        read.csv(RNAi_file) %>% mutate(Gene_name = str_match(Gene_name,"^([:graph:]*?) ")[,2]) %>%
            column_to_rownames(var = "Gene_name") %>% 
            setNames(str_match(colnames(.), "^([:graph:]*?)_")[,2])    ),
    tar_target(
        HUMAN_9606_idmapping_file,#Raw from Uniprot
        here::here("Datasets","Raw","HUMAN_9606_idmapping.dat"),
        format = "file"),
    tar_target(
        CCLE_RNA_seq_mean_df_plot,
        Mean_sd_matrix(CCLE_RNA_seq) %>%
            ggplot(aes(x = Mean, y = SD)) +
            geom_point(size = 0.1,alpha = 0.2)+ 
            xlab("Mean across Cell lines") +
            ylab("SD across Cell lines")+
            ggtitle("CCLE_RNA_seq_mean_df_plot")),
    tar_target(
        CCLE_proteins_mean_df_plot,
        Mean_sd_matrix(CCLE_proteins) %>% 
            ggplot(aes(x = Mean, y = SD)) +
            geom_point(size = 0.1,alpha = 0.2)+ 
            xlab("Mean across Cell lines") +
            ylab("SD across Cell lines")+
            ggtitle("CCLE_proteins_mean_df_plot")),
    tar_target(
        KEGG_path_file,#Raw from Uniprot
        here::here("Datasets","Raw","hsa_pathways.txt"),
        format = "file"),
    tar_target(
        KEGG_pathways,
        read_tsv(KEGG_path_file)),
    tar_target(
        KEGG_genes,
        Retrieve_all_kegg_genes(KEGG_pathways %>% pull(Path_id))),
    tar_target(
        Master_pathways_KEGG,
        KEGG_genes %>% pivot_longer(-pathway,names_to = "Type",values_to = "ID") %>% 
            left_join(dplyr::select(KEGG_pathways, -Path_type), by = c("pathway" = "Path_id")) %>% 
            set_names(c("Path_ID","Type","ID","SubPathway","Pathway"))),
    tar_target(
        HUMAN_9606_idmapping,
        read_tsv(HUMAN_9606_idmapping_file,
                 col_names = FALSE) %>%
            setNames(c("Uniprot", "Type", "ID"))),
    tar_target(
        Metabolite_mapping_file,#Raw from Uniprot
        here::here("Datasets","Raw","CCLE_metabolites_mapped.csv"),
        format = "file"),
    tar_target(
        Metabolite_mapping,
        read_tsv(Metabolite_mapping_file) ),
    tar_target(
        humankegg_df,
        read_tsv(Metabolite_mapping_file) ),
    tar_target(
        CCLE_prot_Intra_omics_corr_matrix,
                   Intra_omic_Feature_correlation(CCLE_proteins %>% 
                                                      rownames_to_column("Uniprot")%>%
                                                      left_join(HUMAN_9606_idmapping%>% 
                                                                    subset(Type == "Gene_Name") %>% 
                                                                    distinct(Uniprot, .keep_all = T) %>% 
                                                                    dplyr::select(-Type)) %>%
                                                      dplyr::select(-Uniprot) %>% 
                                                      distinct(ID, .keep_all = T) %>% 
                                                      subset(!is.na(ID)) %>% 
                                                      remove_rownames()%>%
                                                      column_to_rownames("ID"))),
    tar_target(
        CCLE_RNA_Intra_omics_corr_matrix,
        Intra_omic_Feature_correlation(CCLE_RNA_seq %>%  #removes genes below 25% SD
                                           subset(.,(CCLE_RNA_seq %>% matrixStats::rowSds()) > 
                                                      (CCLE_RNA_seq %>% matrixStats::rowSds() %>% quantile(., c(.25)))))),
    tar_target(
        NCI_Metabo_Intra_omics_corr_matrix,
        Intra_omic_Feature_correlation(NCI_60_metabolites %>%
                                           column_to_rownames('Annotation.ID'))),
    tar_target(
        CCLE_Metabo_Intra_omics_corr_matrix,
        Intra_omic_Feature_correlation(Metab_CCLE)),
    tar_target(
        NCI_RNA_Intra_omics_corr_matrix,
        Intra_omic_Feature_correlation(NCI_60_RNA)),
    tar_target(
        NCI_prot_Intra_omics_corr_matrix,
        Intra_omic_Feature_correlation(NCI_60_proteins)),
    tar_target(
        CCLE_prot_EIFs_genes_heatmap,
        Convert_df_to_heatmap(Data_df = CCLE_prot_Intra_omics_corr_matrix, Heatmap_name = "eTFs Proteome against all proteins",
                              Subset_col = gene_clusters_eIFs$Genes, #vector of features/ samples to keep
                              Subset_row = uniprot_factors %>% subset(Type == "eukaryotic_translation") %>%
                                  pull(`Gene names`) %>% str_split(" ")  %>% unlist()%>% str_subset("^(EEF|EIF)"),
                              show_names = "rows",
                              Main_top_annotation = gene_clusters_eIFs %>% dplyr::select(Uniprot, Pathways) %>%
                                  subset(!duplicated(Uniprot)) %>%
                                  rename(Genes = Uniprot))),
    tar_target(
        CCLE_RNA_EIFs_genes_heatmap,
        Convert_df_to_heatmap(Data_df = CCLE_RNA_Intra_omics_corr_matrix, Heatmap_name = "eTFs Genes against all genes",
                              Subset_col = gene_clusters_eIFs$Genes, #vector of features/ samples to keep
                              Subset_row = uniprot_factors %>% subset(Type == "eukaryotic_translation") %>%
                                  pull(`Gene names`) %>% str_split(" ")  %>% unlist()%>% str_subset("^(EEF|EIF)"),
                              show_names = "rows",
                              Main_top_annotation = gene_clusters_eIFs %>% 
                                  dplyr::select(Genes, Pathways) %>% 
                                  subset(!duplicated(Genes)))),
    tar_target(
        CCLE_RNA_prot_aligned,
        CCLE_proteins %>% 
            rownames_to_column("Uniprot")%>%
            left_join(HUMAN_9606_idmapping%>% 
                          subset(Type == "Gene_Name") %>% 
                          distinct(Uniprot, .keep_all = T) %>% 
                          dplyr::select(-Type)) %>%
            dplyr::select(-Uniprot) %>% 
            distinct(ID, .keep_all = T) %>% 
            subset(!is.na(ID)) %>% 
            remove_rownames()%>%
            column_to_rownames("ID") %>% Aligning_two_matrices(CCLE_RNA_seq)),
    tar_target(
        PTR_CCLE,
        subtract(CCLE_RNA_prot_aligned[[1]],CCLE_RNA_prot_aligned[[2]])),
    tar_target(
        CCLE_RNA_prot_intracorr_aligned,
        Aligning_two_matrices(CCLE_prot_Intra_omics_corr_matrix,CCLE_RNA_Intra_omics_corr_matrix)),
    tar_target(
        CCLE_prot_RNA_intracorr_half_matrix,
        add(get_lower_tri(CCLE_RNA_prot_intracorr_aligned[[1]]) %>% replace_na(0),
            get_upper_tri(CCLE_RNA_prot_intracorr_aligned[[2]]) %>% replace_na(0))),
    tar_target(
        CCLE_prot_RNA_half_heatmap,
        Convert_df_to_heatmap(Data_df = CCLE_prot_RNA_intracorr_half_matrix, Heatmap_name = "Half_matrix",
                              Subset_col = Master_pathways_KEGG$ID, #KEGG metabolism genes
                              Subset_row = Master_pathways_KEGG$ID,
                              cluster = "none",
                              show_names = "none", 
                              row_order = subset(Master_pathways_KEGG, ID %in% rownames(CCLE_prot_RNA_intracorr_half_matrix)) %>% pull(ID) %>% unique(),
                              column_order = subset(Master_pathways_KEGG, ID %in% colnames(CCLE_prot_RNA_intracorr_half_matrix)) %>% pull(ID) %>% unique())),
    
    
    tar_target(
        PTR_mean_df_plot,
        Mean_sd_matrix(PTR_CCLE) %>%
            ggplot(aes(x = Mean, y = SD)) +
            geom_point(size = 0.1,alpha = 0.2)+ 
            xlab("Mean across Cell lines") +
            ylab("SD across Cell lines")+
            ggtitle("PTR_df_plot")),
    tar_target(
        eTFs_PTR_corr,
        Factors_to_PTR_corr(uniprot_factors %>% subset(str_detect(Type,"eukaryotic_translation")) %>%
                                pull(Entry),PTR_CCLE, CCLE_proteins)),
    tar_render(report, here::here("Output","report.Rmd"))
)
