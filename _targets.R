pacman::p_load(targets,openxlsx,tarchetypes,ComplexHeatmap, magrittr,future,readr,furrr,progressr,tidymodels,xaringan, vip,OmnipathR)
source(fs::path_rel(here::here("Codes","functions_minimal.R")))
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("biglm", "tidyverse"))
#in the console run targets::tar_visnetwork()
#and targets::tar_make()
list(
    tar_target(
        sample_info_file, #Raw from https://depmap.org/portal/download/ downloaded 2/11/2020 
        fs::path_rel(here::here("Datasets","Raw","sample_info.csv")),
        format = "file"),
    tar_target(
        uniprot_factors, #Downloaded from uniprot by searching for specific
        purrr::map_df(.x = c("ribosome","proteasome","membrane","spliceosome","eukaryotic_translation"), 
                      ~read_tsv(file = fs::path_rel(here::here("Datasets","Raw",paste0("uniprot_",.x,"_factor.tab")))) %>% 
                          mutate(Type = .x))),
    tar_target(
        humanreactome_df_file,
        fs::path_rel(here::here("Datasets","Raw","Reactome_pathways.csv")),
        format = "file"),
    tar_target(
        humanreactome_df,
        read_csv2(humanreactome_df_file)[,-1]),
    tar_target(
        gene_clusters_eIFs_file,
        fs::path_rel(here::here("Datasets","Processed","eIFS_gene_clusters.csv")),
        format = "file"),
    tar_target(
        gene_clusters_eIFs,
        read_csv(gene_clusters_eIFs_file)),
    tar_target(
        sample_info,
        read.csv(sample_info_file, stringsAsFactors = FALSE)),
    tar_target(
        Metab_CCLE_file, ####Raw from https://portals.broadinstitute.org/ccle/data CCLE_metabolomics_20190502.csv 2/11/2020
        fs::path_rel(here::here("Datasets","Raw","CCLE_metabolites_landscape_of_cancer.xlsx")),
        format = "file"),
    tar_target(
        Metab_CCLE,
        openxlsx::read.xlsx(Metab_CCLE_file,
                            sheet = "1-clean data")[-760,] %>% 
            mutate(X1 = str_match(X1,"^([:graph:]*?)_")[,2]) %>%
            column_to_rownames("X1")%>% t()    ),
    tar_target(
        mRNA_PTR_Protein_Kuster_file, #Raw from https://www.embopress.org/doi/full/10.15252/msb.20188513 3/5/2021
        fs::path_rel(here::here("Datasets","Raw","PTR_Healthy_tissue_Kuster_2019.txt")),
        format = "file"),
    tar_target(
        mRNA_PTR_Protein_Kuster,
        loading_Kuster(mRNA_PTR_Protein_Kuster_file)),
    tar_target(
        mRNA_Kuster,
        mRNA_PTR_Protein_Kuster[[1]]),
    tar_target(
        Protein_Kuster,
        mRNA_PTR_Protein_Kuster[[2]]),
    tar_target(
        PTR_Kuster,
        #mRNA_PTR_Protein_Kuster[[3]] %>% `/`(.,matrixStats::rowMedians(as.matrix(.),na.rm = T)) %>% log2() %>% as.matrix 
        Protein_Kuster - mRNA_Kuster),
    tar_target(
        CCLE_prot_file, #Raw from https://gygi.med.harvard.edu/publications/ccle 2/11/2020
        fs::path_rel(here::here("Datasets","Raw","Table_S2_Protein_Quant_Normalized.xlsx")),
        format = "file"),
    tar_target(
        CCLE_proteins,
        loading_CCLE_prot(x= CCLE_prot_file)),
    tar_target(
        CCLE_RNA_seq_file_original,#Raw from https://depmap.org/portal/download/ 2/11/2020
        fs::path_rel(here::here("Datasets","Raw","CCLE_expression.csv")),
        format = "file"),
    tar_target(
        CCLE_RNA_seq_file,#Raw from https://portals.broadinstitute.org/ccle/data 3/18/2020
        fs::path_rel(here::here("Datasets","Raw","CCLE_RNAseq_rsem_genes_tpm_20180929.txt")),
        format = "file"),
    tar_target(
        CCLE_RNA_seq,#Raw from https://portals.broadinstitute.org/ccle/data 3/18/2020
        loading_CCLE_TPM(CCLE_RNA_seq_file,HUMAN_9606_idmapping)),
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
        fs::path_rel(here::here("Datasets","Raw","41467_2019_9695_MOESM2_ESM.xlsx")),
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
        fs::path_rel(here::here("Datasets","Raw","1-s2.0-S2589004219304407-mmc2.xlsx")),
        format = "file"),
    tar_target(
        NCI_60_proteins,
        openxlsx::read.xlsx(NCI_60_proteins_file, sheet =6) %>%
            .[,-c(2:8)] %>% remove_rownames() %>%
            column_to_rownames(var = "protein.accession.number") %>%
            setNames(str_remove(colnames(.), "^[:graph:]*?_")) %>% log() %>% Z_transf()),
    tar_target(
        NCI_60_RNA_file,#Raw from https://discover.nci.nih.gov/cellminer/loadDownload.do 3/11/2020 - RNAseq Composite expression
        fs::path_rel(here::here("Datasets","Raw","RNA__RNA_seq_composite_expression.xls")),
        format = "file"),
    tar_target(
        NCI_60_RNA,
        readxl::read_xls(NCI_60_RNA_file, skip = 10) %>%
            .[,-c(2:6)] %>% setNames(str_remove_all(colnames(.), "^[:graph:]*:|-| ")) %>%
            column_to_rownames("Genenamed") %>%
            .[!(matrixStats::rowSds(as.matrix(.)) == 0),]),
    tar_target(
        Achilles_file,#Raw from https://depmap.org/portal/download/ 2/11/2020
        fs::path_rel(here::here("Datasets","Raw","Achilles_gene_effect.csv")),
        format = "file"),
    tar_target(
        Achilles,
        inner_join(sample_info[,1:2],read.csv(Achilles_file), by = "DepMap_ID")[,-1] %>%
            column_to_rownames(var = "stripped_cell_line_name") %>% t() %>%
            magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2])),
    tar_target(
        RNAi_file,#Raw from https://depmap.org/portal/download/ 2/11/2020
        fs::path_rel(here::here("Datasets","Raw","D2_combined_gene_dep_scores.csv")),
        format = "file"),
    tar_target(
        RNAi,
        utils::read.csv(RNAi_file) %>% mutate(Gene_name = str_match(Gene_name,"^([:graph:]*?) ")[,2]) %>%
            column_to_rownames(var = "Gene_name") %>% 
            setNames(str_match(colnames(.), "^([:graph:]*?)_")[,2])    ),
    tar_target(
        HUMAN_9606_idmapping_file,#Raw from Uniprot
        fs::path_rel(here::here("Datasets","Raw","HUMAN_9606_idmapping.dat")),
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
        Kuster_PTR_mean_df_plot,
        Mean_sd_matrix(PTR_Kuster) %>%
            ggplot(aes(x = Mean, y = SD)) +
            geom_point(size = 0.1,alpha = 0.2)+ 
            xlab("Mean across Cell lines") +
            ylab("SD across Cell lines")+
            ggtitle("Kuster_PTR_mean_df_plot")),
    tar_target(
        Kuster_RNA_seq_mean_df_plot,
        Mean_sd_matrix(mRNA_Kuster) %>%
            ggplot(aes(x = Mean, y = SD)) +
            geom_point(size = 0.1,alpha = 0.2)+ 
            xlab("Mean across Cell lines") +
            ylab("SD across Cell lines")+
            ggtitle("Kuster_RNA_seq_mean_df_plot")),
    tar_target(
        Kuster_prot_mean_df_plot,
        Mean_sd_matrix(Protein_Kuster) %>%
            ggplot(aes(x = Mean, y = SD)) +
            geom_point(size = 0.1,alpha = 0.2)+ 
            xlab("Mean across Cell lines") +
            ylab("SD across Cell lines")+
            ggtitle("Kuster_prot_mean_df_plot")),
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
        fs::path_rel(here::here("Datasets","Raw","hsa_pathways.txt")),
        format = "file"),
    tar_target(
        KEGG_pathways,
        readr::read_tsv(KEGG_path_file) %>%
            mutate(Path_description = str_remove(Path_description, " - Homo sapiens \\(human\\)"))),
    tar_target(
        KEGG_genes_file,
        fs::path_rel(here::here("Datasets","Raw","KEGG_genes.csv")),
        format = "file"),
    tar_target(
        KEGG_genes,
        read_csv(KEGG_genes_file)),
    tar_target(
        Master_pathways_KEGG,
        KEGG_genes %>% pivot_longer(-pathway,names_to = "Type",values_to = "ID") %>% 
            left_join(dplyr::select(KEGG_pathways, -Path_type), by = c("pathway" = "Path_id")) %>% 
            set_names(c("Path_ID","Type","ID","SubPathway","Pathway"))),
    tar_target(
        HUMAN_9606_idmapping,
        readr::read_tsv(HUMAN_9606_idmapping_file,
                        col_names = FALSE) %>%
            setNames(c("Uniprot", "Type", "ID"))),
    tar_target(
        Metabolite_mapping_file,#Raw from Uniprot
        fs::path_rel(here::here("Datasets","Raw","CCLE_metabolites_mapped.csv")),
        format = "file"),
    tar_target(
        Metabolite_mapping,
        read.csv(Metabolite_mapping_file) ),
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
        mRNA_Kuster_Intra_omics_corr_matrix,
        Intra_omic_Feature_correlation(mRNA_Kuster %>%  #removes genes below 25% SD
                                           subset(.,(mRNA_Kuster %>% matrixStats::rowSds(na.rm = T)) > 
                                                      (mRNA_Kuster %>% matrixStats::rowSds(na.rm = T) %>% quantile(., c(.25),na.rm = T))))),
    tar_target(
        Protein_Kuster_Intra_omics_corr_matrix ,
        Intra_omic_Feature_correlation(Protein_Kuster %>%  #removes genes below 25% SD
                                           subset(.,(Protein_Kuster %>% matrixStats::rowSds(na.rm = T)) > 
                                                      (Protein_Kuster %>% matrixStats::rowSds(na.rm = T) %>% quantile(., c(.25),na.rm = T))))),
    tar_target(
        PTR_Kuster_Intra_omics_corr_matrix ,
        Intra_omic_Feature_correlation(PTR_Kuster %>%  #removes genes below 25% SD
                                           subset(.,(PTR_Kuster %>% matrixStats::rowSds(na.rm = T)) > 
                                                      (PTR_Kuster %>% matrixStats::rowSds(na.rm = T) %>% quantile(., c(.25),na.rm = T))))),
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
        CCLE_RNA_intracor_heatmap,
        Convert_df_to_heatmap(Data_df = CCLE_RNA_Intra_omics_corr_matrix, Heatmap_name = "RNA intracor_matrix",
                              Subset_col = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino")) %>% pull(ID), #KEGG metabolism genes
                              Subset_row = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")& str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino")) %>% pull(ID),
                              cluster = "none",
                              show_names = "none",
                              Main_top_annotation = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")) %>% 
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>% 
                                  distinct(Genes,.keep_all = T),
                              Main_side_annotation = Master_pathways_KEGG %>%
                                  subset(str_detect(Pathway,"Metabolism")) %>%
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>%
                                  distinct(Genes,.keep_all = T),
                              row_order =  Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(CCLE_RNA_Intra_omics_corr_matrix))) %>% 
                                  pull(ID) %>% unique(),
                              column_order = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(CCLE_RNA_Intra_omics_corr_matrix))) %>% 
                                  pull(ID) %>% unique())),
    tar_target(
        CCLE_prot_intracor_heatmap,
        Convert_df_to_heatmap(Data_df = CCLE_prot_Intra_omics_corr_matrix, Heatmap_name = "Proteomic intracor_matrix",
                              Subset_col = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino")) %>% pull(ID), #KEGG metabolism genes
                              Subset_row = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")& str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino")) %>% pull(ID),
                              cluster = "none",
                              show_names = "none",
                              Main_top_annotation = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")) %>% 
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>% 
                                  distinct(Genes,.keep_all = T),
                              Main_side_annotation = Master_pathways_KEGG %>%
                                  subset(str_detect(Pathway,"Metabolism")) %>%
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>%
                                  distinct(Genes,.keep_all = T),
                              row_order =  Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(CCLE_prot_Intra_omics_corr_matrix))) %>% 
                                  pull(ID) %>% unique(),
                              column_order = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(CCLE_prot_Intra_omics_corr_matrix))) %>% 
                                  pull(ID) %>% unique())),
    
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
        PTR_TCA_CCLE,
        PTR_CCLE %>% as.data.frame %>%  rownames_to_column("GeneName") %>% 
            pivot_longer(cols = -GeneName,names_to = "Cell_line",values_to= "PTR") %>% 
            left_join(sample_info[,c("stripped_cell_line_name","lineage")], 
                      by = c("Cell_line" = "stripped_cell_line_name")) %>% 
            group_by(GeneName,lineage) %>% 
            summarise(PTR = median(PTR,na.rm = T)) %>% na.omit() %>% 
            mutate(pathway = case_when(GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00020") %>% 
                                                          pull(ID))~ "TCA",
                                       GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00190") %>% 
                                                          pull(ID))~ "OXPHOS",
                                       TRUE~  "OTHER")) %>% 
            ggplot(aes(x = PTR, fill = pathway))+ geom_density(alpha = 0.5)+
            ggtitle("TCA PTR Across Cancer tissues")+facet_wrap("lineage")),
    tar_target(
        PTR_TCA_Kuster,
        PTR_Kuster %>% as.data.frame %>%  rownames_to_column("GeneName") %>% 
            pivot_longer(cols = -GeneName,names_to = "Tissue",values_to= "PTR") %>% 
            mutate(pathway = case_when(GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00020") %>% 
                                                          pull(ID))~ "TCA",
                                       GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00190") %>% 
                                                          pull(ID))~ "OXPHOS",
                                       TRUE~  "OTHER")) %>% 
            ggplot(aes(x = PTR, fill = pathway))+ geom_density(alpha = 0.5)+
            ggtitle("TCA PTR Across Healthy tissues")+facet_wrap("Tissue")),
    
    tar_target(
        CCLE_RNA_prot_intracorr_aligned,
        Aligning_two_matrices(CCLE_prot_Intra_omics_corr_matrix,CCLE_RNA_Intra_omics_corr_matrix)),
    tar_target(
        CCLE_prot_RNA_intracorr_half_matrix_ordered,
        CCLE_RNA_prot_intracorr_aligned %>% purrr::map(.x = .,~.x[Master_pathways_KEGG %>% 
                                                                      subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Purine|Oxidative|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(.x))) %>% 
                                                                      pull(ID) %>% unique(),Master_pathways_KEGG %>% 
                                                                      subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Purine|Oxidative|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(.x))) %>% 
                                                                      pull(ID) %>% unique()]) %$% 
            add(get_lower_tri(.[[1]]) %>% replace_na(0),
                get_upper_tri(.[[2]]) %>% replace_na(0))),
    tar_target(
        CCLE_prot_RNA_half_heatmap,
        Convert_df_to_heatmap(Data_df = CCLE_prot_RNA_intracorr_half_matrix_ordered, Heatmap_name = "CCLE_Half_matrix",
                              Subset_col = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino")) %>% pull(ID), #KEGG metabolism genes
                              Subset_row = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")& str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino")) %>% pull(ID),
                              cluster = "none",
                              show_names = "none",
                              Main_top_annotation = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")) %>% 
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>% 
                                  distinct(Genes,.keep_all = T),
                              Main_side_annotation = Master_pathways_KEGG %>%
                                  subset(str_detect(Pathway,"Metabolism")) %>%
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>%
                                  distinct(Genes,.keep_all = T),
                              row_order = NULL,
                              column_order = NULL)),
    tar_target(
        Kuster_RNA_prot_intracorr_aligned,
        Aligning_two_matrices(Protein_Kuster_Intra_omics_corr_matrix,mRNA_Kuster_Intra_omics_corr_matrix)),
    
    tar_target(
        Kuster_prot_RNA_intracorr_half_matrix_ordered,
        Kuster_RNA_prot_intracorr_aligned %>% purrr::map(.x = .,~.x[Master_pathways_KEGG %>% 
                                                                        subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(.x))) %>% 
                                                                        pull(ID) %>% unique(),Master_pathways_KEGG %>% 
                                                                        subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Oxidative|Purine|Pyrimidine|Pentose|Glycol|Amino") & (ID %in% rownames(.x))) %>% 
                                                                        pull(ID) %>% unique()]) %$% 
            add(get_lower_tri(.[[1]]) %>% replace_na(0),
                get_upper_tri(.[[2]]) %>% replace_na(0))),
    tar_target(
        Kuster_prot_RNA_half_heatmap,
        Convert_df_to_heatmap(Data_df = Kuster_prot_RNA_intracorr_half_matrix_ordered, Heatmap_name = "Kuster_Half_matrix",
                              Subset_col = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism") & str_detect(SubPathway,"TCA|Purine|Pyrimidine|Pentose|Glycol|Oxidative|Amino")) %>% pull(ID), #KEGG metabolism genes
                              Subset_row = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")& str_detect(SubPathway,"TCA|Purine|Pyrimidine|Pentose|Glycol|Oxidative|Amino")) %>% pull(ID),
                              cluster = "none",
                              show_names = "none",
                              Main_top_annotation = Master_pathways_KEGG %>% 
                                  subset(str_detect(Pathway,"Metabolism")) %>% 
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>% 
                                  distinct(Genes,.keep_all = T),
                              Main_side_annotation = Master_pathways_KEGG %>%
                                  subset(str_detect(Pathway,"Metabolism")) %>%
                                  dplyr::select(ID,SubPathway) %>% set_names(c("Genes","Pathways")) %>%
                                  distinct(Genes,.keep_all = T),
                              row_order = NULL,
                              column_order = NULL)),
    tar_target(
        mRNA_Prot_Cor_Kuster,
        mRNA_Kuster %>%
            as.data.frame() %>%  
            rownames_to_column("GeneName") %>% 
            pivot_longer(-GeneName,names_to = "Tissue",values_to = "mRNA") %>% 
            inner_join(Protein_Kuster %>%
                           as.data.frame() %>%  
                           rownames_to_column("GeneName") %>% 
                           pivot_longer(-GeneName,names_to = "Tissue",
                                        values_to = "Protein"), 
                       by = c("GeneName","Tissue")) %>%
            mutate(GeneName = factor(GeneName,unique(GeneName))) %>% 
            na.omit() %>% 
            group_by(GeneName) %>% 
            summarise(Gene_mRNA_prot = cor(mRNA,Protein)) %>% 
            na.omit() %>% 
            mutate(pathway = case_when(GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00020") %>% 
                                                          pull(ID))~ "TCA",
                                       GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00190") %>% 
                                                          pull(ID))~ "OXPHOS",
                                       TRUE~  "OTHER"))),
    tar_target(
        mRNA_Prot_Cor_Kuster_plot,
        ggplot(mRNA_Prot_Cor_Kuster,aes(x = Gene_mRNA_prot, fill = pathway))+
            geom_density(alpha = 0.5)+
            ggtitle("TCA mRNA protein Correlations Across Healthy tissues")),
    tar_target(
        mRNA_Prot_Cor_CCLE_within,
        CCLE_RNA_prot_aligned[[1]] %>%
            as.data.frame() %>%  
            rownames_to_column("GeneName") %>% 
            pivot_longer(-GeneName,names_to = "Cell_line",values_to = "mRNA") %>% 
            inner_join(CCLE_RNA_prot_aligned[[2]] %>%
                           as.data.frame() %>%  
                           rownames_to_column("GeneName") %>% 
                           pivot_longer(-GeneName,names_to = "Cell_line",
                                        values_to = "Protein"), 
                       by = c("GeneName","Cell_line")) %>%
            mutate(GeneName = factor(GeneName,unique(GeneName))) %>% 
            na.omit() %>% 
            left_join(sample_info[,c("stripped_cell_line_name","lineage")],by = c("Cell_line" = "stripped_cell_line_name")) %>%
            add_count(GeneName,lineage)%>%
            subset(n>5)%>%
            group_by(GeneName,lineage) %>% 
            summarise(Gene_mRNA_prot = cor(mRNA,Protein)) %>% 
            na.omit() %>% 
            mutate(pathway = case_when(GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00020") %>% 
                                                          pull(ID))~ "TCA",
                                       GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00190") %>% 
                                                          pull(ID))~ "OXPHOS",
                                       TRUE~  "OTHER"))),
    tar_target(
        mRNA_Prot_Cor_CCLE_within_plot,
        ggplot(mRNA_Prot_Cor_CCLE_within,aes(x = Gene_mRNA_prot, fill = pathway))+
            geom_density(alpha = 0.5)+
            facet_wrap("lineage")+
            ggtitle("pathway mRNA protein Correlations Within Cancer tissues")),
    tar_target(
        mRNA_Prot_Cor_CCLE,
        CCLE_RNA_prot_aligned[[1]] %>%
            as.data.frame() %>%  
            rownames_to_column("GeneName") %>% 
            pivot_longer(-GeneName,names_to = "Cell_line",values_to = "mRNA") %>% 
            inner_join(CCLE_RNA_prot_aligned[[2]] %>%
                           as.data.frame() %>%  
                           rownames_to_column("GeneName") %>% 
                           pivot_longer(-GeneName,names_to = "Cell_line",
                                        values_to = "Protein"), 
                       by = c("GeneName","Cell_line")) %>%
            mutate(GeneName = factor(GeneName,unique(GeneName))) %>% 
            na.omit() %>% 
            left_join(sample_info[,c("stripped_cell_line_name","lineage")],by = c("Cell_line" = "stripped_cell_line_name")) %>%
            add_count(GeneName,lineage)%>%
            subset(n>5)%>%
            group_by(GeneName,lineage) %>% 
            summarise(Average_mRNA_per_tissue = median(mRNA),
                      Average_protein_per_tissue = median(Protein)) %>%
            ungroup() %>%
            group_by(GeneName) %>%
            summarise(Gene_mRNA_prot = cor(Average_mRNA_per_tissue,Average_protein_per_tissue)) %>% 
            na.omit() %>% 
            mutate(pathway = case_when(GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00020") %>% 
                                                          pull(ID))~ "TCA",
                                       GeneName %in% (KEGG_genes %>%
                                                          subset(pathway == "hsa00190") %>% 
                                                          pull(ID))~ "OXPHOS",
                                       TRUE~  "OTHER"))),
    tar_target(
        mRNA_Prot_Cor_CCLE_plot,
        ggplot(mRNA_Prot_Cor_CCLE,aes(x = Gene_mRNA_prot, fill = pathway))+
            geom_density(alpha = 0.5)+
            #facet_wrap("lineage")+
            ggtitle("pathway mRNA protein Correlations between Cancer tissues")),
    tar_target(
        CCLE_Pathway_Cor_diff,
        Calcu_corr_diff_pathway(mRNA_Prot_Cor_CCLE,KEGG_genes)),
    tar_target(
        Kuster_Pathway_Cor_diff,
        Calcu_corr_diff_pathway(mRNA_Prot_Cor_Kuster,KEGG_genes)),
    tar_target(
        healthy_cancer_corr_pathways,
        purrr::map_dfr(.x = intersect(names(Kuster_Pathway_Cor_diff),
                                      names(CCLE_Pathway_Cor_diff)),
                       ~data.frame(Path_id = .x,
                                   p_value = ks.test(Kuster_Pathway_Cor_diff[[.x]],
                                                     CCLE_Pathway_Cor_diff[[.x]]) %>% 
                                       pluck("p.value"),
                                   Healthy  = median(Kuster_Pathway_Cor_diff[[.x]]),
                                   Cancer = median(CCLE_Pathway_Cor_diff[[.x]])) %>%
                           mutate(In_Cancer = Cancer-Healthy)) %>% left_join(KEGG_pathways) ),
    tar_target(
        healthy_cancer_corr_pathways_plot,
        healthy_cancer_corr_pathways %>%
            subset(Path_type == "metabolic") %>%
            arrange(p_value) %>%
            mutate(Path_id = factor(Path_id,Path_id)) %$%
            ggplot(.,aes(x = Path_id,y = -log10(p_value+0.00000000000000000001), label = Path_description, colour = In_Cancer,))+
            scale_color_distiller(type = "div", limit = max(abs(healthy_cancer_corr_pathways$In_Cancer)) * c(-1, 1))+
            geom_point()+
            ggrepel::geom_text_repel(data = . %>% head(7))+
            ggtitle("Distribution difference of mRNA Prot corr between healthy and cancer",
                    subtitle = "Ks.test of each pathway difference to median, add colour if up or down")),
    
    tar_target(
        PTR_mean_df_plot,
        Mean_sd_matrix(PTR_CCLE) %>%
            ggplot(aes(x = Mean, y = SD)) +
            geom_point(size = 0.1,alpha = 0.2)+ 
            xlab("Mean across Cell lines") +
            ylab("SD across Cell lines")+
            ggtitle("PTR_df_plot")),
    # tar_target(
    #     eTFs_PTR_corr,
    #     Factors_to_PTR_corr(uniprot_factors %>% subset(str_detect(Type,"eukaryotic_translation")) %>%
    #                             pull(Entry),PTR_CCLE, CCLE_proteins)),
    # tar_target(
    #     Ridge_eIFS_prediction,
    #     Running_Ridge_eIFs(CCLE_proteins,PTR_CCLE,uniprot_factors)),
    tar_target(
        PTR_CCLE_Tissue,
        PTR_CCLE%>% t() %>% as.data.frame() %>%
            rownames_to_column("Cell_line") %>%
            left_join(sample_info[,c("stripped_cell_line_name","lineage")],by = c("Cell_line" = "stripped_cell_line_name")) %>%
            dplyr::select(-Cell_line) %>%
            pivot_longer(cols =-lineage, names_to = "GeneName", values_to= "Abundance") %>%
            add_count(GeneName,lineage)%>%
            subset(n>5)%>%
            group_by(GeneName,lineage) %>% 
            summarise(Avg_PTR = median(Abundance,na.rm = T)) %>%
            ungroup()  %>% pivot_wider(names_from = "lineage", values_from = "Avg_PTR") %>%
            column_to_rownames("GeneName") %>% as.matrix()),
    tar_target(
        pathway_PTR_CCLE,
        Calc_path_PTR_Tissues(PTR_CCLE_Tissue , KEGG_genes,KEGG_pathways)),
    tar_target(
        PTR_metabo_pathways_CCLE_plot,
        ggplot(pathway_PTR_CCLE,aes(x = value, y = Tissue,fill = Tissue))+
            ggridges::geom_density_ridges(rel_min_height = 0.01,alpha = 0.5)+
            geom_vline( xintercept = 0) +
            ggtitle("Metabolism PTR Across Cancer tissues", 
                    subtitle = "subsetted for more that 15 genes per pathway")+
            facet_wrap("Path_description") + theme_bw()+
            theme(axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())),
    tar_target(
        PTR_metabo_pathways_Kuster_plot,
        ggplot(pathway_PTR_Kuster,aes(x = value, y = Tissue,fill = Tissue))+
            ggridges::geom_density_ridges(rel_min_height = 0.01,alpha = 0.5)+
            geom_vline( xintercept = 5) +
            ggtitle("Metabolism PTR Across Healthy tissues", 
                    subtitle = "subsetted for more that 15 genes per pathway")+
            facet_wrap("Path_description") + theme_bw()+
            theme(axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())),
    tar_target(
        complexes, 
        import_omnipath_complexes(resources=c("CORUM")) %>%
            mutate(name = janitor::make_clean_names(name))),
    tar_target(
        complexes_units,
        c(Proteosome = "PSMA1",
                        Initiation1 = "EIF1",
                        Initiation2 = "EIF2",
                        Initiation3 = "EIF3",
                        #Initiation4 = "EIF4",
                        #Elongation = "EEF1",
                        Splicing = "ACIN1",
                        Ribosome_60 = "RPL10",
                        Ribosme_40 = "RPS10",
                        Tubulin = "TUBG1")),
    tar_target(
        Complex_complete_units,
        purrr::map(.x = complexes_units,~Extract_complex_uniprots(.x,complexes,complexes_units)) %>%
            unlist(recursive = F)),
    tar_target(
        Kuster_Complex_PCs,
        purrr::imap_dfc(.x=Complex_complete_units, ~Get_PC_complexes(.x,.y,Protein_Kuster))),
    tar_target(
        CCLE_Complex_PCs,
        purrr::imap_dfc(.x=Complex_complete_units, ~Get_PC_complexes(.x,.y,CCLE_RNA_prot_aligned[[1]]))),
    tar_target(
        Kuster_PC_corr,
        pheatmap::pheatmap(cor(Kuster_Complex_PCs), show_colnames = F)),
    tar_target(
        CCLE_PC_corr,
        pheatmap::pheatmap(cor(CCLE_Complex_PCs), show_colnames = F)),
    tar_target(
        pathway_PTR_Kuster,
        Calc_path_PTR_Tissues(PTR_Kuster, KEGG_genes,KEGG_pathways)),
    tar_target(
        PTR_Kuster_ridge_proteins,
        Running_Ridge_eIFs(Protein_Kuster,PTR_Kuster,Kuster_Complex_PCs)),
    tar_render(Presentation, fs::path_rel(here::here("Output","Presentation.Rmd"))),
    tar_render(report, fs::path_rel(here::here("Output","report.Rmd")))
)
