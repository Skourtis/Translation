library(targets)
library(openxlsx)
library(tarchetypes)
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
        sample_info,
        read.csv(sample_info_file, stringsAsFactors = FALSE)),
    tar_target(
        Metab_CCLE_file, ####Raw from https://portals.broadinstitute.org/ccle/data CCLE_metabolomics_20190502.csv 2/11/2020
        here::here("Datasets","Raw","CCLE_metabolites_landscape_of_cancer.xlsx"),
        format = "file"    ),
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
        CCLE_prot,
        loading_CCLE_prot(x= CCLE_prot_file)),
    tar_target(
        CCLE_RNA_seq_file,#Raw from https://depmap.org/portal/download/ 2/11/2020
        here::here("Datasets","Raw","CCLE_expression.csv"),
        format = "file"),
    tar_target(
        CCLE_RNA_seq,
        inner_join(sample_info[,1:2],read.csv(CCLE_RNA_seq_file), by = c("DepMap_ID" = "X"))[,-1] %>%
            column_to_rownames(var = "stripped_cell_line_name") %>% t() %>%
            magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2])),    
    tar_target(
        NCI_60_metabolites_file,#Raw from  https://www.nature.com/articles/s41467-019-09695-9#Sec28 2/11/2020
        here::here("Datasets","Raw","41467_2019_9695_MOESM2_ESM.xlsx"),
        format = "file"),
    tar_target(
        NCI_60_metabolites,
        openxlsx::read.xlsx(NCI_60_metabolites_file, sheet = 3, startRow = 4, na.strings = "NaN") %>%
            .[str_detect(.$`Annotation.ID`,"H|C"),-c(1:3,5)] %>% setNames(str_remove_all(colnames(.), "^[:graph:]*?_"))  %>%
            na.omit() ),
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
            setNames(str_match(colnames(.), "^([:graph:]*?)_")[,2])    )
    #tar_render(report, here::here("Output","report.Rmd")),

 )
