pacman::p_load("biglm",
               "readxl",
               "rmarkdown",
               "tarchetypes",
               "targets",
               "tidyverse")

#' @title Load proteomics data
#' @description Import and clean proteomic CCLE Data
#' @return A Dataframe of TenP protein data, removed highly absent values and isoforms trimmed to master Uniprot
#' @param ##--
#' @examples #x <- data.frame(HEL = runif(5,0,1),MDAMB231 = runif(5,0,1), T47D = runif(5,0,1), Uniprot_Acc = c("O95758-4","P78368","O75145","Q99996","Q99996-2")) # gene_names
loading_CCLE_prot <- function(x){
    CCLE_proteins <- openxlsx::read.xlsx(tar_read(CCLE_prot_file), sheet =2) %>%
        .[,!str_detect(colnames(.),"Peptides|Column")]
    colnames(CCLE_proteins) <- str_remove_all(colnames(CCLE_proteins),"_TenPx..")
    colnames(CCLE_proteins) <- str_match(colnames(CCLE_proteins),"^([:graph:]*?)_")[,2]
    CCLE_proteins <- CCLE_proteins[,which(!duplicated(colnames(CCLE_proteins)))] ##Removing cell_lines with the same name
    CCLE_proteins$Uniprot <- str_remove_all(CCLE_proteins$Uniprot,"-.")
    CCLE_proteins <- CCLE_proteins %>% subset(!duplicated(Uniprot))
    rownames(CCLE_proteins) <- NULL
    CCLE_proteins <- CCLE_proteins %>%
        column_to_rownames(var = "Uniprot") %>%
        dplyr::select(-c(1:5)) 
    index_to_keep <- CCLE_proteins %>% is.na() %>%
    rowSums()<(ncol(CCLE_proteins)*0.75)
    CCLE_proteins %>% .[index_to_keep,] # removing lines with more than 75% NA]
}

#' @title Aligning_matrices
#' @description Take two matrices, Subset for common columns/Rows, reorder and return
#' @return Return a named list of the two matrices, in the same order they were given,
#' @param 2 Matrices
#' @examples  Matrix_1 <- matrix(1:20,nrow = 4,ncol=5, dimnames = list(letters[1:4],LETTERS[1:5]))
#' @examples  Matrix_2 <- matrix(1:30,nrow = 5,ncol=6, dimnames = list(letters[2:6],LETTERS[3:8]))

Aligning_two_matrices <- function(Matrix_1,Matrix_2){
    Names_1 <- deparse(substitute(Matrix_1))
    Names_2 <- deparse(substitute(Matrix_2))
    Common_cols <- intersect(colnames(Matrix_1),colnames(Matrix_2)) #Common Samples
    Common_rows <- intersect(rownames(Matrix_1),rownames(Matrix_2)) #Common Features
    list(as.matrix(Matrix_1) %>%  .[Common_rows,Common_cols],
         as.matrix(Matrix_2) %>% .[Common_rows,Common_cols]) %>%
        set_names(c(Names_1,Names_2))
}

PTR_aligning_matrices <- function(Matrix_1,Matrix_2){
    
    PTR_list <- Aligning_two_matrices(Matrix_1,Matrix_2)
    

}

#' @title Z transformation
#' @description Plot a histogram of ozone concentration.
#' @return A ggplot histogram showing ozone content.
#' @param data 
#' @examples
####
Z_transf <- function(x){
    x <- as.matrix(x)
    (x - mean(x, na.rm = T))/(sd(x, na.rm = T))
}
#' @title Z transformation
#' @description Plot a histogram of ozone concentration.
#' @return A ggplot histogram showing ozone content.
#' @param data 
#' @examples
#pairwise correlation function used in main function
Correlation_pairwise <- function(x,y){
    cor(x,y,use = "pairwise", method = "pearson")
}
#' @title Z transformation
#' @description Plot a histogram of ozone concentration.
#' @return A ggplot histogram showing ozone content.
#' @param data 
#' @examples
# Get lower triangle of the correlation matrix
get_upper_tri<-function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
}
# Get upper triangle of the correlation matrix
#' @title Z transformation
#' @description Plot a histogram of ozone concentration.
#' @return A ggplot histogram showing ozone content.
#' @param data 
#' @examples
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

#' @title Intra_omic Feature correlation
#' @description Creates correlation matrix for features in a single omic
#' @return Symmetrical Correlation Matrix with identical row and col names
#' @param Omic optional(Genes/Proteins of interest) 
#' @examples Omic <- data.frame(HEL = runif(5,0,1),MDAMB231 = runif(5,0,1), T47D = runif(5,0,1), row.names = c("RANBP6","UNC119","HES4")) # gene_names

#compare 2 datasets with common Feature (rows) and Sample names (Columns)
### Change rerun value
Intra_omic_Feature_correlation <- function(Omic = NULL,Genes = NULL){
    #Omic <- data.frame(HEL = runif(3,0,1),MDAMB231 = runif(3,0,1), T47D = runif(3,0,1), row.names = c("RANBP6","UNC119","HES4"))
    if(is.null(Genes)){
        Genes <-  rownames(Omic)
    }
    
    Omic <- Omic %>%
        subset(row.names(.) %in% Genes)
    cor(t(Omic), use = "pairwise.complete.obs") %>%
        round(2)
    
}

###testing function
produce_data <- function() {
    expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
}

#' @title Mean_sd_matrix
#' @description Calculating mean and sd of a matrix row_wise
#' @return Long Data.frame with columns "Gene", "Mean", "SD"
#' @param Data_matrix 
#' @examples Data_matrix <- data.frame(RANBP6 = c(1,0.16,-1),UNC119 = c(0.16,1,-0.25), HES4 = c(-1,0.25,1), row.names = c("RANBP6","UNC119","HES4")) %>% as.matrix()
Mean_sd_matrix <- function(Data_matrix = NULL){
    Data_matrix <- Data_matrix %>% as.matrix()
    data.frame(Genes = row.names(Data_matrix),
               Mean = Data_matrix%>% rowMeans(na.rm = T),
               SD = Data_matrix%>% matrixStats::rowSds(na.rm = T))
}


#' @title df_to_Heatmap
#' @description Creates heatmap from matrix with column and row names
#' @return Heatmap_object
#' @param Data.frame with other optionals optional(Genes/Proteins of interest) 
#' @examples Data_df <- data.frame(RANBP6 = c(1,0.16,-1),UNC119 = c(0.16,1,-0.25), HES4 = c(-1,0.25,1), row.names = c("RANBP6","UNC119","HES4"))
#' @examples Main_top_annotation <- data.frame(Genes = c("UNC119", "HES4"), Pathways = c("b","b"))
#' @examples Subset_col / Subset_row = "RANBP6"

Convert_df_to_heatmap <- function(Data_df = NULL, Main_top_annotation = NULL,#2 column dataframe with colnames Genes and Pathways
                                  Heatmap_name = "Give me a name",
                                  SubAnnotation_top = NULL, legend = c(TRUE,FALSE),
                                  Main_side_annotation = NULL, SubAnnotation_side = NULL, 
                                  Subset_col = NULL, #vector of features/ samples to keep
                                  Subset_row = NULL, 
                                  cluster = c("both",'columns',"rows","none"),
                                  show_names = c("both",'columns',"rows","none"),
                                  row_order = NULL, column_order = NULL,
                                  ...){ ### arguments to pass to heatmap
    # Data_df <- CCLE_prot_Intra_omics_corr_matrix
    # Heatmap_name = "eTFs Proteome against all proteins"
    # Subset_col = NULL #vector of features/ samples to keep
    # Subset_row = tar_read(uniprot_factors) %>% subset(Type == "eukaryotic_translation") %>%
    #     pull(Entry)
    # show_names = "rows"
    # Main_top_annotation = tar_read(gene_clusters_eIFs) %>% dplyr::select(Uniprot, Pathways) %>%
    #     rename(Genes = Uniprot)
    if(is.null(Subset_col)){
        Subset_col <- colnames(Data_df)
    }
    if(is.null(Subset_row)){
        Subset_row <- rownames(Data_df)
    }
    
    Data_df <- Data_df %>% .[rownames(.) %in% Subset_row, 
                             colnames(.) %in% Subset_col] %>%
        as.matrix()
    
    set.seed(30)# to set random generator seed
    cl <- colors(distinct = TRUE)
    Creating_annotations <- function(Annotation,name){
        #Annotation <- Main_top_annotation #########testing
        if(is.null(Annotation)){
            #return(1)
            return(list(NULL,NULL))
        } else if (str_detect(name,"top")) {
            joined <- left_join(as.data.frame(colnames(Data_df)) %>% set_names("Top"), 
                                Annotation, by = c("Top" = "Genes"))
            joined[is.na(joined)] <- "Not in List"
            return(list(Pathways = joined %>% pull(Pathways),
                        Path_colours = sample(cl, length(unique(joined %>% pull(Pathways)))) %>%
                            setNames(unique(joined %>% pull(Pathways)))))
        }else if (str_detect(name,"side")){
            joined <- left_join(as.data.frame(rownames(Data_df)) %>% set_names("Side"), 
                                Annotation, by = c("Side" = "Genes"))
            joined[is.na(joined)] <- "Not in List"
            return(list(Pathways = joined %>% pull(Pathways),
                        Path_colours = sample(cl, length(unique(joined %>% pull(Pathways)))) %>%
                            setNames(unique(joined %>% pull(Pathways)))))
        }
    }
    Annotations <- purrr::imap(.x= list("Main_top_annotation" = Main_top_annotation,
                                        "SubAnnotation_top" = SubAnnotation_top,
                                        "Main_side_annotation" = Main_side_annotation,
                                        "SubAnnotation_side" = SubAnnotation_side),
                               ~Creating_annotations(.x,.y))
    
    if(is.null(Main_top_annotation) & is.null(SubAnnotation_top)){
        ha = NULL
    }else{
        ha = HeatmapAnnotation(Main_top = Annotations[[1]][[1]],
                               Sub_top = Annotations[[2]][[1]],
                               col = list(Main_top = Annotations[[1]][[2]],
                                          Sub_top = Annotations[[2]][[2]]))
    }
    if(is.null(Main_side_annotation) & is.null(SubAnnotation_side)){
        side = NULL
    }else{
        side = rowAnnotation(Main_side = Annotations[[3]][[1]],
                             Sub_side = Annotations[[4]][[1]],
                             col = list(Main_side = Annotations[[3]][[2]],
                                        Sub_side = Annotations[[4]][[2]]))
    }
    
    
    Heatmap(Data_df, name = Heatmap_name, 
            top_annotation = ha  ,
            left_annotation = side,
            cluster_rows = if_else(cluster %in% c("rows","both"),T,F), 
            cluster_columns = if_else(cluster %in% c("columns","both"),T,F),
            show_column_names  = if_else(show_names %in% c("columns","both"),T,F),
            show_row_names = if_else(show_names %in% c("rows","both"),T,F), 
            #row_title = paste0("top = ", deparse(substitute(Dataset_1))), 
            #column_title = paste0("bottom = ", deparse(substitute(Dataset_2))),
            column_title_side = "bottom", row_names_gp = gpar(fontsize = 8),
            heatmap_width  = unit(20,"cm"),
            heatmap_height = unit(20,"cm"),
            raster_quality = 1,
            clustering_distance_columns = "spearman",
            clustering_distance_rows = "spearman",
            row_order = row_order, column_order = column_order)
    
    
}









Dataset_correlation <- function(Dataset_1,Dataset_2,pathways){
    Comparison_type <- paste(deparse(substitute(Dataset_1)),deparse(substitute(Dataset_2)))
    Dataset1 <- as.data.frame(Dataset_1)
    Dataset2 <- as.data.frame(Dataset_2)
    Common_Cell_lines <- intersect(colnames(Dataset1),colnames(Dataset2)) #Common Samples
    Common_proteins <- intersect(rownames(Dataset1),rownames(Dataset2)) #Common Features
    
    Subset_and_matrix <- function(x, triange = c("top","bottom")){
        #x <- Dataset1
        #subsetting and ordering the Datasets
        x_sub <- x %>% .[rownames(.) %in% Common_proteins, colnames(.) %in% Common_Cell_lines]%>% 
            .[match(Common_proteins,rownames(.)),match(Common_Cell_lines,colnames(.))]
        #print(deparse(substitute(x)))
        whole_matrix <- round(cor(t(x_sub[rownames(x_sub) %in% Master_pathways_ordered$ID,])),2) %>% 
            .[match(Master_pathways_ordered$ID,row.names(.)),match(Master_pathways_ordered$ID,colnames(.))]
        cormat <- whole_matrix  %>%
            {if(triange == "top")  get_upper_tri(.)else get_lower_tri(.) }
        cormat[is.na(cormat)] <- 0  
        list(whole_matrix = whole_matrix,
             matrix = cormat,
             df = x_sub)
        #cormat
    }
    
    Master_pathways_ordered <- pathways %>% subset(ID %in% Common_proteins) %>%
        subset(!duplicated(ID))%>% arrange(Pathway,SubPathway)
    ##subsetting for specific master pathways
    # Master_pathways <- Master_pathways %>% 
    #   subset(Pathway %in% c("Amino acid","Nucleotide","TCA cycle"))
    
    #left_join(Gene_pathways[,c("Uniprot","SubPathway")], by = c("ID" = "Uniprot")) %>% subset(!duplicated(ID))
    Master_pathways_ordered[is.na(Master_pathways_ordered)] <- "Not Known"
    
    
    
    #test that Dataset1[1,1] =! Dataset2 [1,1]
    list1 <- Subset_and_matrix(Dataset_1,"top")
    list2 <- Subset_and_matrix(Dataset_2, "bottom")
    Dataset1_sub <- list1$df %>% as.data.frame()
    Dataset2_sub <- list2$df %>% as.data.frame()
    #Combined_df <- rbind(cormat1,cormat2)  
    # #Combined_df$value[Combined_df$Var1 == Combined_df$Var2] <- NA
    # 
    # Feature_order <- Combined_df[!duplicated(Combined_df$Var1),] %>% 
    #   .[order(.$Pathway),"Var1"]
    # 
    # ### ComplexHeatMAP
    # test1 <- Master_pathways_ordered[order(Master_pathways_orderezd$Uniprot),]
    Metabolic_matrix <- as.data.frame(list1$matrix+list2$matrix) %>% as.matrix()
    
    
    set.seed(1587) # to set random generator seed
    cl <- colors(distinct = TRUE)
    
    
    Master_Pathway_colours <- sample(cl, length(unique(Master_pathways_ordered$Pathway))) %>% setNames(unique(Master_pathways_ordered$Pathway))
    Sub_Pathway_colours <- sample(cl, length(unique(Master_pathways_ordered$SubPathway))) %>% setNames(unique(Master_pathways_ordered$SubPathway))
    Metabolic_matrix[Metabolic_matrix ==2] <- 1
    ha = HeatmapAnnotation(Master_path = Master_pathways_ordered$Pathway,
                           #Sub_path = Master_pathways_ordered$SubPathway,
                           col = list(Master_path = Master_Pathway_colours#,Sub_path = Sub_Pathway_colours
                           ))
    side = rowAnnotation(Master_path = Master_pathways_ordered$Pathway,
                         #Sub_path = Master_pathways_ordered$SubPathway,
                         col = list(Master_path = Master_Pathway_colours#,Sub_path = Sub_Pathway_colours
                                    ))
    
    hmap <- Heatmap(Metabolic_matrix, name = Comparison_type, 
                    top_annotation = ha ,
                    left_annotation = side,
                    cluster_rows = F, cluster_columns = F, show_column_names  = F, show_row_names = F, 
                    row_title = paste0("top = ", deparse(substitute(Dataset_1))), 
                    column_title = paste0("bottom = ", deparse(substitute(Dataset_2))),
                    column_title_side = "bottom",
                    heatmap_width  = unit(20,"cm"),
                    heatmap_height = unit(20,"cm"),
                    raster_quality = 1)
    
    # ### GGPLOT
    # Combined_df %>%
    #     subset(Var1 %in% Master_pathways$Uniprot & Var2 %in% Master_pathways$Uniprot) %>% na.omit() %>%
    #     ggplot(aes(factor(Var2, levels = Feature_order), factor(Var1,levels = Feature_order), fill = value))+
    #     geom_tile(color = "white")+
    #     scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "black",
    #                          midpoint = 0, limit = c(-1,1), space = "Lab", 
    #                          name="Pearson\nCorrelation") +
    #     theme_minimal()+ 
    #     theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
    #     coord_fixed() + labs(x = deparse(substitute(Dataset_1)),deparse(substitute(Dataset_2)))#+
    # facet_grid(Pathway~., scales = "free_x")
    
    
    #Correlating corresponding Samples Correlation of only paired values which allows for NA values
    Cell_line_correlation <- map2_dbl(
        Dataset1_sub, 
        Dataset2_sub, 
        Correlation_pairwise) 
    
    Cell_line_perm <- purrr::rerun(100, 
                                   map2_dbl(
                                       Dataset1_sub[,sample(ncol(Dataset1_sub))],#Shuffle column-wise and performing sample correlation as above
                                       Dataset2_sub,
                                       Correlation_pairwise)) %>% unlist()
    
    gene_correlation <- pmap(
        list(pmap(Dataset1_sub, c, use.names = T), # unlisting feature rows into vectors into lists and correlating
             pmap(Dataset2_sub, c, use.names = T)),
        Correlation_pairwise) %>% unlist() %>% setNames(Common_proteins) 
    
    gene_correlation_perm <- purrr::rerun(100, #permutation by row shuflfing and correlation
                                          pmap(list(pmap(Dataset1_sub[sample(nrow(Dataset1_sub)),], 
                                                         c, use.names = T), #Shuffle row-wise
                                                    pmap(Dataset2_sub, c, use.names = T)),
                                               Correlation_pairwise) %>% unlist()) %>% unlist()
    
    # Imputing both datasets with Media to derive overall correlation
    Dataset1_sub_imp <- naniar::impute_median_all(Dataset1_sub)
    Dataset2_sub_imp <- naniar::impute_median_all(Dataset2_sub)
    
    Matri_Corr <- MatrixCorrelation::r1(Dataset1_sub_imp,Dataset2_sub_imp) # whole matrix correlation
    
    Matrix_perm <- purrr::rerun(100, #whole matrix permutation of Dataset 1 col-wise shuffle and Dataset 2 row-wise shuffle and correlation
                                MatrixCorrelation::r1(Dataset1_sub_imp[,sample(ncol(Dataset1_sub_imp))],
                                                      Dataset2_sub_imp[sample(nrow(Dataset2_sub_imp)),])) %>% 
        unlist() # permutation
    
    Hist_plot <- function(Data,Perm){
        #Data = Cell_line_correlation
        #Perm = Cell_line_perm
        Sample_Corr <- data.frame(Samples = c(Data,unlist(Perm)),
                                  Type = rep(c("Data","Perm"), times = c(length(Data),length(unlist(Perm)))))
        p <- ggplot2::ggplot(Sample_Corr,aes(Samples)) +
            geom_histogram(aes(y=..density..,fill = Type ),      # Histogram with density instead of count on y-axis
                           binwidth=.01,alpha=0.1) +
            geom_density(alpha=.5, aes(colour = Type))+
            ggtitle(paste(Comparison_type,deparse(substitute(Data))))
    }
    
    Feature_quartiles <- rowMeans(Dataset2_sub, na.rm = T) %>% data.frame(Feature = names(.),
                                                                          Mean_Abundance = .) %>%
        mutate(quartile = as.factor(ntile(Mean_Abundance, 10))) %>% left_join(stack(gene_correlation), by = c("Feature" = "ind"))
    
    Feature_expression_cor_plot <- ggplot(Feature_quartiles, aes(values,after_stat(density), colour =quartile)) +
        geom_freqpoly(binwidth = 0.1) +
        ggtitle(paste("Feature Expression",Comparison_type))
    
        list(Sample_Corr = list(Cell_line_perm = Cell_line_perm, 
                            Histogram = Hist_plot(Cell_line_correlation,Cell_line_perm)),
         Feature_Corr = list(Features_quart = Feature_quartiles, 
                             Gene_corr_perm = gene_correlation_perm, 
                             Histogram = Hist_plot(gene_correlation, gene_correlation_perm), 
                             Feature_cor_plot = Feature_expression_cor_plot),
         Matrix_Corr= list(Matri_Corr = Matri_Corr, 
                           Matrix_perm = Matrix_perm,
                           Histogram =   Hist_plot(Matri_Corr,Matrix_perm)),
         Heatmap= list(Dataset1_matri = list1$whole_matrix,
                       Dataset2_matri = list2$whole_matrix,
                       Plot_Heatmap = hmap))
    
}



#' @title Graphite pathway Extraction
#' @description Retrieves different databases from Graphite package
#' @return A dataframe
#' @param Database_name,species 
#' @examples species "hsapiens", Database = "reactome"
# Graphite_extract <- function(species = "hsapiens",
#                              Database = c("kegg","reactome")){
#     pacman::p_load(graphite)
#     Pathways <- graphite::pathways(species, Database)
#     df <- data.frame(SubPathway = NULL,
#                      ID = NULL)
#     for(i in names(Pathways)){
#         #print(i)
#         if(length(nodes(Pathways[[i]]))>0){
#             df <- rbind(df,data.frame(SubPathway = i,
#                                                           ID = nodes(Pathways[[i]])))
#         }
#         
#     }
#     if(Database == 'kegg'){
#         df <-    df %>% mutate(ID = str_remove_all(ID,"ENTREZID:")) %>% 
#             left_join(subset(tar_read(HUMAN_9606_idmapping),Type == "GeneID") %>%
#                           dplyr::select(-Type),by = "ID") %>% 
#             left_join(subset(tar_read(HUMAN_9606_idmapping),Type == "Gene_Name") %>%
#                           dplyr::select(-Type),by = "Uniprot") %>%
#             dplyr::select(-ID.x) %>% na.omit()%>% 
#             setNames(c("Kegg_Pathway","Uniprot","Gene_name"))%>%
#             dplyr::add_count(Kegg_Pathway) %>% 
#             subset(n>15) %>% dplyr::select(-n)
#     }
#     
# df
# }
# Graphite_extract(species = "hsapiens",
#                  Database = "kegg")
Retrieve_all_kegg_genes <- function(pathways){
    genes_reactions <- data.frame(Gene_id = NULL,
                                  pathway = NULL,
                                  stringsAsFactors = F)
    Retrieve_Genes_KEGG <- function(pathway, Gene_df =  genes_reactions){
        tmp <- tempfile()
        #pathway = "hsa05416"
        #Gene_df = genes_reactions
        #  Compound_df = compound_reactions
        KEGGgraph::retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
        KGMLGraph <- KEGGgraph::parseKGML2Graph(tmp)
        KGML_n <- KEGGgraph::parseKGML(tmp)
        # mapkG2 <- KEGGpathway2Graph(KGML_n, expandGenes=TRUE)
        # set.seed(124)
        # randomNodes <- sample(nodes(mapkG2), 25)
        # mapkGsub <- subGraph(randomNodes, mapkG2)
        # plot(mapkGsub)
        
        if(length(KGMLGraph@nodeData@defaults[["KEGGNode"]][["nodes"]])>0){
            for (i in 1:length(KGMLGraph@nodeData@defaults[["KEGGNode"]][["nodes"]])){
                df <- data.frame(Gene_id = KGMLGraph@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@entryID,
                                 pathway = pathway,
                                 stringsAsFactors = FALSE)
                Gene_df <- rbind(Gene_df,df)
                
            }
        }
        Gene_df
    }
    Retrieve_Genes_KEGG_safe <- safely(Retrieve_Genes_KEGG)
    
    Gene_pathways <- map(pathways ,Retrieve_Genes_KEGG_safe) %>%
        map("result")%>%
        compact() %>% purrr::reduce(rbind) %>% 
        left_join(subset(tar_read(HUMAN_9606_idmapping),Type == "KEGG") %>%
                      dplyr::select(-Type),by = c("Gene_id" = "ID")) %>%
        dplyr::select(-Gene_id)%>%
        left_join(subset(tar_read(HUMAN_9606_idmapping),Type == "Gene_Name") %>%
                      dplyr::select(-Type),by = "Uniprot") %>% na.omit()
        
    
}
        
#' @title Factors_to_PTR_corr
#' @description Calculates pairwise correlation between all common cell lines between PTRs and Protein abundance
#' @return Long dataframe with the Factor, Gene PTR and Correlation
#' @param Data.frame with other optionals optional(Genes/Proteins of interest) 
#' @examples PTRs_df <- PTR_CCLE [1:10,]
#' @examples Proteins_df <- CCLE_proteins
#' @examples Factors_of_interest <- c("O00570","Q9Y4B5","Fake_uniprot")

Factors_to_PTR_corr <- function(Factors_of_interest,PTRs_df, Proteins_df){
    Protein_df <- Proteins_df[rownames(Proteins_df) %in% Factors_of_interest,] %>%
        rownames_to_column("Uniprot") %>% group_split(Uniprot, .keep = T)
    Protein_df <- Protein_df %>% set_names(.,map_chr(.x = .,~.x[["Uniprot"]])) 
    PTRs <-  PTRs_df %>% as.data.frame() %>%
        rownames_to_column("Uniprot") %>% group_split(Uniprot, .keep = T)
    PTRs <- PTRs %>% set_names(.,map_chr(.x = .,~.x[["Uniprot"]])) 
    
    eIFCorrProtein <- function(x,y){
        #x = map(Combinations$Uniprot_PTR[1:2], ~PTRs[[.x]])[[1]]
        #y = map(Combinations$Factor[1:2], ~Protein_df[[.x]])[[1]]
        common_lines <- intersect(colnames(y)[-1], colnames(x)[-1])
        
        correlation <- cor(x[,common_lines] %>% unlist(),
                           y[,common_lines]%>% unlist(), use = "pairwise.complete.obs")
        
        data.frame(Uniprot_PTR = unique(x$Uniprot), Factor = unique(y$Uniprot), Correlation = correlation)
        #browser()
    }
    Combinations <- Combinations <- tidyr::crossing(names(Protein_df),
                                                    names(PTRs)) %>%
        set_names(c("Factor","Uniprot_PTR"))
    eIFCorrProtein_safe <- safely(eIFCorrProtein)
    library(future)
    plan(multisession, workers = future::availableCores())
    eIF_correlations_f <-  future::future({map2(
        map(Combinations$Uniprot_PTR, ~PTRs[[.x]]),
        map(Combinations$Factor, ~Protein_df[[.x]]),eIFCorrProtein_safe)},
        seed = T)
    eIF_correlations <- value(eIF_correlations_f) %>%
        map("result") %>%  compact()%>% 
        purrr::reduce(rbind)
    plan(sequential)
    eIF_correlations
}
    


