Quality_Control <- function(x){
    ### takes as input the txt output folder of MaxQuant 
    ### in the form of the here::here function
    ### and produces a QC report
    pacman::p_temp("PTXQC")
    require(methods)
    r = createReport(x)
}

load_MaxQuant <- function(txt_folder){
    ### This function takes the txt folder, 
    ### the name of the measure_col,number of samples and 
    ### names of the conditions and output the peptides per condition and 
    ### reads the Maxquant proteinGroups file
    #txt_folder <-  here::here("Datasets", "Raw",
     #                         "2020MQ044txt",
      #                        "txt")
  
    condition_names <-  c(rep("control",3),
                          rep("sh_1",3),
                          rep("sh_2",3))
    samples_ids <-  1:9
    measure_col <-  "Intensity"
    
    measCols <- list(
        HL = measure_col
    )
    
    eviCols <- list(
        sequence = 'Sequence',
        modified_sequence = 'Modified sequence',
        modifications = 'Modifications',
        protein_group = 'Proteins',
        protein = 'Leading razor protein',
        experiment = 'Experiment',
        charge = 'Charge',
        reverse = 'Reverse',
        contaminant = 'Potential contaminant'
    )
    evi <- proteus::readEvidenceFile(paste0(txt_folder,
                                           "/evidence.txt"), 
                                     measure.cols=measCols, 
                                     data.cols=eviCols, 
                                     zeroes.are.missing=FALSE)
    
    
    
    meta <- data.frame(experiment = sort(evi$experiment %>% unique()),
                       measure = measure_col,
                       condition = condition_names,
                       sample = sort(evi$experiment %>% unique()),
                       replicate = c(rep(1,3),rep(2,3),rep(3,3)))
    
    pepdat <- proteus::makePeptideTable(evi, 
                                       meta,
                                       measure.cols=measCols, 
                                       experiment.type="label-free")
    
    measure_columns <- list(paste(measure_col,samples_ids, sep = " " )) %>% unlist %>%
        setNames(paste("Sample", samples_ids, sep = "_"))

    list(proteins = proteus::readProteinGroups(paste0(txt_folder,
                            "/proteinGroups.txt"), 
                      meta,
                      measure.cols = measure_columns),
         peptides = pepdat, meta, measure.cols = measure_columns
         )
}
load_MaxQuanttest <- function(txt_folder){
  ### This function takes the txt folder, 
  ### the name of the measure_col,number of samples and 
  ### names of the conditions and output the peptides per condition and 
  ### reads the Maxquant proteinGroups file
  #txt_folder <-  here::here("Datasets", "Raw",
  #                         "2020MQ044txt",
  #                        "txt")
  
  condition_names <-  c(1:9)
  samples_ids <-  1:9
  measure_col <-  "Intensity"
  
  measCols <- list(
    HL = measure_col
  )
  
  eviCols <- list(
    sequence = 'Sequence',
    modified_sequence = 'Modified sequence',
    modifications = 'Modifications',
    protein_group = 'Proteins',
    protein = 'Leading razor protein',
    experiment = 'Experiment',
    charge = 'Charge',
    reverse = 'Reverse',
    contaminant = 'Potential contaminant'
  )
  evi <- proteus::readEvidenceFile(paste0(txt_folder,
                                          "/evidence.txt"), 
                                   measure.cols=measCols, 
                                   data.cols=eviCols, 
                                   zeroes.are.missing=FALSE)
  
  
  
  meta <- data.frame(experiment = sort(evi$experiment %>% unique()),
                     measure = measure_col,
                     condition = condition_names,
                     sample = sort(evi$experiment %>% unique()),
                     replicate = 1,1,1,2,2,2,3,3,3)
  
  pepdat <- proteus::makePeptideTable(evi, 
                                      meta,
                                      measure.cols=measCols, 
                                      experiment.type="label-free")
  
  measure_columns <- list(paste(measure_col,samples_ids, sep = " " )) %>% unlist %>%
    setNames(paste("Sample", samples_ids, sep = "_"))
  
  list(proteins = proteus::readProteinGroups(paste0(txt_folder,
                                                    "/proteinGroups.txt"), 
                                             meta,
                                             measure.cols = measure_columns),
       peptides = pepdat, meta, measure.cols = measure_columns
  )
}


is_function = function (expr) {
    if (! is_assign(expr))
        return(FALSE)
    value = expr[[3]]
    is.call(value) && as.character(value[[1]]) == 'function'
}

function_name = function (expr){
    as.character(expr[[2]])}

is_assign = function (expr){
    is.call(expr) && as.character(expr[[1]]) %in% c('=', '<-', 'assign')}


produce_ipath_map <- function(df,
                              column_name = NULL){
    ### takes as input a dataframe produced which converts ratios to colours
    ### and the columns to map and produces iPath3 images 
    #df = For_ipath_comparisons
    #column_name = columns_for_mapping[1]
    
    Filtered_proteins <- df %>% 
        dplyr::select(all_of(c("V1","Width",column_name))) %>% 
        na.omit()
    
    #Filtered_proteins$Width <- paste("W",as.character(ntile(Filtered_proteins[,i], 30)),sep="")
    Filtered_proteins$Ipath <- paste(Filtered_proteins$V1, 
                                     Filtered_proteins %>% pull(column_name), 
                                     Filtered_proteins$Width)
    
    file_name <- here::here("Project_Output",
                      paste0(column_name,
                             "_ipath_all_enzymes_bins.tsv"))
    
    selections <- paste('selection=',
                        paste(paste(Filtered_proteins$Ipath,
                                    "%0A ",sep = ""),
                              collapse = " "),
                        collapse = " ")
    
    export_type <- "export_type=svg"
    
    file_ipath <- paste("./../Project_Output/",
                        column_name,"_ipath_all_enzymes_bins.svg",
                        sep="")
    
    ipath_command <- paste(paste("curl -d ",
                                 selections, " -d ", 
                                 export_type, 
                                 " https://pathways.embl.de/mapping.cgi -o ", 
                                 sep = '\"'),
                           file_ipath,
                           sep="")
    system(ipath_command)
    print(column_name)
    write_tsv(Filtered_proteins,file_name,col_names = FALSE)
    #pacman::p_load("rsvg")
    bitmap <- rsvg::rsvg_raw( here::here("Project_Output",
                                   paste0(column_name,"_ipath_all_enzymes_bins.svg")),
                              width = 3600)
    
    cowplot::ggdraw()+cowplot::draw_image(bitmap)+ cowplot::draw_figure_label(column_name,"top", size = 65)
    
    
}

calculateCoveredProtein_sav <- function (sample_id,proteinIDs, markerproteins) 
{
    compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", 
                      "N4", "C1", "C2", "C3", "C4", "C5", "M1", "M2")
    color.code <- c("gold", "orange", "salmon", "tomato2", "grey90", 
                    "grey70", "grey50", "grey30", "lightblue", "aquamarine", 
                    "cyan", "deepskyblue2", "turquoise3", "burlywood4", 
                    "tan4")
    compartment.size <- c(358, 351, 252, 174, 192, 121, 231, 
                          198, 242, 132, 220, 215, 341, 69, 269)
    covered.proteins <- intersect(proteinIDs, markerproteins)
    if (length(covered.proteins) < 1) 
        warning("There is no overlap between marker proteins and data!")
    c.marker.df <- SubCellBarCode::markerProteins[covered.proteins, 
    ]
    coverageCompWise <- lapply(seq_len(length(compartments)), 
                               function(x) {
                                   temp.df <- c.marker.df[c.marker.df$Compartments == 
                                                              compartments[x], ]
                                   values <- list(Compartments = compartments[x], ColorCode = color.code[x], 
                                                  ProteinCoverage = 100 * ((dim(temp.df)[1])/compartment.size[x]))
                               })
    coverage.df <- as.data.frame(do.call("rbind", coverageCompWise))
    non.enriched.loc <- coverage.df[coverage.df$ProteinCoverage < 
                                        20, ]
    if (nrow(non.enriched.loc) == 1) {
        warning("There is not enough enrichment at: ", as.character(non.enriched.loc$Compartments), 
                "\nWe recommend you to perform the fractionation, again.")
    }
    else if (nrow(non.enriched.loc) > 1) {
        comp <- paste(as.character(non.enriched.loc$Compartments), 
                      collapse = ",")
        warning("There are not enough enrichments at: ", comp, 
                "\nWe recommend you to perform the fractionation!")
    }
    coverage.df$ProteinCoverage <- as.numeric(coverage.df$ProteinCoverage)
    coverage.df$Compartments <- as.character(coverage.df$Compartments)
    coverage.df$ColorCode <- as.character(coverage.df$ColorCode)
    ggplot(data = coverage.df, aes(x = coverage.df$Compartments, 
                                         y = coverage.df$ProteinCoverage)) + geom_bar(stat = "identity", 
                                                                                      fill = coverage.df$ColorCode) + 
        scale_x_discrete(limits = c(compartments)) + 
              theme_bw() + theme(text = element_text(size = 5), plot.title = element_text(hjust = 0.5), 
                                 axis.text.x = element_text(face = "bold", color = "black"), 
                                 axis.text.y = element_text(face = "bold", color = "black")) +
        ylim(0, 100)+
        
              labs(title = paste("Marker Protein Coverage Sample",sample_id), 
                   y = "% Protein Coverage", x = "Compartment")

}
