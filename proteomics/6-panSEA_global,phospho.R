# Differential expression & enrichment analyses: global
# Chr8
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-01-02

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)

#### 1. Import metadata & crosstabs ####
setwd(
"~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")
global.df <- read.table(
  "global_data/Chr8_crosstab_global_gene_corrected.txt", 
  sep = "\t")
global.w.mouse.df <- read.table(
  "global_data/Chr8-with-mouse_crosstab_global_gene_corrected.txt", 
  sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
globalw.mouse.df$Gene <- rownames(global.w.mouse.df)

#### 2. Run panSEA across contrasts for each exp, omics type ####
# synapse IDs must match order of omics list
# synapse_id_map <- c("syn" = "global_data/",
#                     "syn" = "global_with_mouse_data/")

## prepare other input parameters
# identify contrasts
treatments <- na.omit(unique(meta.df$Treatment))

# set short-hand names for treatments for better plots
# IMPORTANT: these must have the same order as treatments above
treatment.names <- c("Ref", "WU-487CPP1", "MN-2", "JH-2-055", "JH-2-079",
                     "JH-2-002", "WU-22t Cshim81")

# define types based on short-hand contrasts
# by comparing everything to ref
types <- c()
for (i in 1:(length(treatment.names) - 1)) {
  types[i] <- paste(treatment.names[i + 1], "vs.", treatment.names[1])
}

# get contrasts with full-length names to extract data
contrasts <- c()
for (i in 1:length(types)) {
  # get treatment names (short-hand)
  contrast.treatment.names <- stringr::str_split(types[i], " vs. ")[[1]]
  
  # get full treatment names
  index1 <- treatments[treatment.names == contrast.treatment.names[1]]
  index2 <- treatments[treatment.names == contrast.treatment.names[2]]
  
  # define contrast with full treatment name
  contrasts[i] <- paste(index1, "vs.", index2)
}

# identify column names containing data for each treatment
# add X in front of tube names to match column names 
# after import (e.g., X5 for tube 5)
sample.groups <- list()
for (i in 1:length(treatments)) {
  sample.groups[[treatments[i]]] <- 
    paste0("X", na.omit(meta.df[meta.df$Treatment == treatments[i], ]$Tube))
}

#all.files2 <- c(DEG.files, GSEA.files, DMEA.files)
#subsets <- c("Differential expression", "GSEA", "DMEA")
#setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
#omics <- c("global", "phospho")
omics <- c("global", "global-with-mouse")
for (k in 1:length(omics)) {
  setwd(paste0(base.path, synapse_id_map[k]))
  
  ## prepare set annotations
  # generate gmt.features beforehand to save time
  if (grep("global", omics[k], ignore.case = TRUE)) {
    msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C2", "CP:KEGG")
    
    # extract necessary info into data frame
    msigdb.info <- as.data.frame(msigdb.info[, c(
      "gene_symbol",
      "gs_name",
      "gs_description"
    )])
    
    gmt <- DMEA::as_gmt(
      msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    ) 
  } else if (omics[k] == "phospho") {
    # only create gmt the first time
    # and use PNNL phospho feature IDs instead of ksdb
    SUB_SITE <- phospho.df$SUB_SITE
    phospho.ref <- data.frame(SUB_SITE)
    phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                  remove = FALSE)
    SUB_SITE <- NULL
    gmt <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE", min.per.set = 6)

    # save gmt for future analyses
    saveRDS(gmt, "gmt_PNNL_kinase-substrate_Chr8_5344.rds")
    
    #gmt <- readRDS("gmt_PNNL_kinase-substrate_Chr8_5344.rds")
  }
  
  # run panSEA for each omics type across all contrasts
  # CAUTION: this only works because the # of samples for each treatment type 
  # is equal; otherwise would have to run panSEA for each contrast separately 
  # and perhaps set the group.samples input parameter for panSEA
  
  # assemble inputs
  data.list <- list()
  gmt.features <- list()
  for (i in 1:length(contrasts)) {
    # identify treatments in this contrast
    contrast.treatments <- stringr::str_split(contrasts[i], " vs. ")[[1]]
    
    # assemble each input data set for each contrast
    if (grep("global", omics[k], ignore.case = TRUE)){
      if (grep("mouse", omics[k], ignore.case = TRUE)) {
        data.list[[types[i]]] <- 
          global.w.mouse.df[ , c("Gene", sample.groups[[contrast.treatments[1]]],
                         sample.groups[[contrast.treatments[2]]])]
      } else {
        data.list[[types[i]]] <- 
          global.df[ , c("Gene", sample.groups[[contrast.treatments[1]]],
                         sample.groups[[contrast.treatments[2]]])]
      }
      
      feature.names <- rep("Gene", length(types))
    } else if (omics[k] == "phospho") {
      data.list[[types[i]]] <- 
        phospho.df[ , c("SUB_SITE", sample.groups[[contrast.treatments[1]]],
                       sample.groups[[contrast.treatments[2]]])]
      
      feature.names <- rep("SUB_SITE", length(types))
    }
    
    gmt.features[[types[i]]] <- gmt
  }
  
  # run panSEA by querying CCLE data
  #Sys.setenv('R_MAX_VSIZE'=32000000000)
  #library(usethis) 
  #usethis::edit_r_environ()

  panSEA.results <- panSEA::panSEA(data.list, types, 
                                   feature.names = feature.names, 
                                   gmt.features = gmt.features)
  
  ## save results & upload to Synapse
  # store all results locally
  dir.create("analysis")
  setwd("analysis")
  #saveRDS(panSEA.results, file=paste0("Chr8_", omics[k], "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts
  panSEA.results <- readRDS(paste0("Chr8_", omics[k], "_panSEA_CCLE.rds"))
  
  # set file names
  DEG.files <- list("Differential_expression_results.csv" = 
                      panSEA.results$mDEG.results$compiled.results$results,
                    "Differential_expression_mean_results.csv" =
                      panSEA.results$mDEG.results$compiled.results$mean.results,
                    "Differential_expression_correlation_matrix.pdf" =
                      panSEA.results$mDEG.results$compiled.results$corr.matrix,
                    "Differential_expression_dot_plot.pdf" =
                      panSEA.results$mDEG.results$compiled.results$dot.plot)
  DMEA.files <- list("DMEA_results.csv" =
                       panSEA.results$mDMEA.results$compiled.results$results,
                     "DMEA_mean_results.csv" =
                       panSEA.results$mDMEA.results$compiled.results$mean.results,
                     "DMEA_correlation_matrix.pdf" =
                       panSEA.results$mDMEA.results$compiled.results$corr.matrix,
                     "DMEA_dot_plot.pdf" =
                       panSEA.results$mDMEA.results$compiled.results$dot.plot,
                     "DMEA_interactive_network.graph.html" =
                       panSEA.results$mDMEA.network$interactive)
  if (omics[k] == "phospho") {
    KSEA.files <- list("KSEA_results.csv" =
                         panSEA.results$mKSEA.results$compiled.results$results,
                       "KSEA_mean_results.csv" =
                         panSEA.results$mKSEA.results$compiled.results$mean.results,
                       "KSEA_correlation_matrix.pdf" =
                         panSEA.results$mKSEA.results$compiled.results$corr.matrix,
                       "KSEA_dot_plot.pdf" =
                         panSEA.results$mKSEA.results$compiled.results$dot.plot,
                       "KSEA_interactive_network.graph.html" =
                         panSEA.results$mKSEA.network$interactive)
    all.files <- list('Differential expression' = DEG.files,
                      'KSEA' = KSEA.files,
                      'DMEA' = DMEA.files)
  } else {
    GSEA.files <- list("GSEA_results.csv" =
                         panSEA.results$mGSEA.results$compiled.results$results,
                       "GSEA_mean_results.csv" =
                         panSEA.results$mGSEA.results$compiled.results$mean.results,
                       "GSEA_correlation_matrix.pdf" =
                         panSEA.results$mGSEA.results$compiled.results$corr.matrix,
                       "GSEA_dot_plot.pdf" =
                         panSEA.results$mGSEA.results$compiled.results$dot.plot,
                       "GSEA_interactive_network.graph.html" =
                         panSEA.results$mGSEA.network$interactive)
    all.files <- list('Differential expression' = DEG.files,
                      'GSEA' = GSEA.files,
                      'DMEA' = DMEA.files)
  }
  
  for (i in 1:length(all.files)) {
    # create local folder for subset of results
    setwd(paste0(base.path, synapse_id_map[k], "analysis"))
    dir.create(names(all.files)[i])
    setwd(names(all.files)[i])
    
    
    # save results locally
    temp.files <- all.files[[i]]
    
    CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
    for (j in 1:length(CSV.files)) {
      write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
    }
    
    PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
    for (j in 1:length(PDF.files)) {
      if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
        ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                        device = "pdf", width = 30)
      } else {
        ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                        device = "pdf")
      }
    }
    
    HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
    if (length(HTML.files) > 0) {
      for (j in 1:length(HTML.files)) {
        visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j])
      }
    }
    
    # create folder on Synpase for subset of results
    dataFolder <- 
        synapser::synStore(synapser::Folder(names(all.files)[i],
                                            parent = names(synapse_id_map)[k]))
    
    # upload results to Synapse
    CSVs <- lapply(as.list(CSV.files), synapser::File,
                    parent = dataFolder)
    lapply(CSVs, synapser::synStore)
    
    PDFs <- lapply(as.list(PDF.files), synapser::File,
                   parent = dataFolder)
    lapply(PDFs, synapser::synStore)
    
    if (length(HTML.files) > 0) {
      HTMLs <- lapply(HTML.files, synapser::File,
                                    parent = dataFolder)
      lapply(HTMLs, synapser::synStore) 
    }
  }
  panSEA.results <- NULL # make space to process next omics type
}
