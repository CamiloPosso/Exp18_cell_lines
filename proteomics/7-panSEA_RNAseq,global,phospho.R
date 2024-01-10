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
  "global_with_mouse_data/Chr8_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Chr8_crosstab_phospho_SiteID_corrected.txt", 
  sep = "\t")

synapser::synGet("syn49554376", 
                 downloadLocation = getwd())
RNA.df <- read.table("salmon.merged.gene_counts.tsv", sep = "\t", 
                     header = TRUE)

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
global.w.mouse.df$Gene <- rownames(global.w.mouse.df)
phospho.df$SUB_SITE <- rownames(phospho.df)

# add Tube number to metadata
meta.df$Tube <- NA
tube.num <- 0
for (i in 1:nrow(meta.df)) {
  if (meta.df$SampleID[i] != "reference") {
    tube.num <- tube.num + 1
    meta.df$Tube[i] <- tube.num
  }
}

#### 2. Run panSEA across omics types ####
# synapse IDs must match order of omics list
synapse_id_map <- c("syn" = "multi-omics",
                    "syn" = "Proteomics",
                    "syn" = "Proteomics and RNA-seq")

## prepare other input parameters
# identify samples - only 1 group here so no contrasts (no controls)
treatments <- unique(meta.df$SampleID[meta.df$SampleID != "reference"])

# set names for treatments for plots
meta.df$type <- NA
meta.df <- meta.df[order(meta.df$SampleID), ]

# types <- c("Global", "Global with Mouse", "Phospho", "RNA-seq")
types <- c("Global", "Phospho")

# identify column names containing data for each treatment
# add X in front of tube names to match column names 
# after import (e.g., X5 for tube 5)
sample.groups <- list()
GSEA.rank.var <- c()
for (i in 1:length(types)) {
  sample.groups[[types[i]]] <-
    paste0("X", na.omit(meta.df[meta.df$type == types[i], ]$Tube))
  GSEA.rank.var <- c(GSEA.rank.var, 
                     paste0("X", 
                            na.omit(meta.df[meta.df$type == types[i], ]$Tube)))
}

#all.files2 <- c(DEG.files, GSEA.files, DMEA.files)
#subsets <- c("Differential expression", "GSEA", "DMEA")
#setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
for (k in 1:length(synapse_id_map)) {
  setwd(paste0(base.path, synapse_id_map[k]))
  
  ## prepare set annotations
  # generate gmt.features beforehand to save time
  if (grepl("global", omics[k], ignore.case = TRUE)) {
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
  for (i in 1:length(types)) {
    # assemble each input data set for each treatment
    if (grepl("global", omics[k], ignore.case = TRUE)){
      if (grepl("mouse", omics[k], ignore.case = TRUE)) {
        data.list[[types[i]]] <- 
          global.w.mouse.df[ , c("Gene", sample.groups[[types[i]]])]
        colnames(data.list[[types[i]]])[2] <- "Chr8"
      } else {
        data.list[[types[i]]] <- 
          global.df[ , c("Gene", sample.groups[[types[i]]])]
      }
      
      feature.names <- rep("Gene", length(types))
    } else if (omics[k] == "phospho") {
      data.list[[types[i]]] <- 
        phospho.df[ , c("SUB_SITE", sample.groups[[types[i]]])]
      colnames(data.list[[types[i]]])[2] <- "Chr8"
      
      feature.names <- rep("SUB_SITE", length(types))
    }
    
    gmt.features[[types[i]]] <- gmt
  }

  #Sys.setenv('R_MAX_VSIZE'=32000000000)
  #library(usethis) 
  #usethis::edit_r_environ()

  if (grepl("global", omics[k], ignore.case = TRUE)) {
    if (grepl("mouse", omics[k], ignore.case = TRUE)) {
    # run panSEA by querying CCLE data
    panSEA.results <- panSEA::panSEA(data.list, types,
                                     GSEA.rank.var = rep("Chr8", length(types)),
                                     group.names = "Chr8", group.samples = 2,
                                     feature.names = feature.names,
                                     gmt.features = gmt.features)
    
    # panSEA.results <- panSEA::panSEA(list(global.df), types = "Proteomics",
    #                                  #GSEA.rank.var = GSEA.rank.var,
    #                                  group.names = "Chr8", 
    #                                  group.samples = list(GSEA.rank.var),
    #                                  #feature.names = feature.names,
    #                                  #gmt.features = gmt.features,
    #                                  gmt.features = list(gmt))
    
    DMEA.files <- list("DMEA_results.csv" =
                         panSEA.results$mDMEA.results[["Chr8"]]$compiled.results$results,
                       "DMEA_mean_results.csv" =
                         panSEA.results$mDMEA.results[["Chr8"]]$compiled.results$mean.results,
                       "DMEA_correlation_matrix.pdf" =
                         panSEA.results$mDMEA.results[["Chr8"]]$compiled.results$corr.matrix,
                       "DMEA_dot_plot.pdf" =
                         panSEA.results$mDMEA.results[["Chr8"]]$compiled.results$dot.plot,
                       "DMEA_interactive_network.graph.html" =
                         panSEA.results$mDMEA.network[["Chr8"]]$interactive)
    GSEA.files <- list("GSEA_results.csv" =
                         panSEA.results$mGSEA.results[["Chr8"]]$compiled.results$results,
                       "GSEA_mean_results.csv" =
                         panSEA.results$mGSEA.results[["Chr8"]]$compiled.results$mean.results,
                       "GSEA_correlation_matrix.pdf" =
                         panSEA.results$mGSEA.results[["Chr8"]]$compiled.results$corr.matrix,
                       "GSEA_dot_plot.pdf" =
                         panSEA.results$mGSEA.results[["Chr8"]]$compiled.results$dot.plot,
                       "GSEA_interactive_network.graph.html" =
                         panSEA.results$mGSEA.network[["Chr8"]]$interactive)
    all.files <- list('GSEA' = GSEA.files,
                      'DMEA' = DMEA.files)
    } else {
      # only 53 gene names for "global" (without mouse)! so skipping GSEA
      panSEA.results <- list()
      panSEA.results[["mDMEA.results"]] <- 
        panSEA::mDMEA(weights = data.list, types = types, 
                      weight.values = GSEA.rank.var)
      
      # compile inputs & outputs for network graph
      inputs <- list()
      for (i in 1:length(types)) {
        inputs[[types[i]]] <- 
          panSEA.results[["mDMEA.results"]]$all.results[[types[i]]]$corr.result
      }
      
      outputs <- list()
      for (i in 1:length(types)) {
        outputs[[types[i]]] <- 
          panSEA.results[["mDMEA.results"]]$all.results[[types[i]]]$result
      }
      
      panSEA.results[["mDMEA.network"]] <-
        panSEA::netSEA(
          inputs, outputs, rep("Drug", length(inputs)), 
          rank.var = rep("Pearson.est", length(inputs))
        ) 
      
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
      
      all.files <- list('DMEA' = DMEA.files)
    }
  } else {
    # run mGSEA (can't do DMEA because no baseline phospho CCLE data to query)
    panSEA.results <- list()
    panSEA.results[["mGSEA.results"]] <- panSEA::mGSEA(data.list, 
                                                       gmt = gmt.features, 
                                                       types, feature.names, 
                                    rank.var = rep("Chr8", length(types)))
    
    # compile inputs & results for network graph
    outputs <- list()
    for (i in 1:length(types)) {
      outputs[[types[i]]] <- panSEA.results[["mGSEA.results"]]$all.results$result
    }
    
    ssGSEA.network <- panSEA::netSEA(
      data.list, outputs, feature.names,
      GSEA.rank.var
    )
    
    KSEA.files <- list("KSEA_results.csv" =
                         panSEA.results$mGSEA.results$compiled.results$results,
                       "KSEA_mean_results.csv" =
                         panSEA.results$mGSEA.results$compiled.results$mean.results,
                       "KSEA_correlation_matrix.pdf" =
                         panSEA.results$mGSEA.results$compiled.results$corr.matrix,
                       "KSEA_dot_plot.pdf" =
                         panSEA.results$mGSEA.results$compiled.results$dot.plot,
                       "KSEA_interactive_network.graph.html" =
                         panSEA.results$mGSEA.network$interactive)
    all.files <- list('KSEA' = KSEA.files,
                      'DMEA' = DMEA.files)
  }
  
  ## save results & upload to Synapse
  # store all results locally
  dir.create("analysis")
  setwd("analysis")
  saveRDS(panSEA.results, file=paste0("Chr8_", omics[k], "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts
  # panSEA.results <- readRDS(paste0("Chr8_", omics[k], "_panSEA_CCLE.rds"))
  
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
                        device = "pdf", width = 35)
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
