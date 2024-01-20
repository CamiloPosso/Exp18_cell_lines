# Differential expression & enrichment analyses: RNAseq, global PDX; global, phospho cell lines
# Chr8: amplified vs. not amplified
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-01-19

#### first analysis: only shared MPNST PDX samples w known Chr8 status
#### second analysis: only shared MPNST samples w known Chr8 status
#### third analysis: all MPNST samples w known Chr8 status

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(plyr)

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

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
)
synapser::synGet("syn49554376", 
                 downloadLocation = getwd())
RNA.df <- read.table("salmon.merged.gene_counts.tsv", sep = "\t", 
                     header = TRUE)
RNA.df$gene_id <- NULL
colnames(RNA.df)[1] <- "Gene"

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


# add Chr8 status to metadata
meta.df$Chr8_amp <- TRUE
meta.df[meta.df$Sample == "Reference" | 
          meta.df$SampleID %in% c("JH-2-055", "WU-487CPP1"), ]$Chr8_amp <- FALSE

#### first analysis: only shared MPNST PDX samples w known Chr8 status
# identify global samples in RNAseq
meta.df$RNAseq <- FALSE
prot.samples.in.RNAseq <- c("WU-487CPP1", "MN-2", "JH-2-055", "JH-2-079")
meta.df[meta.df$SampleID %in% prot.samples.in.RNAseq, ]$RNAseq <- TRUE
prot.tubes.in.RNAseq <- paste0("X", meta.df[meta.df$RNAseq, ]$Tube)
#global.w.mouse.df1 <- global.w.mouse.df[ , c("Gene", prot.tubes.in.RNAseq)]
global.w.mouse.df1 <- plyr::ddply(global.w.mouse.df, .(Gene), summarize,
                                 'WU-487' = mean(c(X1, X9), na.rm = TRUE),
                                 'MN-2' = mean(c(X2, X8), na.rm = TRUE),
                                 'JH-2-055' = mean(c(X3, X6), na.rm = TRUE),
                                 'JH-2-079' = mean(c(X4, X12), na.rm = TRUE))

# identify RNAseq samples in global
RNAseq.samples.in.prot <- c("WU.487_pdx", "mn2_pdx", "X2.055_pdx", "X2.079_pdx")
RNA.df1 <- RNA.df[ , c("Gene", RNAseq.samples.in.prot)]
RNA.df1 <- plyr::ddply(RNA.df1, .(Gene), summarize,
                       'WU-487' = mean(WU.487_pdx, na.rm = TRUE),
                       'MN-2' = mean(mn2_pdx, na.rm = TRUE),
                       'JH-2-055' = mean(X2.055_pdx, na.rm = TRUE),
                       'JH-2-079' = mean(X2.079_pdx, na.rm = TRUE))
# colnames(RNA.df1) <- colnames(global.w.mouse.df1)
# RNA.df1 <- plyr::ddply(RNA.df1, .(Gene), summarize,
#                                   'WU-487' = mean('WU-487', na.rm = TRUE),
#                                   'MN-2' = mean('MN-2', na.rm = TRUE),
#                                   'JH-2-055' = mean('JH-2-055', na.rm = TRUE),
#                                   'JH-2-079' = mean('JH-2-079', na.rm = TRUE))

# run panSEA
data.list <- list(RNA.df1, global.w.mouse.df1)
types <- c("Transcriptomics", "Proteomics")
library(Biobase)
panSEA1 <- panSEA::panSEA(data.list, types, 
                          group.names=c("Chr8_amp", "not_amp"), 
                          group.samples = list(c("JH-2-079", "MN-2"),
                                               c("JH-2-055", "WU-487")))

DEG.files <- list("DEG_results.csv" =
                     panSEA1$mDEG.results$compiled.results$results,
                   "DEG_mean_results.csv" =
                     panSEA1$mDEG.results$compiled.results$mean.results,
                   "DEG_correlation_matrix.pdf" =
                     panSEA1$mDEG.results$compiled.results$corr.matrix,
                   "DEG_venn_diagram.pdf" =
                     panSEA1$mDEG.results$compiled.results$venn.diagram,
                   "DEG_dot_plot.pdf" =
                     panSEA1$mDEG.results$compiled.results$dot.plot)
DMEA.files <- list("DMEA_results.csv" =
                     panSEA1$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     panSEA1$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     panSEA1$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     panSEA1$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     panSEA1$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     panSEA1$mDMEA.network$interactive,
                   "DMEA_transcriptomics_volcano_plot.pdf" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$volcano.plot,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_proteomics_mountain_plot_adenosine_receptor_antagonist.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$mtn.plot[["adenosine receptor antagonist"]],
                   "DMEA_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_transcriptomics_mountain_plot_CHK_inhibitor.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$mtn.plot[["CHK inhibitor"]],
                   "DMEA_proteomics_scatter_plots.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$corr.scatter.plots,
                   "DMEA_transcriptomics_scatter_plots.pdf" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$corr.scatter.plots,
                   "DMEA_transcriptomics_correlations.csv" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$corr.result,
                   "DMEA_proteomics_correlations.csv" =
                     panSEA1$mDMEA.results$all.results$Proteomics$corr.result)
GSEA.files <- list("GSEA_results.csv" =
                     panSEA1$mGSEA.results$compiled.results$results,
                   "GSEA_mean_results.csv" =
                     panSEA1$mGSEA.results$compiled.results$mean.results,
                   "GSEA_correlation_matrix.pdf" =
                     panSEA1$mGSEA.results$compiled.results$corr.matrix,
                   "GSEA_venn_diagram.pdf" =
                     panSEA1$mGSEA.results$compiled.results$venn.diagram,
                   "GSEA_dot_plot.pdf" =
                     panSEA1$mGSEA.results$compiled.results$dot.plot,
                   "GSEA_interactive_network.graph.html" =
                     panSEA1$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     panSEA1$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_proteomics_volcano_plot.pdf" =
                     panSEA1$mGSEA.results$all.results$Proteomics$volcano.plot,
                   "GSEA_transcriptomics_mountain_plot_KEGG_RIBOSOME.pdf" =
                     panSEA1$mGSEA.results$all.results$Transcriptomics$mtn.plot[["KEGG_RIBOSOME"]])

# generate gmt.features beforehand to save time
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "H")

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
gmt.features <- list(gmt, gmt)
panSEA1_hallmarkGSEA <- panSEA::panSEA(data.list, types, DMEA = FALSE,
                              group.names=c("Chr8_amp", "not_amp"), 
                              group.samples = list(c("JH-2-079", "MN-2"), 
                                                   c("JH-2-055", "WU-487")),
                              gmt.features = gmt.features,
                              scatter.plots = FALSE)
panSEA1_hallmarkGSEA <- panSEA1_C1GSEA
# generate gmt.features beforehand to save time
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

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
gmt.features <- list(gmt, gmt)
panSEA1_C1GSEA <- panSEA::panSEA(data.list, types, DMEA = FALSE,
                                       group.names=c("Chr8_amp", "not_amp"), 
                                       group.samples = list(c("JH-2-079", "MN-2"), 
                                                            c("JH-2-055", "WU-487")),
                                       gmt.features = gmt.features, 
                                 scatter.plots = FALSE)
# proteomics didn't have coverage for 2+ gene sets
panSEA1_C1GSEA <- panSEA::panSEA(list(RNA.df1), "Transcriptomics", DMEA = FALSE,
                                 group.names=c("Chr8_amp", "not_amp"), 
                                 group.samples = list(c("JH-2-079", "MN-2"), 
                                                      c("JH-2-055", "WU-487")),
                                 gmt.features = list(gmt), 
                                 scatter.plots = FALSE)

# try w Chr8 cancer-associated genes
Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"
gmt.features <- list(gmt, gmt)
panSEA1_C1GSEA_custom <- panSEA::panSEA(data.list, types, DMEA = FALSE,
                                 group.names=c("Chr8_amp", "not_amp"), 
                                 group.samples = list(c("JH-2-079", "MN-2"), 
                                                      c("JH-2-055", "WU-487")),
                                 gmt.features = gmt.features, 
                                 scatter.plots = FALSE)
panSEA1_C1GSEA_custom <- panSEA::panSEA(list(RNA.df1), "Transcriptomics", DMEA = FALSE,
                                        group.names=c("Chr8_amp", "not_amp"), 
                                        group.samples = list(c("JH-2-079", "MN-2"), 
                                                             c("JH-2-055", "WU-487")),
                                        gmt.features = list(gmt), 
                                        scatter.plots = FALSE)

GSEA.H.files <- list("GSEA_results.csv" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$results,
                   "GSEA_mean_results.csv" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$mean.results,
                   "GSEA_correlation_matrix.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$corr.matrix,
                   "GSEA_venn_diagram.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$venn.diagram,
                   "GSEA_dot_plot.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$dot.plot,
                   "GSEA_interactive_network.graph.html" =
                     panSEA1_hallmarkGSEA$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_proteomics_volcano_plot.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$all.results$Proteomics$volcano.plot,
                   "GSEA_transcriptomics_mountain_plot_HALLMARK_MYC_TARGETS_V2.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$all.results$Transcriptomics$mtn.plots[["HALLMARK_MYC_TARGETS_V2"]])

GSEA.C1.files <- list("GSEA_transcriptomics_results.csv" =
                       panSEA1_C1GSEA$mGSEA.results$all.results$Transcriptomics$result,
                      "GSEA_interactive_network.graph.html" =
                        panSEA1_C1GSEA$mGSEA.network$interactive,
                     "GSEA_transcriptomics_volcano_plot.pdf" =
                       panSEA1_C1GSEA$mGSEA.results$all.results$Transcriptomics$volcano.plot)

GSEA.C1.custom.files <- list("GSEA_transcriptomics_results.csv" =
                        panSEA1_C1GSEA_custom$mGSEA.results$all.results$Transcriptomics$result,
                      "GSEA_interactive_network.graph.html" =
                        panSEA1_C1GSEA_custom$mGSEA.network$interactive,
                      "GSEA_transcriptomics_volcano_plot.pdf" =
                        panSEA1_C1GSEA_custom$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                      "GSEA_transcriptomics_mountain_plot_Chr8_cancer-associated_genes.pdf" =
                        panSEA1_C1GSEA_custom$mGSEA.results$all.results$Transcriptomics$mtn.plots[["Chr8 cancer-associated genes"]])

all.files <- list('GSEA_KEGG_LEGACY' = GSEA.files,
                  'DMEA' = DMEA.files,
                  'GSEA_HALLMARK' = GSEA.H.files,
                  'GSEA_POSITIONAL' = GSEA.C1.files
                  #,
                  #'GSEA_POSITIONAL_CUSTOM' = GSEA.C1.custom.files
                  )
all.files <- list('GSEA_POSITIONAL' = GSEA.C1.files,
                  'GSEA_POSITIONAL_CUSTOM' = GSEA.C1.custom.files
)
all.files <- list('Differential expression' = DEG.files)
all.files <- list('DMEA' = DMEA.files)

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/global_with_mouse_data/analysis/"
)
saveRDS(panSEA1_C1GSEA_custom, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional_custom.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1_C1GSEA, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1_hallmarkGSEA, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_GSEA_hallmark.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
synapse_id_map <- c("syn53354456" = "global_with_mouse_data/")
k <- 1
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
    # if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf", width = 35)
    # } else {
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf")
    # }
    ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                    device = "pdf")
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
#panSEA.results <- NULL # make space to process next omics type


#### second analysis: only shared MPNST samples w known Chr8 status
#### third analysis: all MPNST samples w known Chr8 status



#### 2. Run panSEA across omics types ####
# synapse IDs must match order of omics list
# synapse_id_map <- c("syn" = "multi-omics",
#                     "syn" = "Proteomics",
#                     "syn" = "Proteomics and RNA-seq")

## prepare other input parameters
# # identify samples - only 1 group here so no contrasts (no controls)
# treatments <- unique(meta.df$SampleID[meta.df$SampleID != "reference"])
# 
# # set names for treatments for plots
# meta.df$type <- NA
# meta.df <- meta.df[order(meta.df$SampleID), ]

# types <- c("Global", "Global with Mouse", "Phospho", "RNA-seq")
types <- c("Transcriptomics", "Proteomics")

# identify column names containing data for each treatment
# add X in front of tube names to match column names 
# after import (e.g., X5 for tube 5)
sample.groups <- colnames(global.w.mouse.df)

setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/")
dir.create("multi-omics")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/multi-omics/"

setwd(paste0(base.path, synapse_id_map[k]))

## prepare set annotations
# generate gmt.features beforehand to save time
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
      data.list[[types[i]]] <- 
        global.w.mouse.df[ , c("Gene", sample.groups[[types[i]]])]
      colnames(data.list[[types[i]]])[2] <- "Chr8"
    
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

