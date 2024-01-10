# Differential expression & enrichment analyses: global
# Chr8
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-01-02

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)

#### 1. Import metadata & crosstabs ####
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
)

synapser::synGet("syn49554376", 
                 downloadLocation = getwd())
RNA.df <- read.table("salmon.merged.gene_counts.tsv", sep = "\t", 
                     header = TRUE)
RNA.df$gene_id <- NULL
colnames(RNA.df)[1] <- "Gene"

#### 2. Run panSEA across contrasts for each exp, omics type ####
# synapse IDs must match order of omics list
# synapse_id_map <- c("syn53217098" = "global_data/",
#                     "syn53217121" = "global_with_mouse_data/",
#                     "syn53217102" = "phospho_data/")

## prepare other input parameters
# identify samples - only 1 group here so no contrasts (no controls)
group.samples <- list(2:ncol(RNA.df))

# set types as treatment names (revised column names) for plots
# must match order of column names
types <- c("JH-2-009", "JH-2-055", "JH-2-079", 
           "WU-487", "WU-536", "WU-561", 
           "MN-2", "MN-3")

# identify column names containing data for each treatment
sample.groups <- list()
GSEA.rank.var <- c()
for (i in 1:length(types)) {
  sample.groups[[types[i]]] <- colnames(RNA.df)[i + 1]
  GSEA.rank.var <- c(GSEA.rank.var, colnames(RNA.df)[i + 1])
}

#all.files2 <- c(DEG.files, GSEA.files, DMEA.files)
#subsets <- c("Differential expression", "GSEA", "DMEA")
#setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
setwd(base.path)

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

## run panSEA across all samples
# average across duplicate gene names
RNA.df <- plyr::ddply(RNA.df, .(Gene), summarize,
                      X2.009_pdx = mean(X2.009_pdx, na.rm = TRUE),
                      X2.055_pdx = mean(X2.055_pdx, na.rm = TRUE),
                      X2.079_pdx = mean(X2.079_pdx, na.rm = TRUE),
                      WU.487_pdx = mean(WU.487_pdx, na.rm = TRUE),
                      WU.536_pdx = mean(WU.536_pdx, na.rm = TRUE),
                      WU.561_pdx = mean(WU.561_pdx, na.rm = TRUE),
                      mn2_pdx = mean(mn2_pdx, na.rm = TRUE),
                      mn3_pdx_rna = mean(mn3_pdx_rna, na.rm = TRUE)) # down to 60347 genes from 61598

# assemble inputs
data.list <- list()
gmt.features <- list()
for (i in 1:length(types)) {
  # assemble each input data set for each treatment
  data.list[[types[i]]] <- 
    RNA.df[ , c("Gene", sample.groups[[types[i]]])]
  colnames(data.list[[types[i]]])[2] <- "Chr8"
  
  feature.names <- rep("Gene", length(types))
  
  gmt.features[[types[i]]] <- gmt
}

#Sys.setenv('R_MAX_VSIZE'=32000000000)
#library(usethis) 
#usethis::edit_r_environ()

# run panSEA by querying CCLE data
panSEA.results <- panSEA::panSEA(data.list, types,
                                 GSEA.rank.var = rep("Chr8", length(types)),
                                 group.names = "Chr8", group.samples = 2,
                                 feature.names = feature.names,
                                 gmt.features = gmt.features)

# Running ssGSEA using JH-2-009 data
# Running enrichment analysis...
# Running ssGSEA using JH-2-055 data
# Running enrichment analysis...
# Error in if (temp.NES >= 0) { : missing value where TRUE/FALSE needed
#   In addition: Warning messages:
#     1: In GSEA_custom(input, gmt, num.permutations, stat.type, min.per.set,  :
#                         Removing drug sets with less than 6 drugs observed in data set
#                       2: In DMEA::drugSEA(data, gmt, feature.names, rank.var, "gs_name",  :
#                                             No enrichments met the FDR cut-off to produce mountain plots

# remove JH-2-055
data.list <- list(data.list[[1]], data.list[[3]], data.list[[4]],
                  data.list[[5]], data.list[[6]], data.list[[7]],
                  data.list[[8]])
types <- c(types[1], types[3:length(types)])
feature.names <- rep("Gene", 7)
gmt.features <- list(gmt.features[[1]], gmt.features[[3]], gmt.features[[4]],
                     gmt.features[[5]], gmt.features[[6]], gmt.features[[7]],
                     gmt.features[[8]])

panSEA.results <- panSEA::panSEA(data.list, types,
                                 GSEA.rank.var = rep("Chr8", length(types)),
                                 group.names = "Chr8", group.samples = 2,
                                 feature.names = feature.names,
                                 gmt.features = gmt.features)

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

## save results & upload to Synapse
# store all results locally
dir.create("analysis")
setwd("analysis")
#saveRDS(panSEA.results, file=paste0("Chr8_", "RNA-seq", "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts
# panSEA.results <- readRDS(paste0("Chr8_", omics[k], "_panSEA_CCLE.rds"))

for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path, "analysis"))
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
    ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], device = "pdf")
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
    }
  }
  
  # # create folder on Synpase for subset of results
  # dataFolder <- 
  #   synapser::synStore(synapser::Folder(names(all.files)[i],
  #                                       parent = names(synapse_id_map)[k]))
  # 
  # # upload results to Synapse
  # CSVs <- lapply(as.list(CSV.files), synapser::File,
  #                parent = dataFolder)
  # lapply(CSVs, synapser::synStore)
  # 
  # PDFs <- lapply(as.list(PDF.files), synapser::File,
  #                parent = dataFolder)
  # lapply(PDFs, synapser::synStore)
  # 
  # if (length(HTML.files) > 0) {
  #   HTMLs <- lapply(HTML.files, synapser::File,
  #                   parent = dataFolder)
  #   lapply(HTMLs, synapser::synStore) 
  # }
}
panSEA.results <- NULL # make space for next analysis

