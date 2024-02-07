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

DMEA1_prot <- panSEA::panSEA(data.list, types, GSEA = FALSE,
                          group.names=c("Chr8_amp", "not_amp"), 
                          group.samples = list(c("JH-2-079", "MN-2"),
                                               c("JH-2-055", "WU-487")),
                          expression = list("adherent CCLE", prot.df.noNA))

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
                   "DMEA_proteomics_mountain_plot_TGF_beta_receptor_inhibitor.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$mtn.plot[["TGF beta receptor inhibitor"]],
                   "DMEA_transcriptomics_mountain_plot_CHK_inhibitor.pdf" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$mtn.plot[["CHK inhibitor"]],
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

DMEA.prot.files <- list("DMEA_results.csv" =
                     DMEA1_prot$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     DMEA1_prot$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     DMEA1_prot$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     DMEA1_prot$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     DMEA1_prot$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     DMEA1_prot$mDMEA.network$interactive,
                   "DMEA_transcriptomics_volcano_plot.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Transcriptomics$volcano.plot,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_proteomics_mountain_plot_adenosine_receptor_antagonist.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["adenosine receptor antagonist"]],
                   "DMEA_proteomics_mountain_plot_FGFR_inhibitor.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["FGFR inhibitor"]],
                   "DMEA_proteomics_mountain_plot_RET_tyrosine_kinase_inhibitor.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["RET tyrosine kinase inhibitor"]],
                   "DMEA_proteomics_mountain_plot_ALK_tyrosine_kinase_receptor_inhibitor.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["ALK tyrosine kinase receptor inhibitor"]],
                   "DMEA_transcriptomics_scatter_plots.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Transcriptomics$corr.scatter.plots,
                   "DMEA_transcriptomics_correlations.csv" =
                     DMEA1_prot$mDMEA.results$all.results$Transcriptomics$corr.result,
                   "DMEA_proteomics_correlations.csv" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$corr.result)

all.files <- list('DMEA_prot' = DMEA.prot.files)

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/global_with_mouse_data/analysis/"
)
saveRDS(panSEA1_C1GSEA_custom, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional_custom.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1_C1GSEA, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1_hallmarkGSEA, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_GSEA_hallmark.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts
saveRDS(DMEA1_prot, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_DMEA_CCLE_proteomics.rds")) # 8.5 GB for 7 contrasts
panSEA1 <- readRDS(file=paste0("Chr8_", "RNAseq_and_global_PDX", "_panSEA_CCLE.rds"))

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
synapse_id_map <- c("syn53354456" = "global_with_mouse_data/")
k <- 1
for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path, synapse_id_map[k], "analysis/PDX_RNAseq_and_proteomics"))
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
### DMEA
data.list <- list(global.df, global.w.mouse.df)
types <- c("Proteomics", "PDX Proteomics")
Chr8.prot.samples <- paste0("X", meta.df[meta.df$Chr8_amp & 
                                           meta.df$Sample == "Sample", ]$Tube)
non.Chr8.prot.samples <- paste0("X", meta.df[!meta.df$Chr8_amp & 
                                           meta.df$Sample == "Sample", ]$Tube)
DMEA3 <- panSEA::panSEA(data.list, types, GSEA = FALSE, 
                        group.names = c("Chr8_amp","not_amp"),
                        group.samples = list(Chr8.prot.samples,
                                             non.Chr8.prot.samples))

data.list <- list(global.w.mouse.df, phospho.df)
types <- c("PDX Proteomics", "Phospho-proteomics")
gmt.ksea <- readRDS(paste0(
  base.path,"phospho_data/gmt_PNNL_kinase-substrate_Chr8_5344.rds"))
GSEA3 <- panSEA::panSEA(data.list, types, DMEA = FALSE, 
                        feature.names = c("Gene", "SUB_SITE"),
                        group.names = c("Chr8_amp","not_amp"),
                        group.samples = list(Chr8.prot.samples,
                                             non.Chr8.prot.samples),
                        gmt.features = list("msigdb_Homo sapiens_C2_CP:KEGG",
                                            gmt.ksea))
GSEA3_hallmark <- panSEA::panSEA(list(global.w.mouse.df), "PDX Proteomics", 
                                 DMEA = FALSE,
                                 group.names = c("Chr8_amp","not_amp"),
                                 group.samples = list(Chr8.prot.samples,
                                                      non.Chr8.prot.samples),
                                 gmt.features = list("msigdb_Homo sapiens_H"))
# GSEA3_positional <- panSEA::panSEA(list(global.w.mouse.df), "PDX Proteomics", 
#                                  DMEA = FALSE,
#                                  group.names = c("Chr8_amp","not_amp"),
#                                  group.samples = list(Chr8.prot.samples,
#                                                       non.Chr8.prot.samples),
#                                  gmt.features = list("msigdb_Homo sapiens_C1"))

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

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"

GSEA3_positional_custom <- panSEA::panSEA(list(global.w.mouse.df), "PDX Proteomics", 
                                   DMEA = FALSE,
                                   group.names = c("Chr8_amp","not_amp"),
                                   group.samples = list(Chr8.prot.samples,
                                                        non.Chr8.prot.samples),
                                   gmt.features = list(gmt))

GSEA.files <- list("GSEA_results.csv" =
                       GSEA3$mGSEA.results$compiled.results$results,
                     "GSEA_mean_results.csv" =
                       GSEA3$mGSEA.results$compiled.results$mean.results,
                     "GSEA_correlation_matrix.pdf" =
                       GSEA3$mGSEA.results$compiled.results$corr.matrix,
                     "GSEA_venn_diagram.pdf" =
                       GSEA3$mGSEA.results$compiled.results$venn.diagram,
                     "GSEA_dot_plot.pdf" =
                       GSEA3$mGSEA.results$compiled.results$dot.plot,
                     "GSEA_interactive_network.graph.html" =
                       GSEA3$mGSEA.network$interactive,
                     "GSEA_PDX_proteomics_volcano_plot.pdf" =
                       GSEA3$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot,
                     "GSEA_phospho-proteomics_volcano_plot.pdf" =
                       GSEA3$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot)

GSEA.H.files <- list("GSEA_PDX_proteomics_results.csv" =
                        GSEA3_hallmark$mGSEA.results$all.results$'PDX Proteomics'$result,
                      "GSEA_interactive_network.graph.html" =
                        GSEA3_hallmark$mGSEA.network$interactive,
                      "GSEA_PDX_proteomics_volcano_plot.pdf" =
                        GSEA3_hallmark$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot)

# GSEA.C1.files <- list("GSEA_PDX_proteomics_results.csv" =
#                         GSEA3_positional$mGSEA.results$all.results$'PDX Proteomics'$result,
#                       "GSEA_interactive_network.graph.html" =
#                         GSEA3_positional$mGSEA.network$interactive,
#                       "GSEA_PDX_proteomics_volcano_plot.pdf" =
#                         GSEA3_positional$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot)

# GSEA.C1.custom.files <- list("GSEA_PDX_proteomics_results.csv" =
#                                GSEA3_positional_custom$mGSEA.results$all.results$'PDX Proteomics'$result,
#                              "GSEA_interactive_network.graph.html" =
#                                GSEA3_positional_custom$mGSEA.network$interactive,
#                              "GSEA_PDX_proteomics_volcano_plot.pdf" =
#                                GSEA3_positional_custom$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot,
#                              "GSEA_PDX_proteomics_mountain_plot_Chr8_cancer-associated_genes.pdf" =
#                                GSEA3_positional_custom$mGSEA.results$all.results$'PDX Proteomics'$mtn.plots[["Chr8 cancer-associated genes"]])

# try using CCLE prot for prot DMEA
#### proteomics (CCLE from IMPROVE)
download.file("https://figshare.com/ndownloader/files/41466702", "proteomics.csv.gz")
prot.df <- read.csv(gzfile("proteomics.csv.gz"),fileEncoding="UTF-16LE")

allgenes = readr::read_csv("https://figshare.com/ndownloader/files/40576109")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()
#genes <- genes[genes$gene_symbol %in% colnames(RNA.df)[2:ncol(RNA.df)], ]

allsamples = readr::read_csv('https://figshare.com/ndownloader/files/40576103')
CCLE.samples <- dplyr::distinct(allsamples[allsamples$id_source == "CCLE",
                                           c("other_id","improve_sample_id")])

# merge prot.df with genes, samples to stop using improve IDs
prot.df <- merge(prot.df, CCLE.samples)
prot.df <- merge(prot.df, genes)
prot.df$entrez_id <- NULL
prot.df <- dplyr::distinct(prot.df)

# convert to wide format for DMEA
prot.df <- reshape2::dcast(prot.df, other_id ~ gene_symbol, mean,
                           value.var = "proteomics")
colnames(prot.df)[1] <- "CCLE_ID"
prot.df.noNA <- prot.df[, colSums(is.na(prot.df)) == 0] # 23304 gene names
# DMEA
data.list <- list(global.df, global.w.mouse.df)
types <- c("Proteomics", "PDX Proteomics")
Chr8.prot.samples <- paste0("X", meta.df[meta.df$Chr8_amp & 
                                           meta.df$Sample == "Sample", ]$Tube)
non.Chr8.prot.samples <- paste0("X", meta.df[!meta.df$Chr8_amp & 
                                               meta.df$Sample == "Sample", ]$Tube)
DMEA3_prot <- panSEA::panSEA(data.list, types, GSEA = FALSE, 
                        group.names = c("Chr8_amp","not_amp"),
                        group.samples = list(Chr8.prot.samples,
                                             non.Chr8.prot.samples),
                        expression = list(prot.df.noNA, prot.df.noNA))

DEG.files <- list("DEG_results.csv" =
                    DMEA3$mDEG.results$compiled.results$results,
                  "DEG_mean_results.csv" =
                    DMEA3$mDEG.results$compiled.results$mean.results,
                  "DEG_correlation_matrix.pdf" =
                    DMEA3$mDEG.results$compiled.results$corr.matrix,
                  "DEG_venn_diagram.pdf" =
                    DMEA3$mDEG.results$compiled.results$venn.diagram,
                  "DEG_dot_plot.pdf" =
                    DMEA3$mDEG.results$compiled.results$dot.plot,
                  "DEG_results_phospho-proteomics.csv" = 
                    GSEA3$mDEG.results$all.results$`Phospho-proteomics`,
                  "DEG_dot_plot_phospho-proteomics.pdf" = 
                    GSEA3$mDEG.results$compiled.results$dot.plot)
DMEA.files <- list("DMEA_results.csv" =
                     DMEA3$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     DMEA3$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     DMEA3$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     DMEA3$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     DMEA3$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     DMEA3$mDMEA.network$interactive,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     DMEA3$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_PDX_proteomics_volcano_plot.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$volcano.plot,
                   "DMEA_PDX_proteomics_mountain_plot_TGF_beta_receptor_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["TGF beta receptor inhibitor"]],
                   "DMEA_PDX_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_PDX_proteomics_mountain_plot_RAF_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["RAF inhibitor"]],
                   "DMEA_proteomics_mountain_plot_TGF_beta_receptor_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'Proteomics'$mtn.plot[["TGF beta receptor inhibitor"]],
                   "DMEA_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_proteomics_mountain_plot_RAF_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'Proteomics'$mtn.plot[["RAF inhibitor"]],
                   "DMEA_proteomics_correlations.csv" =
                     DMEA3$mDMEA.results$all.results$Proteomics$corr.result,
                   "DMEA_PDX_proteomics_correlations.csv" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$corr.result)
DMEA.prot.files <- list("DMEA_results.csv" =
                     DMEA3_prot$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     DMEA3_prot$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     DMEA3_prot$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     DMEA3_prot$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     DMEA3_prot$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     DMEA3_prot$mDMEA.network$interactive,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_PDX_proteomics_volcano_plot.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$volcano.plot,
                   "DMEA_PDX_proteomics_mountain_plot_FGFR_inhibitor.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["FGFR inhibitor"]],
                   "DMEA_PDX_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_proteomics_correlations.csv" =
                     DMEA3_prot$mDMEA.results$all.results$Proteomics$corr.result,
                   "DMEA_PDX_proteomics_correlations.csv" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$corr.result)

all.files <- list('DMEA' = DMEA.files,
                  'DMEA_proteomics' = DMEA.prot.files,
                  'GSEA' = GSEA.files,
                  'GSEA_hallmark' = GSEA.H.files,
                  'Differential expression' = DEG.files)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
synapse_id_map <- c("syn53357947" = "global_with_mouse_data/")
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

# RNAseq only
RNA.df3 <- plyr::ddply(RNA.df, .(Gene), summarize,
                       'WU-487' = mean(WU.487_pdx, na.rm = TRUE),
                       'MN-2' = mean(mn2_pdx, na.rm = TRUE),
                       'JH-2-055' = mean(X2.055_pdx, na.rm = TRUE),
                       'JH-2-079' = mean(X2.079_pdx, na.rm = TRUE),
                       'WU-561' = mean(WU.561_pdx, na.rm = TRUE))
RNA.chr8 <- c("JH-2-079", "WU-561", "MN-2")
RNA.non.chr8 <- c("JH-2-055", "WU-487")
library(Biobase)
RNApanSEA <- panSEA::panSEA(list(RNA.df3), "Transcriptomics", 
                            group.samples = list(RNA.chr8, RNA.non.chr8))

DEG.files <- list("DEG_results.csv" =
                    RNApanSEA$mDEG.results$all.results$Transcriptomics)
DMEA.files <- list("DMEA_results.csv" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$result,
                   "DMEA_interactive_network.graph.html" =
                     RNApanSEA$mDMEA.network$interactive,
                   "DMEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$volcano.plot,
                   "DMEA_transcriptomics_mountain_plot_HMGCR_inhibitor.pdf" =
                     RNApanSEA$mDMEA.results$all.results$'Transcriptomics'$mtn.plot[["HMGCR inhibitor"]],
                   "DMEA_transcriptomics_mountain_plot_protein_synthesis_inhibitor.pdf" =
                     RNApanSEA$mDMEA.results$all.results$'Transcriptomics'$mtn.plot[["protein synthesis inhibitor"]],
                   "DMEA_transcriptomics_scatter_plots.pdf" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$corr.scatter.plots,
                   "DMEA_transcriptomics_correlations.csv" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$corr.result)
GSEA.files <- list("GSEA_results.csv" =
                     RNApanSEA$mGSEA.results$all.results$Transcriptomics$result,
                   "GSEA_interactive_network.graph.html" =
                     RNApanSEA$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_transcriptomics_correlations.csv" =
                     RNApanSEA$mGSEA.results$all.results$Transcriptomics$corr.result)
RNApanSEA_H <- panSEA::panSEA(list(RNA.df3), "Transcriptomics", DMEA = FALSE,
                            group.samples = list(RNA.chr8, RNA.non.chr8),
                            gmt.features = list("msigdb_Homo sapiens_H"))

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

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"

RNApanSEA_pos_custom <- panSEA::panSEA(list(RNA.df3), "Transcriptomics", DMEA = FALSE,
                              group.samples = list(RNA.chr8, RNA.non.chr8),
                              gmt.features = list(gmt))

GSEA.H.files <- list("GSEA_results.csv" =
                     RNApanSEA_H$mGSEA.results$all.results$Transcriptomics$result,
                   "GSEA_interactive_network.graph.html" =
                     RNApanSEA_H$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA_H$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_transcriptomics_correlations.csv" =
                     RNApanSEA_H$mGSEA.results$all.results$Transcriptomics$corr.result,
                   "GSEA_transcriptomics_mountain_plot_HALLMARK_MYC_TARGETS_V2.pdf" =
                     RNApanSEA_H$mGSEA.results$all.results$'Transcriptomics'$mtn.plot[["HALLMARK_MYC_TARGETS_V2"]],
                   "GSEA_transcriptomics_mountain_plot_HALLMARK_MYC_TARGETS.pdf" =
                     RNApanSEA_H$mGSEA.results$all.results$'Transcriptomics'$mtn.plot[["HALLMARK_MYC_TARGETS"]])
GSEA.pos.custom.files <- list("GSEA_results.csv" =
                     RNApanSEA_pos_custom$mGSEA.results$all.results$Transcriptomics$result,
                   "GSEA_interactive_network.graph.html" =
                     RNApanSEA_pos_custom$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA_pos_custom$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_transcriptomics_correlations.csv" =
                     RNApanSEA_pos_custom$mGSEA.results$all.results$Transcriptomics$corr.result)

all.files <- list('DMEA' = DMEA.files,
                  'GSEA' = GSEA.files,
                  'GSEA_HALLMARK' = GSEA.H.files,
                  'GSEA_POSITIONAL_CUSTOM' = GSEA.pos.custom.files,
                  'Differential expression' = DEG.files)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
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
}
setwd(paste0(base.path, "analysis"))
saveRDS(RNApanSEA_pos_custom, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional_custom.rds")) # 8.5 GB for 7 contrasts
saveRDS(RNApanSEA_H, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_hallmark.rds")) # 8.5 GB for 7 contrasts
saveRDS(RNApanSEA, file=paste0("Chr8_", "RNAseq_PDX", "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts
