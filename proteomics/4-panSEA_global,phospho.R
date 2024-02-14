# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 24 TMT
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-02-09

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)
library(dplyr)

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
source("panSEA_helper.R")

# to do:
# revise run_contrasts to make sure everything is saved correctly (incl mtn plots)

# overview
#### 1. Import metadata & crosstabs ####
#### 2. Import BeatAML data formatted for DMEA ####
#### 3. Run panSEA across omics for each contrast ####

#### 1. Import metadata & crosstabs ####
setwd("data")
meta.df <- readxl::read_excel("Exp24metadataTable_TMT.xlsx") 
meta.df$MeasurementName <- as.character(rownames(meta.df))

# add other drug info & make sure sensitivity is correctly labeled
sens.info <- read.csv("Exp24_drug_sensitivity_20240209.csv")
meta.df <- merge(meta.df, sens.info)
meta.df$X <- NULL

global.df <- read.table(
  "global_data/Exp24_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Exp24_crosstab_phospho_SiteID_corrected.txt", 
  sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
phospho.df$SUB_SITE <- rownames(phospho.df)

#### 2. Import BeatAML data formatted for DMEA ####
# import drug MOA annotations
moa.BeatAML <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

# login to Synapse
#synLogin()

# load data from Synapse
BeatAML.data <- load_BeatAML_for_DMEA("BeatAML_DMEA_inputs")

#### 3. Run panSEA across sort.types for each omics, sens.type, drug.type ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/"
# synapse IDs must match order of omics list
synapse_id_map <- c("syn53180717" = "global_data/",
                    "syn53180718" = "phospho_data/")

# set up comparisons
#sort.types <- c("CD14+", "CD34+", "MSC")
sort.types <- unique(na.omit(meta.df$SampleType))
#sens.types <- c("sensitive", "resistant")
sens.types <- unique(na.omit(meta.df$Drug))
drug.types <- c("Aza", "Ven", "Aza.Ven")

# sort.type contrasts
sort.contrasts <- list()
counter <- 0
for (i in 1:length(sort.types)) {
  for (j in 1:length(sort.types)) {
    if (i != j) {
      counter <- counter + 1
      alpha.pair <- sort(c(sort.types[i], sort.types[j]))
      sort.contrasts[[counter]] <- alpha.pair
    }
  }
}
sort.contrasts <- unique(sort.contrasts)

# drug.sens contrasts
sens.contrasts <- list()
counter <- 0
for (i in drug.types) {
  sens.contrasts[[i]] <- c("Sensitive", "Resistant")
}

# get set annotations for pathway analyses
## prepare set annotations
# generate gmt.features beforehand to save time
gmt <- readRDS("gmt_MSigDB_Homo-sapiens_C2_CP_KEGG.rds")

msigdb.info <- msigdbr::msigdbr("Homo sapiens", "H")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt.H <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 
saveRDS(gmt.H, "gmt_MSigDB_Homo-sapiens_H.rds")

gmt.ksdb <- readRDS("gmt_ksdb_human_PNNL.rds")

# only create substrate gmt the first time
SUB_SITE <- phospho.df$SUB_SITE
phospho.ref <- data.frame(SUB_SITE)
phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                              remove = FALSE)
SUB_SITE <- NULL
gmt.sub <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE", min.per.set = 6)
saveRDS(gmt.sub, "gmt_PNNL_kinase-substrate_PTRC2_exp24.rds")
#gmt <- readRDS("gmt_PNNL_kinase-substrate_PTRC2_exp24.rds")

gmt.drug <- DMEA::as_gmt(moa.BeatAML, min.per.set = 6)
saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA.rds")

omics <- list("global" = global.df,
              "phospho" = phospho.df)
beatAML <- list("global" = BeatAML.data$global,
                "phospho" = BeatAML.data$phospho)

# prepare universal inputs
gmt.features <- list(gmt, gmt.ksdb)
gmt.features2 <- list(gmt.H, gmt.sub)
gmt.features3a <- list(gmt, gmt) # for phospho GSEA^2
gmt.features3b <- list(gmt.H, gmt.H) # for phospho GSEA^2

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/"
)
dir.create("analysis")
setwd("analysis")
for (j in drug.types) {
  for (i in sens.types) {
    # filter meta.df
    meta.filtered <- meta.df[meta.df[ , j] == i, ]
    
    # run cell type contrasts for each drug sensitivity
    run_contrasts_global_phospho(sort.contrasts, "SampleType", meta.filtered,
                                 omics, beatAML, gmt.features, gmt.features2, 
                                 gmt.features3a, gmt.features3b, 
                                 gmt.drug, BeatAML.data$drug, 
                                 synapse_id = "syn53521265")
  } 
}
# run cell type contrasts without drug sensitivity filter
run_contrasts_global_phospho(sort.contrasts, "SampleType", meta.df,
                             omics, beatAML, gmt.features, gmt.features2, 
                             gmt.features3a, gmt.features3b, 
                             gmt.drug, BeatAML.data$drug, 
                             synapse_id = "syn53521265")

for (i in sort.types) {
  # filter meta.df
  meta.filtered <- meta.df[meta.df$SampleType == i, ]
  
  # run sens vs. res contrasts for each cell type
  for (j in drug.types) {
    run_contrasts_global_phospho(list(c("Sensitive", "Resistant")), j, meta.filtered, 
                                 omics, beatAML, gmt.features, gmt.features2, 
                                 gmt.features3, gmt.drug, BeatAML.data$drug, 
                                 synapse_id = "syn53521265")
  }
}
# run sens vs. res contrasts without cell type filter
for (j in drug.types) {
  run_contrasts_global_phospho(list(c("Sensitive", "Resistant")), j, meta.df,
                               omics, beatAML, gmt.features, gmt.features2, 
                               gmt.features3, gmt.drug, BeatAML.data$drug, 
                               synapse_id = "syn53521265")
}