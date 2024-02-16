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
## set up comparisons
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

# re-do sort.type contrasts without '+' at end of names
sort.types2 <- sort.types
for (i in 1:length(sort.types2)) {
  if (substr(sort.types2[i], nchar(sort.types2[i]), nchar(sort.types2[i])) == "+") {
    sort.types2[i] <- paste(substr(sort.types2[i], 1, nchar(sort.types2[i])-1), "Pos")
  }
}

meta.df$SampleType2 <- meta.df$SampleType
for (i in 1:nrow(meta.df)) {
  temp.type <- meta.df$SampleType[i]
  if (substr(temp.type, nchar(temp.type), nchar(temp.type)) == "+") {
    meta.df$SampleType2[i] <- paste(substr(temp.type, 1, nchar(temp.type)-1), "Pos")
  }
}

# re-do sort.type contrasts without spaces in names
sort.types3 <- sort.types
for (i in 1:length(sort.types3)) {
  if (grepl(" ", sort.types3[i])) {
    sort.types3[i] <- sub(" ", "_", sort.types3[i])
  }
}

meta.df$SampleType3 <- meta.df$SampleType
for (i in 1:nrow(meta.df)) {
  temp.type <- meta.df$SampleType[i]
  if (grepl(" ", temp.type)) {
    meta.df$SampleType3[i] <- sub(" ", "_", temp.type)
  }
}

sort.contrasts3 <- list()
counter <- 0
for (i in 1:length(sort.types3)) {
  for (j in 1:length(sort.types3)) {
    if (i != j) {
      counter <- counter + 1
      alpha.pair <- sort(c(sort.types3[i], sort.types3[j]))
      sort.contrasts3[[counter]] <- alpha.pair
    }
  }
}
sort.contrasts3 <- unique(sort.contrasts3)

# re-do sort.type contrasts without spaces or + in names
sort.types4 <- sort.types2
for (i in 1:length(sort.types4)) {
  if (grepl(" ", sort.types4[i])) {
    sort.types4[i] <- sub(" ", "_", sort.types4[i])
  }
}

meta.df$SampleType4 <- meta.df$SampleType
meta.df[meta.df$SampleType == "CD14+",]$SampleType4 <- "CD14_Pos"
meta.df[meta.df$SampleType == "CD34+",]$SampleType4 <- "CD34_Pos"
meta.df[meta.df$SampleType == "CD14+ Flow",]$SampleType4 <- "CD14_Pos_Flow"
meta.df[meta.df$SampleType == "CD34+ Flow",]$SampleType4 <- "CD34_Pos_Flow"
meta.df[meta.df$SampleType == "MSC Flow",]$SampleType4 <- "MSC_Flow"
sort.types4 <- na.omit(unique(meta.df$SampleType4))
sort.contrasts4 <- list()
counter <- 0
for (i in 1:length(sort.types4)) {
  for (j in 1:length(sort.types4)) {
    if (i != j) {
      counter <- counter + 1
      alpha.pair <- sort(c(sort.types4[i], sort.types4[j]))
      sort.contrasts4[[counter]] <- alpha.pair
    }
  }
}
sort.contrasts4 <- unique(sort.contrasts4)

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
# SUB_SITE <- phospho.df$SUB_SITE
# phospho.ref <- data.frame(SUB_SITE)
# phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
#                                               remove = FALSE)
# SUB_SITE <- NULL
# gmt.sub <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE", min.per.set = 6)
# saveRDS(gmt.sub, "gmt_PNNL_kinase-substrate_PTRC2_exp24.rds")
gmt.sub <- readRDS("gmt_PNNL_kinase-substrate_PTRC2_exp24.rds")

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

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
# add 'X' to MeasurementName to match colnames of data
meta.df$id <- paste0('X', meta.df$MeasurementName)

synapse_id <- "syn53606822"

# run sens vs. res contrasts without cell type filter
for (j in drug.types) {
  setwd(base.path)
  run_contrasts_global_phospho(list(c("Sensitive", "Resistant")), j, 'id',
                               meta.df, omics, beatAML, gmt.features, 
                               gmt.features2, gmt.features3a, gmt.features3b, 
                               gmt.drug, BeatAML.data$drug, base.path,
                               synapse_id)
}
for (i in sort.types) {
  # filter meta.df
  meta.filtered <- meta.df[meta.df$SampleType == i, ]
  
  # create folder for each cell type
  setwd(base.path)
  dir.create(file.path(i))
  setwd(file.path(i))
  temp.base.path <- file.path(base.path, i)
  dataFolder <- 
    synapser::synStore(synapser::Folder(file.path(i),
                                        parent = synapse_id))
  
  # run sens vs. res contrasts for each cell type
  for (j in drug.types) {
    if (nrow(meta.filtered[meta.filtered[ , j] == "Sensitive", ]) > 0 &
        nrow(meta.filtered[meta.filtered[ , j] == "Resistant", ]) > 0) {
      run_contrasts_global_phospho(list(c("Sensitive", "Resistant")), j, 'id', 
                                   meta.filtered, omics, beatAML, gmt.features, 
                                   gmt.features2, gmt.features3a, gmt.features3b, 
                                   gmt.drug, BeatAML.data$drug, temp.base.path,
                                   synapse_id = dataFolder) 
    }
  } 
}

# need to change names w "+" at the end because it gives this error during limma:
# Error in parse(text = e[j]) : <text>:2:0: unexpected end of input
# 1: CD14+-CD34+
#   ^
# run cell type contrasts without drug sensitivity filter
setwd(base.path)
sort.contrasts5 <- sort.contrasts4[3:length(sort.contrasts4)]
# skip CD14+ vs CD14+ Flow for now bc need to add checking of ksdb gmt coverage
# also skip CD14_Pos_Flow vs CD34_Pos for same reason
# also skip CD14+ Flow vs MSC Flow for same reason
# also skip CD14+ Flow vs CD34+ Flow for same reason
sort.contrasts6 <- list(c("CD14_Pos", "CD14_Pos_Flow"), 
                        c("CD14_Pos_Flow", "CD34_Pos"),
                        c("CD14_Pos_Flow", "MSC_Flow"),
                        c("CD14_Pos_Flow", "CD34_Pos_Flow"))
run_contrasts_global_phospho(sort.contrasts6, "SampleType4",'id', meta.df,
                             omics, beatAML, gmt.features, gmt.features2, 
                             gmt.features3a, gmt.features3b, 
                             gmt.drug, BeatAML.data$drug, base.path,
                             synapse_id, KSEA = FALSE)
for (j in drug.types) {
  j = "Ven"
  #j = "Aza.Ven"
  for (i in sens.types) {
    #i = "Resistant"
    # filter meta.df
    meta.filtered <- meta.df[meta.df[ , j] == i, ]
    
    # create folder for each drug sensitivity type
    drug.sens.type <- paste0(j, "_", i)
    setwd(base.path)
    dir.create(drug.sens.type)
    setwd(drug.sens.type)
    temp.base.path <- file.path(base.path, drug.sens.type)
    dataFolder <- 
      synapser::synStore(synapser::Folder(drug.sens.type,
                                          parent = synapse_id))
    
    # run cell type contrasts for each drug sensitivity
    run_contrasts_global_phospho(sort.contrasts4[1:length(sort.contrasts4)], "SampleType4", 'id', 
                                 meta.filtered, omics, beatAML, gmt.features, 
                                 gmt.features2, gmt.features3a, gmt.features3b, 
                                 gmt.drug, BeatAML.data$drug, temp.base.path,
                                 synapse_id = dataFolder)
    
    # no KSEA for Aza_Resistant sort.contrasts4 items 3, 6, 8, 10
    # no KSEA for Aza.Ven_Resistant sort.contrasts4 items 3, 6, 8
    # no KSEA for Ven_Resistant sort.contrasts4 items 3, 6, 8
    # not enough samples for Aza.Ven_Resistant sort.contrasts4 item 10 or Ven_Resistant sort.contrasts4 item 10:
    # Error in .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  : 
    #                    No residual degrees of freedom in linear model fits
    # run_contrasts_global_phospho(sort.contrasts4[8], "SampleType4", 'id',
    #                              meta.filtered, omics, beatAML, gmt.features,
    #                              gmt.features2, gmt.features3a, gmt.features3b,
    #                              gmt.drug, BeatAML.data$drug, temp.base.path,
    #                              synapse_id = dataFolder, KSEA = FALSE)
  } 
}

