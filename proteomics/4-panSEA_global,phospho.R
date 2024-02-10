# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 24 TMT
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-02-09

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)

# overview
#### 1. Import metadata & crosstabs ####
#### 2. Import BeatAML data ####
#### 3. Format BeatAML data for DMEA ####
#### 4. Run panSEA across omics for each contrast ####

# load analysis function
run_contrasts_global_phospho <- function(contrasts, contrast.type, meta.df,
                                         omics, beatAML, gmt.features, 
                                         gmt.features2, gmt.features3,
                                         gmt.drug, drug.BeatAML) {
  for (k in 1:length(contrasts)) {
    setwd(paste0(base.path, synapse_id_map[k]))
    
    # identify samples for each side of contrast
    c1 <- contrasts[[k]][1]
    c2 <- contrasts[[k]][2]
    group.names <- c(c1, c2)
    group.samples <- list(meta.df[meta.df[,contrast.type] == c1, ]$MeasurementName,
                          meta.df[meta.df[,contrast.type] == c2, ]$MeasurementName)
    
    # run panSEA across omics types
    # CAUTION: this only works because the # of samples for each treatment type 
    # is equal; otherwise would have to run panSEA for each contrast separately 
    # and perhaps set the group.samples input parameter for panSEA
    panSEA1 <- panSEA::panSEA(omics, c("global", "phospho"), 
                              feature.names = c("Gene", "SUB_SITE"), 
                              gmt.features = gmt.features,
                              gmt.drugs = gmt.drug,
                              drug.sensitivity = drug.BeatAML,
                              expression = beatAML, 
                              min.per.set = 6)
    
    panSEA2 <- panSEA::panSEA(omics, c("global", "phospho"), DMEA = FALSE,
                              feature.names = c("Gene", "SUB_SITE"), 
                              gmt.features = gmt.features2,
                              min.per.set = 6)
    
    phospho.results <- list("Kinase activity" = panSEA1$mGSEA.results$phospho$result,
                            "Substrate activity" = panSEA2$mGSEA.results$phospho$result)
    panSEA3 <- panSEA::panSEA(phospho.results, names(phospho.results), DMEA = FALSE,
                              feature.names = rep("Feature_set", 2), 
                              GSEA.rank.var = rep("NES", 2),
                              gmt.features = gmt.features3,
                              min.per.set = 6)
    
    ## save results & upload to Synapse
    # store all results locally
    dir.create("analysis")
    setwd("analysis")
    #saveRDS(panSEA.BeatAML, file=paste0("exp24_", omics[k], "_panSEA_BeatAML.rds")) # 8.5 GB for 7 contrasts
    panSEA.BeatAML <- readRDS(paste0("exp24_", omics[k], "_panSEA_BeatAML.rds"))
    
    # set file names
    DEG.files <- list("Differential_expression_results.csv" = 
                        panSEA.BeatAML$mDEG.results$compiled.results$results,
                      "Differential_expression_mean_results.csv" =
                        panSEA.BeatAML$mDEG.results$compiled.results$mean.results,
                      "Differential_expression_correlation_matrix.pdf" =
                        panSEA.BeatAML$mDEG.results$compiled.results$corr.matrix,
                      "Differential_expression_dot_plot.pdf" =
                        panSEA.BeatAML$mDEG.results$compiled.results$dot.plot)
    DMEA.files <- list("DMEA_results.csv" =
                         panSEA.BeatAML$mDMEA.results$compiled.results$results,
                       "DMEA_mean_results.csv" =
                         panSEA.BeatAML$mDMEA.results$compiled.results$mean.results,
                       "DMEA_correlation_matrix.pdf" =
                         panSEA.BeatAML$mDMEA.results$compiled.results$corr.matrix,
                       "DMEA_dot_plot.pdf" =
                         panSEA.BeatAML$mDMEA.results$compiled.results$dot.plot,
                       "DMEA_interactive_network.graph.html" =
                         panSEA.BeatAML$mDMEA.network$interactive)
    if (omics[k] == "phospho") {
      KSEA.files <- list("KSEA_results.csv" =
                           panSEA.BeatAML$mKSEA.results$compiled.results$results,
                         "KSEA_mean_results.csv" =
                           panSEA.BeatAML$mKSEA.results$compiled.results$mean.results,
                         "KSEA_correlation_matrix.pdf" =
                           panSEA.BeatAML$mKSEA.results$compiled.results$corr.matrix,
                         "KSEA_dot_plot.pdf" =
                           panSEA.BeatAML$mKSEA.results$compiled.results$dot.plot,
                         "KSEA_interactive_network.graph.html" =
                           panSEA.BeatAML$mKSEA.network$interactive)
      all.files <- list('Differential expression' = DEG.files,
                        'KSEA' = KSEA.files,
                        'DMEA' = DMEA.files)
    } else {
      GSEA.files <- list("GSEA_results.csv" =
                           panSEA.BeatAML$mGSEA.results$compiled.results$results,
                         "GSEA_mean_results.csv" =
                           panSEA.BeatAML$mGSEA.results$compiled.results$mean.results,
                         "GSEA_correlation_matrix.pdf" =
                           panSEA.BeatAML$mGSEA.results$compiled.results$corr.matrix,
                         "GSEA_dot_plot.pdf" =
                           panSEA.BeatAML$mGSEA.results$compiled.results$dot.plot,
                         "GSEA_interactive_network.graph.html" =
                           panSEA.BeatAML$mGSEA.network$interactive)
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
    panSEA.BeatAML <- NULL # make space to process next omics type
  }
}

#### 1. Import metadata & crosstabs ####
setwd(
"~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/"
)
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

#### 2. Import BeatAML data ####
# import drug MOA annotations
moa.BeatAML <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

# login to Synapse
#synLogin()

# load data from Synapse
BeatAML.path <- "BeatAML_DMEA_inputs"
BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                        "Ex10_metadata.txt" = "syn25807733",
                        "ptrc_ex10_crosstab_global_gene_corrected.txt" = "syn25714248",
                        "ptrc_ex10_crosstab_phospho_siteID_corrected(1).txt" = "syn25714921")
lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)

setwd(BeatAML.path)
drug.BeatAML <- read.csv(names(BeatAML_synapse_id)[1])
meta.BeatAML <- read.table(names(BeatAML_synapse_id)[2], sep = "\t", 
                           header = TRUE)
global.BeatAML <- read.table(names(BeatAML_synapse_id)[3], sep = "\t", 
                             header = TRUE)
phospho.BeatAML <- read.table(names(BeatAML_synapse_id)[4], sep = "\t", 
                              header = TRUE)

#### 3. Format BeatAML data for DMEA ####
sample.names <- "Barcode.ID"

## format drug sensitivity data frame
# format drug.BeatAML wide (samples in first column, drug names for rest of columns)
drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, 
                                value.var = "auc", fill = NA)

# remove drugs without moa annotations and drug combos
valid.drugs <- 
  names(drug.BeatAML)[names(drug.BeatAML) %in% 
                        moa.BeatAML[!is.na(moa.BeatAML),]$Drug] # 167 drugs
drug.BeatAML <- drug.BeatAML[ , c("sample_id", valid.drugs)] # 167 drugs
moa.BeatAML <- 
  moa.BeatAML[moa.BeatAML$Drug %in% names(drug.BeatAML)[2:ncol(drug.BeatAML)], ]

# change sample column name to match expression data
names(drug.BeatAML)[1] <- sample.names

## format global proteomics data frame
# change global.BeatAML column names from SampleID.abbrev to 
# Barcode.ID to match drug.BeatAML
global.ids <- names(global.BeatAML)

# remove X and any 0's from start of each column name and then
# replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
for(i in seq_len(length(global.ids))){
  global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  
  if(substring(global.ids[i], 1, 1) == 0){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  }
  
  if(global.ids[i] %in% meta.BeatAML$SampleID.abbrev){
    global.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == global.ids[i], ]$Barcode.ID
  }
}

# replace global.BeatAML column names 
names(global.BeatAML) <- global.ids

# transpose global.BeatAML so that first column is Barcode.ID and 
# rest of columns are gene symbols
global.BeatAML <- as.data.frame(t(global.BeatAML))

# make first column Barcode.ID
global.BeatAML[, sample.names] <- rownames(global.BeatAML)
global.BeatAML <- 
  global.BeatAML[ , c(sample.names, 
                      names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]

## format phospho-proteomics data frame
# change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
phospho.ids <- names(phospho.BeatAML)

# remove X and any 0's from start of each column name and then
# replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
for(i in seq_len(length(phospho.ids))){
  phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))

  if(substring(phospho.ids[i], 1, 1) == 0){
    phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
  }

  if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
    phospho.ids[i] <- meta.BeatAML[
      meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
  }
}

# replace phospho.BeatAML column names
names(phospho.BeatAML) <- phospho.ids

# transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))

# make first column Barcode.ID
phospho.BeatAML[, sample.names] <- rownames(phospho.BeatAML)
phospho.BeatAML <- phospho.BeatAML[ , c(sample.names, names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]

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
beatAML <- list("global" = global.BeatAML,
                "phospho" = phospho.BeatAML)

# prepare universal inputs
gmt.features <- list(gmt, gmt.ksdb)
gmt.features2 <- list(gmt.H, gmt.sub)
gmt.features3 <- list(gmt, gmt.H) # for phospho GSEA^2

for (j in drug.types) {
  for (i in sens.types) {
    # filter meta.df
    meta.filtered <- meta.df[meta.df[ , j] == i, ]
    
    # run cell type contrasts for each drug sensitivity
    run_contrasts_global_phospho(sort.contrasts, "SampleType", meta.filtered,
                                 omics, beatAML, gmt.features, gmt.features2, 
                                 gmt.features3, gmt.drug, drug.BeatAML)
  } 
}
# run cell type contrasts without drug sensitivity filter
run_contrasts_global_phospho(sort.contrasts, "SampleType", meta.df,
                             omics, beatAML, gmt.features, gmt.features2, 
                             gmt.features3, gmt.drug, drug.BeatAML)

for (i in sort.types) {
  # filter meta.df
  meta.filtered <- meta.df[meta.df$SampleType == i, ]
  
  # run sens vs. res contrasts for each cell type
  for (j in drug.types) {
    run_contrasts_global_phospho(c("Sensitive", "Resistant"), j, meta.filtered, 
                                 omics, beatAML, gmt.features, gmt.features2, 
                                 gmt.features3, gmt.drug, drug.BeatAML)
  }
}
# run sens vs. res contrasts without cell type filter
for (j in drug.types) {
  run_contrasts_global_phospho(c("Sensitive", "Resistant"), j, meta.df,
                               omics, beatAML, gmt.features, gmt.features2, 
                               gmt.features3, gmt.drug, drug.BeatAML)
}