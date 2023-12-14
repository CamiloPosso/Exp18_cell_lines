# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 21
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2023-12-14

library(readxl);library(panSEA); library(synapser)

#### 1. Import metadata & crosstabs ####
setwd(
"~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/"
)
meta.df <- readxl::read_excel("Experiment 21 - Patient Samples_F.xlsx", 
                              sheet = 2)
global.df <- read.table(
  "data/global_data/ptrc_ex21_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "data/phospho_data/ptrc_ex21_crosstab_phospho_SiteID_corrected.txt", 
  sep = "\t")

# add column for feature names and then make it the first column
global.df$Gene <- rownames(global.df)

#### 2. Import BeatAML data ####
# import drug MOA annotations
moa.BeatAML <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

# login to Synapse
synLogin()

# load data from Synapse
BeatAML_synapse_id <- c("BeatAML_drug_response_AUC.csv" = "syn51674470", 
                        "BeatAML_proteomics_metadata.txt" = "syn25807733",
                        "BeatAML_global_proteomics_corrected.txt" = "syn25714248",
                        "BeatAML_phospho_proteomics_corrected.txt" = "syn25714921")

for (i in 1:length(synapse_id_map)) {
  lapply(BeatAML_synapse_id, synapser::synGet, 
         downloadLocation = names(BeatAML_synapse_id))
}

drug.BeatAML <- read.csv(names(BeatAML_synapse_id)[1])
meta.BeatAML <- read.table(names(BeatAML_synapse_id)[2])
global.BeatAML <- read.table(names(BeatAML_synapse_id)[3])
phospho.BeatAML <- read.table(names(BeatAML_synapse_id)[4])

#### Step 4. Format BeatAML data for DMEA ####
## format drug sensitivity data frame
# format drug.BeatAML wide (samples in first column, drug names for rest of columns)
drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, value.var = "auc", fill = NA)

# remove drugs without moa annotations and drug combos
valid.drugs <- names(drug.BeatAML)[names(drug.BeatAML) %in% moa.BeatAML[!is.na(moa.BeatAML),]$Drug] # 167 drugs
drug.BeatAML <- drug.BeatAML[ , c("sample_id", valid.drugs)] # 167 drugs
moa.BeatAML <- moa.BeatAML[moa.BeatAML$Drug %in% names(drug.BeatAML)[2:ncol(drug.BeatAML)], ]

## format global proteomics data frame
# change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
global.ids <- names(global.BeatAML)

# remove X and any 0's from start of each column name
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

# transpose global.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
global.BeatAML <- as.data.frame(t(global.BeatAML))

# make first column Barcode.ID
global.BeatAML$Barcode.ID <- rownames(global.BeatAML)
global.BeatAML <- global.BeatAML[ , c("Barcode.ID", names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]

# make sure columns with sample IDs have same names between global.BeatAML and drug.BeatAML
names(drug.BeatAML)[1] <- names(global.BeatAML)[1]

# ## format phospho-proteomics data frame
# # change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
# phospho.ids <- names(phospho.BeatAML)
# 
# # remove X and any 0's from start of each column name
# # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
# for(i in seq_len(length(phospho.ids))){
#   phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
#   
#   if(substring(phospho.ids[i], 1, 1) == 0){
#     phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
#   }
#   
#   if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
#     phospho.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
#   }
# }
# 
# # replace phospho.BeatAML column names 
# names(phospho.BeatAML) <- phospho.ids
# 
# # transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
# phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))
# 
# # make first column Barcode.ID
# phospho.BeatAML$Barcode.ID <- rownames(phospho.BeatAML)
# phospho.BeatAML <- phospho.BeatAML[ , c("Barcode.ID", names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]
# 
# # make sure columns with sample IDs have same names between phospho.BeatAML and drug.BeatAML
# names(drug.BeatAML)[1] <- names(phospho.BeatAML)[1]

#### 2. Run panSEA across contrasts for each exp, omics type ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/data/")
omics <- c("global", "phospho")
for (k in 1:length(omics)) {
  ## prepare set annotations
  # generate gmt.features beforehand to save time
  if (omics[k] == "global") {
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
    ksdb <- read.csv("ksdb_human_20231101.csv")
    ksdb$SUB_SITE <- paste0(ksdb$SUBSTRATE, ksdb$SUB_MOD_RSD, collapse = "_")
    gmt <- DMEA::as_gmt(ksdb, "SUB_SITE", "KINASE", min.per.set = 6, 
                        descriptions = "KIN_ACC_ID")
  }
  
  ## prepare other input parameters
  # identify contrasts
  treatments <- na.omit(unique(meta.df$Treatment))
  contrasts <- c()
  for (i in 1:length(treatments)) {
    for (j in 1:length(treatments)) {
      # only pair different treatments
      if (i != j) {
        # put treatments in alphabetical order
        alpha.pair <- sort(c(treatments[i], treatments[j]))
        
        # keep track of all treatment contrasts
        contrasts <- c(contrasts,
                       paste(alpha.pair[1], "vs.", alpha.pair[2]))
      }
    }
  }
  contrasts <- unique(contrasts)
  
  # identify column names containing data for each treatment
  # add X in front of tube names to match column names 
  # after import (e.g., X5 for tube 5)
  sample.groups <- list()
  for (i in 1:length(treatments)) {
    sample.groups[[treatments[i]]] <- 
      paste0("X", na.omit(meta.df[meta.df$Treatment == treatments[i], ]$Tube))
  }
  
  
  # run panSEA for each omics type across all contrasts
  # CAUTION: this only works because the # of samples for each treatment type is equal
  # otherwise would have to run panSEA for each contrast separately and 
  # perhaps set the group.samples input parameter for panSEA
  
  # assemble inputs
  data.list <- list()
  expr.BeatAML <- list()
  gmt.features <- list()
  types <- contrasts
  for (i in 1:length(contrasts)) {
    # identify treatments in this contrast
    contrast.treatments <- stringr::str_split(contrasts[i], " vs. ")[[1]]
    
    # assemble each input data set for each contrast
    data.list[[contrasts[i]]] <- 
      global.df[ , c("Gene", sample.groups[[contrast.treatments[1]]],
                     sample.groups[[contrast.treatments[3]]])]
    
    expr.BeatAML[[contrasts[i]]] <- global.BeatAML
    
    gmt.features[[i]] <- gmt
  }
  
  # run panSEA by querying BeatAML data
  global.panSEA.BeatAML <- panSEA::panSEA(data.list, types, 
                                          gmt.features = gmt.features,
                                          gmt.drugs = DMEA::as_gmt(moa.BeatAML),
                                          drug.sensitivity = drug.BeatAML,
                                          expression = expr.BeatAML)
  
  ## upload results to Synapse
  saveRDS(global.panSEA.BeatAML, file="exp21_global_panSEA_BeatAML.rds")
  synapse_id_map <- c("syn51409382" = "data/global_data/", 
                      "syn51409382" = "data/phospho_data/")
  # instead of crosstabs, change to upload panSEA results
  for (i in 1:length(synapse_id_map)) {
    crosstabs <- list.files(path = synapse_id_map[i],
                            pattern = "crosstab.*.txt",
                            full.names = T)
    files <- lapply(crosstabs, synapser::File,
                    parent = names(synapse_id_map)[i])
    lapply(files, synapser::synStore)
  }
}
