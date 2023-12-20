# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 21
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2023-12-14

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)

#### 1. Import metadata & crosstabs ####
setwd(
"~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/data/"
)
meta.df <- readxl::read_excel("Experiment 21 - Patient Samples_F.xlsx", 
                              sheet = 2)
global.df <- read.table(
  "global_data/ptrc_ex21_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/ptrc_ex21_crosstab_phospho_SiteID_corrected.txt", 
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
synLogin()

# load data from Synapse
BeatAML.path <- "BeatAML_DMEA_inputs"
BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                        "Ex10_metadata.txt" = "syn25807733",
                        "ptrc_ex10_crosstab_global_gene_corrected.txt" = "syn25714248",
                        "ptrc_ex10_crosstab_phospho_siteID_corrected.txt" = "syn25714921")
lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)

setwd(BeatAML.path)
drug.BeatAML <- read.csv(names(BeatAML_synapse_id)[1])
meta.BeatAML <- read.table(names(BeatAML_synapse_id)[2], sep = "\t", 
                           header = TRUE)
global.BeatAML <- read.table(names(BeatAML_synapse_id)[3], sep = "\t", 
                             header = TRUE)
phospho.BeatAML <- read.table(names(BeatAML_synapse_id)[4], sep = "\t", 
                              header = TRUE)

#### Step 4. Format BeatAML data for DMEA ####
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
global.BeatAML$Barcode.ID <- rownames(global.BeatAML)
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
phospho.BeatAML$Barcode.ID <- rownames(phospho.BeatAML)
phospho.BeatAML <- phospho.BeatAML[ , c(sample.names, names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]

#### 2. Run panSEA across contrasts for each exp, omics type ####
# synapse IDs must match order of omics list
# synapse_id_map <- c("syn51409382" = "global_data/", 
#                     "syn51409382" = "phospho_data/")
synapse_id_map <- c("syn53180717" = "global_data/",
                    "syn53180718" = "phospho_data/")

## prepare other input parameters
# identify contrasts
treatments <- na.omit(unique(meta.df$Treatment))

# set short-hand names for treatments for better plots
# IMPORTANT: these must have the same order as treatments above
treatment.names <- c("Ctl", "ASO", "ASO + Gilt", 
                     "NRAS ASO", "NRAS ASO + Gilt")

# define types based on short-hand contrasts
types <- c(
  # compare everything to control
  "ASO vs. Ctl", "ASO + Gilt vs. Ctl", 
  "NRAS ASO vs. Ctl", "NRAS ASO + Gilt vs. Ctl",
  
  # compare to ASO
  "ASO + Gilt vs. ASO", "NRAS ASO vs. ASO", 
  
  # compare to ASO + Gilt
  "NRAS ASO + Gilt vs. ASO + Gilt"
)

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
#setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/data/")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/data/"
omics <- c("global", "phospho")
for (k in 1:length(omics)) {
  setwd(paste0(base.path, synapse_id_map[k]))
  
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
    # only create gmt the first time
    # ksdb <- read.csv(paste0("https://raw.githubusercontent.com/BelindaBGarana/",
    #                         "panSEA/shiny-app/data/ksdb_20231101.csv"))
    # ksdb.human <- ksdb[
    #   ksdb$KIN_ORGANISM == "human" & ksdb$SUB_ORGANISM == "human", ]
    # ksdb <- NULL # save space
    # # ksdb.human$SUB_SITE <- paste(ksdb.human$SUBSTRATE, ksdb.human$SUB_MOD_RSD,
    # #                               collapse = "_")
    # ksdb.human <- ksdb.human %>% unite("SUB_SITE", 
    #                                    c("SUBSTRATE", "SUB_MOD_RSD"),
    #                                    sep = "-", remove = FALSE)
    # gmt <- DMEA::as_gmt(ksdb.human, "SUB_SITE", "KINASE", min.per.set = 6,
    #                     descriptions = "KIN_ACC_ID")
    # ksdb.human <- NULL # save space
    # 
    # # save gmt for future analyses
    # saveRDS(gmt, "gmt_ksdb_human_20231101.rds")
    
    gmt <- readRDS("gmt_ksdb_human_20231101.rds")
  }
  
  # run panSEA for each omics type across all contrasts
  # CAUTION: this only works because the # of samples for each treatment type 
  # is equal; otherwise would have to run panSEA for each contrast separately 
  # and perhaps set the group.samples input parameter for panSEA
  
  # assemble inputs
  data.list <- list()
  expr.BeatAML <- list()
  gmt.features <- list()
  for (i in 1:length(contrasts)) {
    # identify treatments in this contrast
    contrast.treatments <- stringr::str_split(contrasts[i], " vs. ")[[1]]
    
    # assemble each input data set for each contrast
    if (omics[k] == "global"){
      data.list[[types[i]]] <- 
        global.df[ , c("Gene", sample.groups[[contrast.treatments[1]]],
                       sample.groups[[contrast.treatments[2]]])]
      
      expr.BeatAML[[types[i]]] <- global.BeatAML 
      
      feature.names <- rep("Gene", length(types))
    } else if (omics[k] == "phospho") {
      data.list[[types[i]]] <- 
        phospho.df[ , c("SUB_SITE", sample.groups[[contrast.treatments[1]]],
                       sample.groups[[contrast.treatments[2]]])]
      
      expr.BeatAML[[types[i]]] <- phospho.BeatAML
      
      feature.names <- rep("SUB_SITE", length(types))
    }
    
    gmt.features[[types[i]]] <- gmt
  }
  
  # run panSEA by querying BeatAML data
  panSEA.BeatAML <- panSEA::panSEA(data.list, types, 
                                   feature.names = feature.names, 
                                   gmt.features = gmt.features,
                                   gmt.drugs = DMEA::as_gmt(moa.BeatAML),
                                   drug.sensitivity = drug.BeatAML,
                                   expression = expr.BeatAML)
  
  
  # phospho didn't have enough kinase-substrate sets
  # Error in GSEA_custom(input, gmt, num.permutations, stat.type, min.per.set,  : 
  #                        annotations for 2+ drug sets are required
  
  ## save results & upload to Synapse
  # store all results locally
  dir.create("analysis")
  setwd("analysis")
  saveRDS(panSEA.BeatAML, file=paste0("exp21_", omics[k], "_panSEA_BeatAML.rds")) # 8.5 GB for 7 contrasts
  #panSEA.BeatAML <- readRDS(paste0("exp21_", omics[k], "_panSEA_BeatAML.rds"))
  
  # set file names
  DEG.files <- list("Differential_expression_results.csv" = 
                      panSEA.BeatAML$mDEG.results$compiled.results$results,
                    "Differential_expression_mean_results.csv" =
                      panSEA.BeatAML$mDEG.results$compiled.results$mean.results,
                    # "Differential_expression_venn_diagram.pdf" =
                    #   panSEA.BeatAML$mDEG.results$compiled.results$venn.diagram,
                    "Differential_expression_correlation_matrix.pdf" =
                      panSEA.BeatAML$mDEG.results$compiled.results$corr.plot,
                    "Differential_expression_dot_plot.pdf" =
                      panSEA.BeatAML$mDEG.results$compiled.results$dot.plot)
  GSEA.files <- list("GSEA_results.csv" =
                       panSEA.BeatAML$mGSEA.results$compiled.results$results,
                     "GSEA_mean_results.csv" =
                       panSEA.BeatAML$mGSEA.results$compiled.results$mean.results,
                     # "GSEA_venn_diagram.pdf" =
                     #   panSEA.BeatAML$mGSEA.results$compiled.results$venn.diagram,
                     "GSEA_correlation_matrix.pdf" =
                       panSEA.BeatAML$mGSEA.results$compiled.results$corr.plot,
                     "GSEA_dot_plot.pdf" =
                       panSEA.BeatAML$mGSEA.results$compiled.results$dot.plot,
                     "GSEA_interactive_network.graph.html" =
                       panSEA.BeatAML$mGSEA.network$interactive)
  DMEA.files <- list("DMEA_results.csv" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$results,
                     "DMEA_mean_results.csv" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$mean.results,
                     # "DMEA_venn_diagram.pdf" =
                     #   panSEA.BeatAML$mDMEA.results$compiled.results$venn.diagram,
                     "DMEA_correlation_matrix.pdf" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$corr.matrix,
                     "DMEA_dot_plot.pdf" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$dot.plot,
                     "DMEA_interactive_network.graph.html" =
                       panSEA.BeatAML$mDMEA.network$interactive)
  all.files <- list('Differential expression' = DEG.files,
                    'GSEA' = GSEA.files,
                    'DMEA' = DMEA.files)
  
  # instead of crosstabs, change to upload panSEA results
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
    dataFolder <- synapser::synStore(synapser::Folder(names(all.files)[i],
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
