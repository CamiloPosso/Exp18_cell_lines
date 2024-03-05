# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 24 TMT
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-02-09

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)
library(dplyr)

## get list of top n mountain plots (if any)
# input: EA output (e.g., mGSEA.result or mDMEA.result)
# output: list of top n mountain plots
get_top_mtn_plots <- function(base.result, n.top = 10, EA.type = "GSEA", 
                              sets = "Feature_set") {
  all.mtn <- base.result$mtn.plots
  temp.results <- base.result$result
  if (length(all.mtn) > 0) {
    if (length(all.mtn) > n.top) {
      # identify top significant enrichments
      mtn.results <- temp.results[temp.results[ , sets] %in% names(all.mtn), ]
      top.mtn.results <- mtn.results %>% slice_max(abs(NES), n = n.top)
      all.top.mtn <- all.mtn[names(all.mtn) %in% top.mtn.results[ , sets]]
    } else {
      all.top.mtn <- all.mtn
    }
    
    # create list of mtn plots for files
    mtn.file.names <- c()
    for (i in 1:length(all.top.mtn)) {
      mtn.file.names <- c(mtn.file.names, 
                          paste0(EA.type, "_mtn_plot_", names(all.top.mtn)[i], ".pdf"))
    }
    names(all.top.mtn) <- mtn.file.names
  } else {
    all.top.mtn <- list()
  }
  return(all.top.mtn)
}

save_to_synapse <- function(temp.files, resultsFolder) {
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
  if (length(CSV.files) > 0) {
    # save locally
    for (j in 1:length(CSV.files)) {
      write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
    }
    
    # save to synapse
    CSVs <- lapply(as.list(CSV.files), synapser::File,
                   parent = resultsFolder)
    lapply(CSVs, synapser::synStore)
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  if (length(PDF.files) > 0) {
    # save locally
    for (j in 1:length(PDF.files)) {
      if (is.list(temp.files[[PDF.files[j]]])) {
        ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                        device = "pdf") 
      }
    }
    
    # save to synapse
    PDF.files <- list.files(pattern = ".*.pdf", full.names = TRUE)
    PDFs <- lapply(as.list(PDF.files), synapser::File,
                   parent = resultsFolder)
    lapply(PDFs, synapser::synStore)
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    # save locally
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j])
    }
    
    # save to synapse
    HTMLs <- lapply(HTML.files, synapser::File,
                    parent = resultsFolder)
    lapply(HTMLs, synapser::synStore) 
  }
}

## run panSEA across global & phospho types for multiple contrasts
# input: vector of contrasts, contrast.type column for meta.df extraction, 
#   omics data.list, beatAML expression list, gmt set information, 
#   BeatAML drug AUC, base.path, Synapse ID info
# output: panSEA result files are saved locally & uploaded to Synapse
run_contrasts_global_phospho <- function(contrasts, contrast.type, 
                                         id.type, meta.df,
                                         omics, beatAML, gmt.features, 
                                         gmt.features2, gmt.features3a,
                                         gmt.features3b,
                                         gmt.drug, drug.BeatAML, 
                                         base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                         synapse_id, KSEA = TRUE) {
  # check if sufficient coverage for enrichment analyses
  # get KEGG info
  msigdb.KEGG <- msigdbr::msigdbr("Homo sapiens", "C2", "CP")
  msigdb.KEGG <- as.data.frame(msigdb.info[, c(
    "gene_symbol",
    "gs_name",
    "gs_description"
  )])
  
  msigdb.H <- msigdbr::msigdbr("Homo sapiens", "H")
  msigdb.H <- as.data.frame(msigdb.info[, c(
    "gene_symbol",
    "gs_name",
    "gs_description"
  )])
  
  global.KEGG <- msigdb.KEGG[msigdb.KEGG$gene_symbol %in% global.df$Gene, ]
  global.H <- msigdb.H[msigdb.H$gene_symbol %in% global.df$Gene, ]
  
  if (KSEA) {
    library(plyr); library(dplyr); library(tidyr)
    ksdb <- read.csv(paste0("https://raw.githubusercontent.com/BelindaBGarana/",
                            "panSEA/shiny-app/data/ksdb_20231101.csv"))
    ksdb.human <- ksdb[
      ksdb$KIN_ORGANISM == "human" & ksdb$SUB_ORGANISM == "human", ]
    ksdb <- NULL # save space
    # ksdb.human$SUB_SITE <- paste(ksdb.human$SUBSTRATE, ksdb.human$SUB_MOD_RSD,
    #                               collapse = "_")
    ksdb.human$SUB_MOD_RSD_PNNL <- ksdb.human$SUB_MOD_RSD
    ksdb.human$SUB_MOD_RSD_PNNL <- paste0(ksdb.human$SUB_MOD_RSD_PNNL, tolower(substr(ksdb.human$SUB_MOD_RSD_PNNL,1,1)))
    ksdb.human <- ksdb.human %>% unite("SUB_SITE",
                                       c("SUBSTRATE", "SUB_MOD_RSD_PNNL"),
                                       sep = "-", remove = FALSE)
    phospho.ksdb <- ksdb.human[ksdb.human$SUB_SITE %in% phospho.df$SUB_SITE, ]
  } else {
    phospho.ksdb <- data.frame()
  }
  
  if (nrow(global.KEGG) > 0) {
    gmt.global.KEGG <- DMEA::as_gmt(
      global.KEGG, "gene_symbol", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    )
  } else {
    gmt.global.KEGG <- list("genesets" = list())
  }
  
  if (nrow(global.H) > 0) {
    gmt.global.H <- DMEA::as_gmt(
      global.H, "gene_symbol", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    )
  } else {
    gmt.global.H <- list("genesets" = list())
  }
  
  if (nrow(phospho.ksdb) > 0) {
    gmt.phospho.ksdb <- DMEA::as_gmt(
      phospho.ksdb, "SUB_SITE", "KINASE", min.per.set = 6,
      descriptions = "KIN_ACC_ID"
    )
  } else {
    gmt.phospho.ksdb <- list("genesets" = list())
  }
  
  for (k in 1:length(contrasts)) {
    # identify samples for each side of contrast
    c1 <- contrasts[[k]][1]
    c2 <- contrasts[[k]][2]
    group.names <- c(c1, c2)
    contrast.name <- paste0(contrast.type, "_", contrasts[[k]][1], "_vs_", contrasts[[k]][2])
    group.samples <- list(meta.df[meta.df[,contrast.type] == c1, id.type],
                          meta.df[meta.df[,contrast.type] == c2, id.type])
    
    if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0) {
      # run panSEA across omics types
      # CAUTION: this only works because the # of samples for each treatment type 
      # is equal; otherwise would have to run panSEA for each contrast separately 
      # and perhaps set the group.samples input parameter for panSEA
      
      # run global against KEGG sets, phospho against kinase sets from ksdb
      if (length(gmt.global.KEGG$genesets) > 1 &
          length(gmt.phospho.ksdb$genesets) > 1 & KSEA) {
        panSEA1 <- panSEA::panSEA(omics, c("global", "phospho"), 
                                  feature.names = c("Gene", "SUB_SITE"), 
                                  group.names = group.names,
                                  group.samples = group.samples,
                                  gmt.features = gmt.features,
                                  gmt.drugs = gmt.drug,
                                  drug.sensitivity = drug.BeatAML,
                                  expression = beatAML, 
                                  min.per.set = 6)
        
        
      } else if (length(gmt.global.KEGG$genesets) > 1) {
        panSEA1 <- panSEA::panSEA(omics, c("global", "phospho"), 
                                  GSEA = FALSE,
                                  feature.names = c("Gene", "SUB_SITE"), 
                                  group.names = group.names,
                                  group.samples = group.samples,
                                  gmt.features = gmt.features,
                                  gmt.drugs = gmt.drug,
                                  drug.sensitivity = drug.BeatAML,
                                  expression = beatAML, 
                                  min.per.set = 6)
        
        global1.input <- list("global" = panSEA1$mDEG.results$all.results$global)
        global1 <- panSEA::panSEA(global1.input, c("global"), 
                                  feature.names = c("Gene"), DMEA = FALSE,
                                  group.names = contrast.name,
                                  group.samples = 2, # "Log2FC" is in column 2
                                  gmt.features = list(gmt.features[[1]]),
                                  min.per.set = 6)
        
        # add phospho data to global
        panSEA1$mGSEA.results <- list("compiled.results" = NA, 
                                      "all.results" = list(
                                        "global" = global1$mGSEA.results$Log2FC$all.results$global,
                                        "phospho" = NA))
      } else if (length(gmt.phospho.ksdb$genesets) > 1 & KSEA) {
        panSEA1 <- panSEA::panSEA(omics, c("global", "phospho"), 
                                  GSEA = FALSE,
                                  feature.names = c("Gene", "SUB_SITE"), 
                                  group.names = group.names,
                                  group.samples = group.samples,
                                  gmt.features = gmt.features,
                                  gmt.drugs = gmt.drug,
                                  drug.sensitivity = drug.BeatAML,
                                  expression = beatAML, 
                                  min.per.set = 6)
        
        phospho1.input <- list("phospho" = panSEA1$mDEG.results$all.results$phospho)
        phospho1 <- panSEA::panSEA(phospho1.input, c("phospho"), 
                                   feature.names = c("SUB_SITE"), DMEA = FALSE,
                                   group.names = contrast.name,
                                   group.samples = 2, # "Log2FC" is in column 2
                                   gmt.features = list(gmt.features[[2]]), 
                                   min.per.set = 6)
        
        # add phospho data to global
        panSEA1$mGSEA.results <- list("compiled.results" = NA, 
                                      "all.results" = list(
                                        "global" = NA,
                                        "phospho" = phospho1$mGSEA.results$Log2FC$all.results$global))
      } else {
        panSEA1 <- NULL
      }
      
      
      # use DEG results to avoid re-running DEG analysis
      # run global against hallmark sets, phospho against substrate sets
      deg.results <- list("global" = panSEA1$mDEG.results$all.results$global,
                          "phospho" = panSEA1$mDEG.results$all.results$phospho)
      
      # check if substrate sets are sufficiently covered
      phospho.degs <- panSEA1$mDEG.results$all.results$phospho
      phospho.degs <- phospho.degs %>% tidyr::extract(SUB_SITE, "KINASE", remove = FALSE)
      phospho.degs <- na.omit(phospho.degs[phospho.degs$Log2FC != 0, ])
      gmt.sub <- DMEA::as_gmt(phospho.degs, "SUB_SITE", "KINASE", min.per.set = 6)
      
      if (length(gmt.sub$genesets) > 1 & length(gmt.global.H$genesets) > 1) {
        panSEA2 <- panSEA::panSEA(deg.results, names(deg.results), DMEA = FALSE,
                                  group.names = contrast.name,
                                  group.samples = 2, # "Log2FC" is in column 2
                                  feature.names = c("Gene", "SUB_SITE"), 
                                  gmt.features = gmt.features2,
                                  min.per.set = 6) 
        sub.results <- panSEA2$mGSEA.results$Log2FC$all.results$phospho$result
      } else if (length(gmt.global.H$genesets) > 1) {
        panSEA2 <- panSEA::panSEA(list("global" = deg.results[[1]]), "global", 
                                  DMEA = FALSE,
                                  group.names = contrast.name,
                                  group.samples = 2, # "Log2FC" is in column 2
                                  feature.names = c("Gene"), 
                                  gmt.features = list(gmt.features2[[1]]),
                                  min.per.set = 6) 
        
        sub.results <- data.frame()
      } else if (length(gmt.sub$genesets) > 1) {
        panSEA2 <- panSEA::panSEA(list("phospho" = deg.results[[2]]), "phospho", 
                                  DMEA = FALSE,
                                  group.names = contrast.name,
                                  group.samples = 2, # "Log2FC" is in column 2
                                  feature.names = c("SUB_SITE"), 
                                  gmt.features = list(gmt.features2[[2]]),
                                  min.per.set = 6) 
        
        sub.results <- panSEA2$mGSEA.results$Log2FC$all.results$phospho$result
      } else {
        panSEA2 <- NULL
        sub.results <- data.frame()
      }
      
      if (length(gmt.phospho.ksdb$genesets) > 1) {
        ksea.results <- panSEA1$mGSEA.results$all.results$phospho$result
      } else {
        ksea.results <- data.frame()
      }
      
      # check if 2+ gene sets are covered by KSEA or substrate enrichment results
      # merge with kinase, substrate results
      # create gmts to determine if we can do further GSEA
      if (nrow(ksea.results) > 0) {
        kin.KEGG <- msigdb.KEGG[msigdb.KEGG$gene_symbol %in% ksea.results$Feature_set, ]
        if (nrow(kin.KEGG) > 0) {
          gmt.kin.KEGG <- DMEA::as_gmt(
            kin.KEGG, "gene_symbol", "gs_name", min.per.set = 6,
            descriptions = "gs_description"
          )  
        } else {
          gmt.kin.KEGG <- list('genesets' = list())
        }
        
        kin.H <- msigdb.H[msigdb.H$gene_symbol %in% ksea.results$Feature_set, ]
        if (nrow(kin.H) > 0) {
          gmt.kin.H <- DMEA::as_gmt(
            kin.H, "gene_symbol", "gs_name", min.per.set = 6,
            descriptions = "gs_description"
          )  
        } else {
          gmt.kin.H <- list('genesets' = list())
        }
      } else {
        gmt.kin.KEGG <- list('genesets' = list())
        gmt.kin.H <- list('genesets' = list())
      }
      
      if (nrow(sub.results) > 0) {
        sub.KEGG <- msigdb.KEGG[msigdb.KEGG$gene_symbol %in% sub.results$Feature_set, ]
        if (nrow(sub.KEGG) > 0) {
          gmt.sub.KEGG <- DMEA::as_gmt(
            sub.KEGG, "gene_symbol", "gs_name", min.per.set = 6,
            descriptions = "gs_description" 
          )
        } else {
          gmt.sub.KEGG <- list('genesets' = list())
        }
        
        sub.H <- msigdb.H[msigdb.H$gene_symbol %in% sub.results$Feature_set, ]
        if (nrow(sub.H) > 0) {
          gmt.sub.H <- DMEA::as_gmt(
            sub.H, "gene_symbol", "gs_name", min.per.set = 6,
            descriptions = "gs_description"
          )   
        } else {
          gmt.sub.H <- list('genesets' = list())
        }
      } else {
        gmt.sub.KEGG <- list('genesets' = list())
        gmt.sub.H <- list('genesets' = list())
      }
      
      if (length(gmt.kin.KEGG$genesets) > 1 & length(gmt.sub.KEGG$genesets) > 1) {
        phospho.a <- list("KSEA" = ksea.results,
                          "Substrate enrichment" = sub.results)
      } else if (length(gmt.sub.KEGG$genesets) > 1) {
        phospho.a <- list("Substrate enrichment" = sub.results)
      } else if (length(gmt.kin.KEGG$genesets) > 1) {
        phospho.a <- list("KSEA" = ksea.reuslts)
      } else {
        phospho.a <- NULL
      }
      
      if (length(gmt.kin.H$genesets) > 1 & length(gmt.sub.H$genesets) > 1) {
        phospho.b <- list("KSEA" = ksea.results,
                          "Substrate enrichment" = sub.results)
      } else if (length(gmt.sub.H$genesets) > 1) {
        phospho.b <- list("Substrate enrichment" = sub.results)
      } else if (length(gmt.kin.H$genesets) > 1) {
        phospho.b <- list("KSEA" = ksea.reuslts)
      } else {
        phospho.b <- NULL
      } 
      
      if (is.null(phospho.a) & is.null(phospho.b)) {
        panSEA3a <- NULL
        panSEA3b <- NULL
      } else if (is.null(phospho.a)) {
        panSEA3a <- NULL
        
        # run kinase and substrate results against hallmark sets
        panSEA3b <- 
          panSEA::panSEA(phospho.b, names(phospho.b), DMEA = FALSE,
                         group.names = contrast.name,
                         group.samples = 4, # "NES" is in column 4
                         feature.names = rep("Feature_set", 
                                             length(phospho.b)), 
                         GSEA.rank.var = rep("NES", length(phospho.b)),
                         gmt.features = gmt.features3b[1:length(phospho.b)],
                         min.per.set = 6) 
      } else if (is.null(phospho.b)) {
        # run kinase and substrate results against KEGG sets
        panSEA3a <- 
          panSEA::panSEA(phospho.a, names(phospho.a), DMEA = FALSE,
                         group.names = contrast.name,
                         group.samples = 4, # "NES" is in column 4
                         feature.names = rep("Feature_set", 
                                             length(phospho.a)), 
                         GSEA.rank.var = rep("NES", length(phospho.a)),
                         gmt.features = gmt.features3a[1:length(phospho.a)],
                         min.per.set = 6)
        
        panSEA3b <- NULL
      } else {
        # run kinase and substrate results against KEGG sets
        panSEA3a <- 
          panSEA::panSEA(phospho.a, names(phospho.a), DMEA = FALSE,
                         group.names = contrast.name,
                         group.samples = 4, # "NES" is in column 4
                         feature.names = rep("Feature_set", 
                                             length(phospho.a)), 
                         GSEA.rank.var = rep("NES", length(phospho.a)),
                         gmt.features = gmt.features3a[1:length(phospho.a)],
                         min.per.set = 6)
        
        # run kinase and substrate results against hallmark sets
        panSEA3b <- 
          panSEA::panSEA(phospho.b, names(phospho.b), DMEA = FALSE,
                         group.names = contrast.name,
                         group.samples = 4, # "NES" is in column 4
                         feature.names = rep("Feature_set", 
                                             length(phospho.b)), 
                         GSEA.rank.var = rep("NES", length(phospho.b)),
                         gmt.features = gmt.features3b[1:length(phospho.b)],
                         min.per.set = 6) 
      }
      
      ## save results & upload to Synapse
      # get GSEA mtn plots
      if (is.null(panSEA3a) & is.null(panSEA3b)) {
        phospho.KSEA.KEGG.files <- list()
        phospho.KSEA.hallmark.files <- list()
        
        phospho.sub.KEGG.files <- list()
        phospho.sub.hallmark.files <- list()
      } else if (is.null(panSEA3a)) {
        kin.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$NES$all.results$'KSEA',
                                       EA.type = "KSEA_GSEA_hallmark")
        phospho.KSEA.KEGG.files <- list()
        phospho.KSEA.hallmark.files <- list("KSEA_GSEA_hallmark_results.csv" =
                                              panSEA3b$mGSEA.results$NES$all.results$phospho$result,
                                            "KSEA_GSEA_hallmark_volcano_plot.pdf" =
                                              panSEA3b$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                            "mtn_plots" = kin.mtn3b)
        
        sub.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$NES$all.results$'Substrate enrichment',
                                       EA.type = "Substrate_enrichment_GSEA_hallmark")
        phospho.sub.KEGG.files <- list()
        phospho.sub.hallmark.files <- list("Substrate_enrichment_GSEA_hallmark_results.csv" =
                                             panSEA3b$mGSEA.results$NES$all.results$phospho$result,
                                           "Substrate_enrichment_GSEA_hallmark_volcano_plot.pdf" =
                                             panSEA3b$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                           "mtn_plots" = sub.mtn3b)
      } else if (is.null(panSEA3b)) {
        kin.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$NES$all.results$'KSEA', 
                                       EA.type = "KSEA_GSEA_KEGG")
        phospho.KSEA.KEGG.files <- list("KSEA_GSEA_KEGG_results.csv" =
                                          panSEA3a$mGSEA.results$NES$all.results$phospho$result,
                                        "KSEA_GSEA_KEGG_volcano_plot.pdf" =
                                          panSEA3a$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                        "mtn_plots" = kin.mtn3a)
        phospho.KSEA.hallmark.files <- list()
        
        sub.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$NES$all.results$'Substrate enrichment',
                                       EA.type = "Substrate_enrichment_GSEA_KEGG")
        phospho.sub.KEGG.files <- list("Substrate_enrichment_GSEA_KEGG_results.csv" =
                                         panSEA3a$mGSEA.results$NES$all.results$phospho$result,
                                       "Substrate_enrichment_GSEA_KEGG_volcano_plot.pdf" =
                                         panSEA3a$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                       "mtn_plots" = sub.mtn3a)
        phospho.sub.hallmark.files <- list()
      } else {
        kin.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$NES$all.results$'KSEA', 
                                       EA.type = "KSEA_GSEA_KEGG")
        kin.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$NES$all.results$'KSEA',
                                       EA.type = "KSEA_GSEA_hallmark")
        phospho.KSEA.KEGG.files <- list("KSEA_GSEA_KEGG_results.csv" =
                                          panSEA3a$mGSEA.results$NES$all.results$phospho$result,
                                        "KSEA_GSEA_KEGG_volcano_plot.pdf" =
                                          panSEA3a$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                        "mtn_plots" = kin.mtn3a)
        phospho.KSEA.hallmark.files <- list("KSEA_GSEA_hallmark_results.csv" =
                                              panSEA3b$mGSEA.results$NES$all.results$phospho$result,
                                            "KSEA_GSEA_hallmark_volcano_plot.pdf" =
                                              panSEA3b$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                            "mtn_plots" = kin.mtn3b)
        
        sub.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$NES$all.results$'Substrate enrichment',
                                       EA.type = "Substrate_enrichment_GSEA_KEGG")
        sub.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$NES$all.results$'Substrate enrichment',
                                       EA.type = "Substrate_enrichment_GSEA_hallmark")
        phospho.sub.KEGG.files <- list("Substrate_enrichment_GSEA_KEGG_results.csv" =
                                         panSEA3a$mGSEA.results$NES$all.results$phospho$result,
                                       "Substrate_enrichment_GSEA_KEGG_volcano_plot.pdf" =
                                         panSEA3a$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                       "mtn_plots" = sub.mtn3a)
        phospho.sub.hallmark.files <- list("Substrate_enrichment_GSEA_hallmark_results.csv" =
                                             panSEA3b$mGSEA.results$NES$all.results$phospho$result,
                                           "Substrate_enrichment_GSEA_hallmark_volcano_plot.pdf" =
                                             panSEA3b$mGSEA.results$NES$all.results$phospho$volcano.plot,
                                           "mtn_plots" = sub.mtn3b)
      }
      
      # set file names
      if (!is.null(panSEA1)) {
        global.DEG.files <- list("Differential_expression_results.csv" = 
                                   panSEA1$mDEG.results$all.results$global)
        phospho.DEG.files <- list("Differential_expression_results.csv" = 
                                    panSEA1$mDEG.results$all.results$phospho)
        
        combo.DMEA.files <- list("DMEA_results.csv" =
                                   panSEA1$mDMEA.results$compiled.results$results,
                                 "DMEA_mean_results.csv" =
                                   panSEA1$mDMEA.results$compiled.results$mean.results,
                                 "DMEA_correlation_matrix.pdf" =
                                   panSEA1$mDMEA.results$compiled.results$corr.matrix,
                                 "DMEA_dot_plot.pdf" =
                                   panSEA1$mDMEA.results$compiled.results$dot.plot)
        
        # get DMEA mtn plots
        DMEA.global.mtn <- get_top_mtn_plots(panSEA1$mDMEA.results$all.results$global,
                                             EA.type = "DMEA")
        DMEA.phospho.mtn <- get_top_mtn_plots(panSEA1$mDMEA.results$all.results$phospho,
                                              EA.type = "DMEA")
        global.DMEA.files <- list("DMEA_results.csv" =
                                    panSEA1$mDMEA.results$all.results$global$result,
                                  "DMEA_correlation_results.csv" = 
                                    panSEA1$mDMEA.results$all.results$global$corr.result,
                                  "DMEA_correlation_scatter_plots.pdf" = 
                                    panSEA1$mDMEA.results$all.results$global$corr.scatter.plots,
                                  "DMEA_volcano_plot.pdf" =
                                    panSEA1$mDMEA.results$all.results$global$volcano.plot,
                                  "mtn_plots" = DMEA.global.mtn)
        phospho.DMEA.files <- list("DMEA_results.csv" =
                                     panSEA1$mDMEA.results$all.results$phospho$result,
                                   "DMEA_correlation_results.csv" = 
                                     panSEA1$mDMEA.results$all.results$phospho$corr.result,
                                   "DMEA_correlation_scatter_plots.pdf" = 
                                     panSEA1$mDMEA.results$all.results$phospho$corr.scatter.plots,
                                   "DMEA_volcano_plot.pdf" =
                                     panSEA1$mDMEA.results$all.results$phospho$volcano.plot,
                                   "mtn_plots" = DMEA.phospho.mtn)
        
      }
      
      if (length(gmt.global.KEGG$genesets) > 1) {
        global.mtn1 <- get_top_mtn_plots(panSEA1$mGSEA.results$all.results$global,
                                         EA.type = "GSEA_KEGG")
        global.GSEA.KEGG.files <- list("GSEA_KEGG_results.csv" =
                                         panSEA1$mGSEA.results$all.results$global$result,
                                       "GSEA_KEGG_volcano_plot.pdf" =
                                         panSEA1$mGSEA.results$all.results$global$volcano.plot,
                                       "mtn_plots" = global.mtn1) 
      } else {
        global.GSEA.KEGG.files <- list()
      }
      
      if (KSEA & length(gmt.phospho.ksdb$genesets) > 1) {
        phospho.mtn1 <- get_top_mtn_plots(panSEA1$mGSEA.results$all.results$phospho,
                                          EA.type = "KSEA")
        
        kin.files <- list('GSEA_KEGG' = phospho.KSEA.KEGG.files,
                          'GSEA_hallmark' = phospho.KSEA.hallmark.files)
        phospho.KSEA.files <- list("KSEA_results.csv" =
                                     panSEA1$mGSEA.results$all.results$phospho$result,
                                   "KSEA_volcano_plot.pdf" =
                                     panSEA1$mGSEA.results$all.results$phospho$volcano.plot,
                                   "mtn_plots" = phospho.mtn1,
                                   "GSEA" = kin.files) 
      } else {
        phospho.KSEA.files <- list()
      }
      
      if (length(gmt.global.H$genesets) > 1) {
        global.mtn2 <- get_top_mtn_plots(panSEA2$mGSEA.results$Log2FC$all.results$global,
                                         EA.type = "GSEA_hallmark")
        global.GSEA.hallmark.files <- list("GSEA_hallmark_results.csv" =
                                             panSEA2$mGSEA.results$Log2FC$all.results$global$result,
                                           "GSEA_hallmark_volcano_plot.pdf" =
                                             panSEA2$mGSEA.results$Log2FC$all.results$global$volcano.plot,
                                           "mtn_plots" = global.mtn2)
      }
      
      if (length(gmt.sub$genesets) > 1) {
        phospho.mtn2 <- get_top_mtn_plots(panSEA2$mGSEA.results$Log2FC$all.results$phospho,
                                          EA.type = "Substrate_enrichment")
        sub.files <- list('GSEA_KEGG' = phospho.sub.KEGG.files,
                          'GSEA_hallmark' = phospho.sub.hallmark.files)
        phospho.sub.files <- list("Substrate_enrichment_results.csv" =
                                    panSEA2$mGSEA.results$Log2FC$all.results$phospho$result,
                                  "Substrate_enrichment_volcano_plot.pdf" =
                                    panSEA2$mGSEA.results$Log2FC$all.results$phospho$volcano.plot,
                                  "mtn_plots" = phospho.mtn2,
                                  "GSEA" = sub.files) 
      } else {
        phospho.sub.files <- list()
      }
      
      global.files <- list('Differential_expression' = global.DEG.files, 
                           'DMEA' = global.DMEA.files,
                           'GSEA_KEGG' = global.GSEA.hallmark.files,
                           'GSEA_hallmark' = global.GSEA.KEGG.files)
      phospho.files <- list('Differential_expression' = phospho.DEG.files, 
                            'DMEA' = phospho.DMEA.files,
                            'KSEA' = phospho.KSEA.files,
                            'Substrate_enrichment' = phospho.sub.files)
      combo.files <- list('DMEA' = combo.DMEA.files)
      all.files <- list('global_and_phospho' = combo.files,
                        'global' = global.files,
                        'phospho' = phospho.files)
      
      # create folder for contrast
      setwd(base.path)
      dir.create(contrast.name)
      setwd(contrast.name)
      contrastFolder <- 
        synapser::synStore(synapser::Folder(contrast.name,
                                            parent = synapse_id))
      
      for (i in 1:length(all.files)) {
        # create folder for omics type (e.g., global)
        setwd(file.path(base.path, contrast.name))
        dir.create(names(all.files)[i])
        setwd(names(all.files)[i])
        omics.files <- all.files[[i]]
        omicsFolder <- 
          synapser::synStore(synapser::Folder(names(all.files)[i],
                                              parent = contrastFolder))
        for (j in 1:length(omics.files)) {
          # create folder for result type (e.g., Differential expression)
          setwd(file.path(base.path, contrast.name, names(all.files)[i]))
          dir.create(names(omics.files)[j])
          setwd(names(omics.files)[j])
          temp.files <- omics.files[[j]]
          resultsFolder <- 
            synapser::synStore(synapser::Folder(names(omics.files)[j],
                                                parent = omicsFolder))
          
          # save to synapse
          save_to_synapse(temp.files, resultsFolder)
          
          # save mtn plots if relevant
          mtn.files <- names(temp.files)[grepl("mtn_plots", names(temp.files))]
          if (length(mtn.files) > 0) {
            temp.mtn.files <- temp.files[["mtn_plots"]]
            if (length(temp.mtn.files) > 0) {
              # create folder for mtn plots
              dir.create("mtn_plots")
              setwd("mtn_plots")
              mtnFolder <- 
                synapser::synStore(synapser::Folder("mountain_plots",
                                                    parent = resultsFolder))
              
              # save mtn plots locally
              for (m in 1:length(temp.mtn.files)) {
                ggplot2::ggsave(names(temp.mtn.files)[m], temp.mtn.files[[m]], 
                                device = "pdf")
              } 
              
              # upload mtn plots to synapse
              PDFs <- lapply(as.list(names(temp.mtn.files)), synapser::File,
                             parent = mtnFolder)
              lapply(PDFs, synapser::synStore)
            }
          }
          setwd(file.path(base.path, contrast.name, names(all.files[i]), names(omics.files)[j]))
          
          # save KSEA_GSEA or sub_GSEA if relevant
          gsea.files <- names(temp.files)[grepl("GSEA", names(temp.files))]
          if (length(gsea.files) > 0) {
            temp.gsea.files <- temp.files[["GSEA"]]
            if (length(temp.gsea.files) > 0) {
              for (m in 1:length(temp.gsea.files)) {
                # create folder for mtn plots
                dir.create(names(temp.gsea.files)[m])
                setwd(names(temp.gsea.files)[m])
                gseaFolder <- 
                  synapser::synStore(synapser::Folder(paste0(names(omics.files)[j], "_", names(temp.gsea.files)[m]),
                                                      parent = resultsFolder))
                
                save_to_synapse(temp.gsea.files[[m]], gseaFolder)
                setwd(file.path(base.path, contrast.name, names(all.files[i]), names(omics.files)[j]))
              } 
            }
          }
        }
      }
      # make space to process next contrast
      panSEA1 <- NULL
      panSEA2 <- NULL
      panSEA3a <- NULL
      panSEA3b <- NULL 
    }
  }
}

## load BeatAML data formatted for DMEA
# input: file path where BeatAML data should be saved
# output: list of BeatAML meta, drug AUC, global, and phospho data frames
load_BeatAML_for_DMEA <- function(BeatAML.path) {
  BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                             "Ex10_metadata.txt" = "syn25807733",
                             "ptrc_ex10_crosstab_global_gene_corrected.txt" = "syn25714248",
                             "ptrc_ex10_crosstab_phospho_siteID_corrected(1).txt" = "syn25714921")
  
  ### download files if any not already downloaded
  if (!file.exists(BeatAML.path)) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  } else if (!any(FALSE %in% lapply(names(BeatAML_synapse_id), file.exists))) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  }
  
  ### load files
  drug.BeatAML <- read.csv(file.path(BeatAML.path, names(BeatAML_synapse_id)[1]))
  meta.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[2]), 
                             sep = "\t", header = TRUE)
  global.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[3]),
                               sep = "\t", header = TRUE)
  phospho.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[4]),
                                sep = "\t", header = TRUE)
  
  ### format BeatAML data for DMEA
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
  
  return(list(meta = meta.BeatAML, drug = drug.BeatAML,
              global = global.BeatAML, phospho = phospho.BeatAML))
}
