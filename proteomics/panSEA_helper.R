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
  temp.result <- base.result$result
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
      mtn.file.names <- c(mtn.files, 
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
      ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                      device = "pdf")
    }
    
    # save to synapse
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
run_contrasts_global_phospho <- function(contrasts, contrast.type, meta.df,
                                         omics, beatAML, gmt.features, 
                                         gmt.features2, gmt.features3a,
                                         gmt.features3b,
                                         gmt.drug, drug.BeatAML, 
                                         base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                         synapse_id) {
  for (k in 1:length(contrasts)) {
    # identify samples for each side of contrast
    c1 <- contrasts[[k]][1]
    c2 <- contrasts[[k]][2]
    group.names <- c(c1, c2)
    contrast.name <- paste0(contrasts[[k]][1], "_vs_", contrasts[[k]][2])
    group.samples <- list(meta.df[meta.df[,contrast.type] == c1, ]$MeasurementName,
                          meta.df[meta.df[,contrast.type] == c2, ]$MeasurementName)
    
    # run panSEA across omics types
    # CAUTION: this only works because the # of samples for each treatment type 
    # is equal; otherwise would have to run panSEA for each contrast separately 
    # and perhaps set the group.samples input parameter for panSEA
    
    # run global against KEGG sets, phospho against kinase sets from ksdb
    panSEA1 <- panSEA::panSEA(omics, c("global", "phospho"), 
                              feature.names = c("Gene", "SUB_SITE"), 
                              group.names = group.names,
                              group.samples = group.samples,
                              gmt.features = gmt.features,
                              gmt.drugs = gmt.drug,
                              drug.sensitivity = drug.BeatAML,
                              expression = beatAML, 
                              min.per.set = 6)
    
    # use DEG results to avoid re-running DEG analysis
    # run global against hallmark sets, phospho against substrate sets
    deg.results <- list("global" = panSEA1$mDEG.results$all.results$global,
                        "phospho" = panSEA1$mDEG.results$all.results$phospho)
    panSEA2 <- panSEA::panSEA(deg.results, names(deg.results), DMEA = FALSE,
                              group.names = contrast.name,
                              group.samples = 2, # "Log2FC" is in column 2
                              feature.names = c("Gene", "SUB_SITE"), 
                              gmt.features = gmt.features2,
                              min.per.set = 6)
    
    ksea.results <- panSEA1$mGSEA.results$phospho$result
    sub.results <- panSEA2$mGSEA.results$phospho$result
    
    # check if 2+ gene sets are covered by KSEA or substrate enrichment results
    # get KEGG info
    msigdb.KEGG <- msigdbr::msigdbr("Homo sapiens", "C2", "CP")
    msigdb.KEGG <- as.data.frame(msigdb.info[, c(
      "gene_symbol",
      "gs_name",
      "gs_description"
    )])
    colnames(msigdb.KEGG) <- c("Feature_set", "gs_name", "gs_description")
    
    msigdb.H <- msigdbr::msigdbr("Homo sapiens", "H")
    msigdb.H <- as.data.frame(msigdb.info[, c(
      "gene_symbol",
      "gs_name",
      "gs_description"
    )])
    colnames(msigdb.H) <- c("Feature_set", "gs_name", "gs_description")
    
    # merge with kinase, substrate results
    kin.KEGG <- merge(ksea.results, msigdb.KEGG)
    kin.H <- merge(ksea.results, msigdb.H)
    
    sub.KEGG <- merge(sub.results, msigdb.KEGG)
    sub.H <- merge(sub.results, msigdb.H)
    
    # create gmts to determine if we can do further GSEA
    gmt.kin.KEGG <- DMEA::as_gmt(
      kin.KEGG, "Feature_set", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    ) 
    gmt.kin.H <- DMEA::as_gmt(
      kin.H, "Feature_set", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    ) 
    
    gmt.sub.KEGG <- DMEA::as_gmt(
      sub.KEGG, "Feature_set", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    ) 
    gmt.sub.H <- DMEA::as_gmt(
      sub.H, "Feature_set", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    ) 
    
    if (length(gmt.kin.KEGG) > 1 & length(gmt.sub.KEGG) > 1) {
      phospho.a <- list("KSEA" = ksea.results,
                        "Substrate enrichment" = sub.results)
    } else if (length(gmt.sub.KEGG) > 1) {
      phospho.a <- list("Substrate enrichment" = sub.results)
    } else if (length(gmt.kin.KEGG) > 1) {
      phospho.a <- list("KSEA" = ksea.reuslts)
    } else {
      phospho.a <- NULL
    }
    
    if (length(gmt.kin.H) > 1 & length(gmt.sub.H) > 1) {
      phospho.b <- list("KSEA" = ksea.results,
                        "Substrate enrichment" = sub.results)
    } else if (length(gmt.sub.H) > 1) {
      phospho.b <- list("Substrate enrichment" = sub.results)
    } else if (length(gmt.kin.H) > 1) {
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
    # get DMEA mtn plots
    DMEA.global.mtn <- get_top_mtn_plots(panSEA1$mDMEA.results$all.results$global,
                                         EA.type = "DMEA")
    DMEA.phospho.mtn <- get_top_mtn_plots(panSEA1$mDMEA.results$all.results$phospho,
                                          EA.type = "DMEA")
    
    # get GSEA mtn plots
    global.mtn1 <- get_top_mtn_plots(panSEA1$mGSEA.results$all.results$global,
                                     EA.type = "GSEA_KEGG")
    global.mtn2 <- get_top_mtn_plots(panSEA2$mGSEA.results$all.results$global,
                                     EA.type = "GSEA_hallmark")
    phospho.mtn1 <- get_top_mtn_plots(panSEA1$mGSEA.results$all.results$phospho,
                                      EA.type = "KSEA")
    phospho.mtn2 <- get_top_mtn_plots(panSEA2$mGSEA.results$all.results$phospho,
                                      EA.type = "Substrate_enrichment")
    
    if (is.null(panSEA3a)) {
      phospho.KSEA.KEGG.files <- list()
      phospho.KSEA.hallmark.files <- list()
      
      phospho.sub.KEGG.files <- list()
      phospho.sub.hallmark.files <- list()
    } else if (is.na(panSEA1$mGSEA.network$phospho$static)) {
      phospho.KSEA.KEGG.files <- list()
      phospho.KSEA.hallmark.files <- list()
      
      sub.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$all.results$'Substrate enrichment',
                                     EA.type = "Substrate_enrichment_GSEA_KEGG")
      sub.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$all.results$'Substrate enrichment',
                                     EA.type = "Substrate_enrichment_GSEA_hallmark")
      phospho.sub.KEGG.files <- list("Substrate_enrichment_GSEA_KEGG_results.csv" =
                                       panSEA3a$mGSEA.results$all.results$phospho$result,
                                     "Substrate_enrichment_GSEA_KEGG_volcano_plot.pdf" =
                                       panSEA3a$mGSEA.results$all.results$phospho$volcano.plot,
                                     "mtn_plots" = sub.mtn3a)
      phospho.sub.hallmark.files <- list("Substrate_enrichment_GSEA_hallmark_results.csv" =
                                           panSEA3b$mGSEA.results$all.results$phospho$result,
                                         "Substrate_enrichment_GSEA_hallmark_volcano_plot.pdf" =
                                           panSEA3b$mGSEA.results$all.results$phospho$volcano.plot,
                                         "mtn_plots" = sub.mtn3b)
    } else if (is.na(panSEA2$mGSEA.network$phospho$static)) {
      kin.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$all.results$'KSEA', 
                                     EA.type = "KSEA_GSEA_KEGG")
      kin.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$all.results$'KSEA',
                                     EA.type = "KSEA_GSEA_hallmark")
      phospho.KSEA.KEGG.files <- list("KSEA_GSEA_KEGG_results.csv" =
                                        panSEA3a$mGSEA.results$all.results$phospho$result,
                                      "KSEA_GSEA_KEGG_volcano_plot.pdf" =
                                        panSEA3a$mGSEA.results$all.results$phospho$volcano.plot,
                                      "mtn_plots" = kin.mtn3a)
      phospho.KSEA.hallmark.files <- list("KSEA_GSEA_hallmark_results.csv" =
                                            panSEA3b$mGSEA.results$all.results$phospho$result,
                                          "KSEA_GSEA_hallmark_volcano_plot.pdf" =
                                            panSEA3b$mGSEA.results$all.results$phospho$volcano.plot,
                                          "mtn_plots" = kin.mtn3b)
      
      phospho.sub.KEGG.files <- list()
      phospho.sub.hallmark.files <- list()
    } else {
      kin.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$all.results$'KSEA', 
                                     EA.type = "KSEA_GSEA_KEGG")
      kin.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$all.results$'KSEA',
                                     EA.type = "KSEA_GSEA_hallmark")
      phospho.KSEA.KEGG.files <- list("KSEA_GSEA_KEGG_results.csv" =
                                        panSEA3a$mGSEA.results$all.results$phospho$result,
                                      "KSEA_GSEA_KEGG_volcano_plot.pdf" =
                                        panSEA3a$mGSEA.results$all.results$phospho$volcano.plot,
                                      "mtn_plots" = kin.mtn3a)
      phospho.KSEA.hallmark.files <- list("KSEA_GSEA_hallmark_results.csv" =
                                            panSEA3b$mGSEA.results$all.results$phospho$result,
                                          "KSEA_GSEA_hallmark_volcano_plot.pdf" =
                                            panSEA3b$mGSEA.results$all.results$phospho$volcano.plot,
                                          "mtn_plots" = kin.mtn3b)
      
      sub.mtn3a <- get_top_mtn_plots(panSEA3a$mGSEA.results$all.results$'Substrate enrichment',
                                     EA.type = "Substrate_enrichment_GSEA_KEGG")
      sub.mtn3b <- get_top_mtn_plots(panSEA3b$mGSEA.results$all.results$'Substrate enrichment',
                                     EA.type = "Substrate_enrichment_GSEA_hallmark")
      phospho.sub.KEGG.files <- list("Substrate_enrichment_GSEA_KEGG_results.csv" =
                                       panSEA3a$mGSEA.results$all.results$phospho$result,
                                     "Substrate_enrichment_GSEA_KEGG_volcano_plot.pdf" =
                                       panSEA3a$mGSEA.results$all.results$phospho$volcano.plot,
                                     "mtn_plots" = sub.mtn3a)
      phospho.sub.hallmark.files <- list("Substrate_enrichment_GSEA_hallmark_results.csv" =
                                           panSEA3b$mGSEA.results$all.results$phospho$result,
                                         "Substrate_enrichment_GSEA_hallmark_volcano_plot.pdf" =
                                           panSEA3b$mGSEA.results$all.results$phospho$volcano.plot,
                                         "mtn_plots" = sub.mtn3b)
    }
    
    # set file names
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
                               panSEA1$mDMEA.results$compiled.results$dot.plot,
                             "DMEA_interactive_network_graph.html" =
                               panSEA1$mDMEA.network$interactive,
                             "DMEA_static_network_graph.html" =
                               panSEA1$mDMEA.network$static)
    global.DMEA.files <- list("DMEA_results.csv" =
                                panSEA1$mDMEA.results$all.results$global$result,
                              "DMEA_volcano_plot.pdf" =
                                panSEA1$mDMEA.results$all.results$global$volcano.plot,
                              "mtn_plots" = DMEA.global.mtn)
    phospho.DMEA.files <- list("DMEA_results.csv" =
                                 panSEA1$mDMEA.results$all.results$phospho$result,
                               "DMEA_volcano_plot.pdf" =
                                 panSEA1$mDMEA.results$all.results$phospho$volcano.plot,
                               "mtn_plots" = DMEA.phospho.mtn)
    
    global.GSEA.KEGG.files <- list("GSEA_KEGG_results.csv" =
                                     panSEA1$mGSEA.results$all.results$global$result,
                                   "GSEA_KEGG_volcano_plot.pdf" =
                                     panSEA1$mGSEA.results$all.results$global$volcano.plot,
                                   "mtn_plots" = global.mtn1)
    kin.files <- list('GSEA_KEGG' = phospho.kin.KEGG.files,
                      'GSEA_hallmark' = phospho.kin.hallmark.files)
    phospho.KSEA.files <- list("KSEA_results.csv" =
                                 panSEA1$mGSEA.results$all.results$phospho$result,
                               "KSEA_volcano_plot.pdf" =
                                 panSEA1$mGSEA.results$all.results$phospho$volcano.plot,
                               "mtn_plots" = phospho.mtn1,
                               "GSEA" = kin.files)
    
    global.GSEA.hallmark.files <- list("GSEA_hallmark_results.csv" =
                                         panSEA2$mGSEA.results$all.results$global$result,
                                       "GSEA_hallmark_volcano_plot.pdf" =
                                         panSEA2$mGSEA.results$all.results$global$volcano.plot,
                                       "mtn_plots" = global.mtn2)
    sub.files <- list('GSEA_KEGG' = phospho.sub.KEGG.files,
                      'GSEA_hallmark' = phospho.sub.hallmark.files)
    phospho.sub.files <- list("Substrate_enrichment_results.csv" =
                                panSEA2$mGSEA.results$all.results$phospho$result,
                              "Substrate_enrichment_volcano_plot.pdf" =
                                panSEA2$mGSEA.results$all.results$phospho$volcano.plot,
                              "mtn_plots" = phospho.mtn2,
                              "GSEA" = sub.files)
    
    global.files <- list('Differential_expression' = global.DEG.files, 
                         'DMEA' = global.DMEA.files,
                         'GSEA_KEGG' = global.GSEA.hallmark.files,
                         'GSEA_hallmark' = global.GSEA.KEGG.files)
    phospho.files <- list('Differential_expression' = phospho.DEG.files, 
                          'DMEA' = phospho.DMEA.files,
                          'KSEA' = phospho.KSEA.files,
                          'Substrate_enrichment' = phospho.sub.files)
    all.files <- list('global_and_phospho' = combo.DMEA.files,
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
        dir.create(names(omics.files)[i])
        setwd(names(omics.files)[i])
        temp.files <- all.files[[i]]
        resultsFolder <- 
          synapser::synStore(synapser::Folder(names(omics.files)[i],
                                              parent = omicsFolder))
        
        # save to synapse
        save_to_synapse(temp.files, resultsFolder)
        
        # save mtn plots if relevant
        mtn.files <- names(temp.files)[grepl("mtn_plots", names(temp.files))]
        if (length(mtn.files) > 0) {
          # create folder for mtn plots
          dir.create("mtn_plots")
          setwd("mtn_plots")
          mtnFolder <- 
            synapser::synStore(synapser::Folder("mountain_plots",
                                                parent = resultsFolder))
          
          temp.mtn.files <- temp.files[["mtn_plots"]]
          for (m in 1:length(temp.mtn.files)) {
            ggplot2::ggsave(names(temp.mtn.files)[m], temp.mtn.files[[m]], 
                            device = "pdf")
          }
          
          # upload mtn plots to synapse
          PDFs <- lapply(as.list(names(temp.mtn.files)), synapser::File,
                         parent = mtnFolder)
          lapply(PDFs, synapser::synStore)
        }
        setwd(file.path(base.path, names(all.files[i]), names(omics.files)[j]))
        
        # save KSEA_GSEA or sub_GSEA if relevant
        gsea.files <- names(temp.files)[grepl("GSEA", names(temp.files))]
        if (length(gsea.files) > 0) {
          for (m in 1:length(gsea.files)) {
            # create folder for mtn plots
            dir.create(names(gsea.files)[m])
            setwd(names(gsea.files)[m])
            gseaFolder <- 
              synapser::synStore(synapser::Folder(names(gsea.files)[m],
                                                  parent = resultsFolder))
            
            save_to_synapse(gsea.files[[m]], gseaFolder)
            setwd(file.path(base.path, contrast.name, names(all.files[i]), names(omics.files)[j]))
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
  } else if (any(!lapply(names(BeatAML_synapse_id), file.exists))) {
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
