library(MSnSet.utils)
library(clusterProfiler)
library(dplyr)
library(KSEAapp)


GSEA_helper_single <- function(m, contrasts, t2g, t2g_name, coef.str){
  limma_res <- limma_contrasts(m, model.str = paste("~0 +", coef.str), 
                               coef.str = coef.str, contrasts = contrasts)
  
  plot_title <- gsub(coef.str, "", contrasts) %>% gsub("_", " ", .) %>% gsub("-", " vs ", .)
  group1 <- sub(" vs .*$", "", plot_title)
  group2 <- sub("^.* vs ", "", plot_title)
  plot_path <- gsub(" ", "_", plot_title) %>% paste0("GSEA_", ., paste0("_", t2g_name, ".png"))
  tbl_path <- sub(".png", ".txt", plot_path)
  
  ## For some pathways, gsea is unable to assess a pvalue. These are removed by the function
  ## Internally. This is why some pathways (like "GOBP NCRNA METABOLISM") are excluded from
  ## the result.
  
  if (file.exists(tbl_path)){
    print("using saved results")
    gsea_res <- read.table(tbl_path, sep = "\t")
  } else {
    fold_change <- limma_res$logFC
    names(fold_change) <- limma_res$feature
    fold_change <- sort(fold_change, decreasing = TRUE)
    set.seed(69)
    
    gsea_res <- GSEA(fold_change, eps = 1e-16, minGSSize = 10, 
                     pvalueCutoff = 1, TERM2GENE = t2g)@result %>%
      dplyr::select(Description, setSize, NES, pvalue, p.adjust, core_enrichment) %>%
      dplyr::mutate(contrast = contrasts)
    write.table(x = gsea_res, file = tbl_path, 
                sep = "\t", quote = F)
  }
  
  # yy <- gsea_res # %>%
  # #   compress_enrichment(colname = "p.adjust", threshold = 0.75) %>%
  # #   filter(p.adjust < 0.05)
  # plot_df <- yy %>% select(`Biological Process` = Description, NES) %>%
  #   mutate(Group = sign(NES)) %>%
  #   group_by(Group) %>%
  #   top_n(n = 10, wt = abs(NES)) %>%
  #   arrange(NES) %>%
  #   mutate(`Biological Process` = str_wrap(`Biological Process`, width = 50)) %>%
  #   mutate(`Biological Process` = factor(`Biological Process`, levels = `Biological Process`),
  #          Group = as.factor(Group))
  # 
  # p <- ggplot(plot_df, aes(x = NES, y = `Biological Process`, fill = Group)) + 
  #   geom_bar(stat = 'identity', color = "gray30") + 
  #   scale_fill_manual(values = c("#619CFF", "#F8766D"), 
  #                     labels = c(group2, group1)) + 
  #   ggtitle(paste0("GSEA - ", plot_title)) + theme_bw() +
  #   theme(axis.text = element_text(size = 18),
  #         axis.title = element_text(size = 18),
  #         legend.text = element_text(size = 18),
  #         title = element_text(size = 19))
  # 
  # ggsave(filename = paste0("Figures/", plot_path), plot = p, height = 12, width = 16)
}




GSEA_helper <- function(m, contrasts, t2g, t2g_name, coef.str){
  for (contrast in contrasts){
    GSEA_helper_single(m, contrast, t2g, t2g_name, coef.str)
  }

  file_pattern <- paste0("^GSEA_.*_vs_.*_", t2g_name, ".txt$")
  combined <- data.frame()
  for (file_path in list.files(pattern = file_pattern)){
    xx <- read.table(file_path, sep = "\t")
    combined <- rbind(combined, xx)
  }
  write.table(x = combined, file = paste0("GSEA_", t2g_name, "_combined.txt"), 
              sep = "\t", quote = F)
}






KSEA_helper_single <- function(m, contrasts, coef.str){
  plot_title <- gsub(coef.str, "", contrasts) %>% gsub("_", " ", .) %>% gsub("-", " vs ", .)
  group1 <- sub(" vs .*$", "", plot_title)
  group2 <- sub("^.* vs ", "", plot_title)
  plot_path <- gsub(" ", "_", plot_title) %>% paste0("KSEA_", ., paste0("_", "NetworKIN_5", ".png"))
  tbl_path <- sub(".png", ".txt", plot_path)
  
  KSDB <- read.csv(system.file('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                               package='amlresistancenetworks'),stringsAsFactors = FALSE)
  
  limma_res <- limma_contrasts(m, model.str = paste("~0 +", coef.str), 
                               coef.str = coef.str, contrasts = contrasts)
  
  fold_change <- limma_res$logFC
  fold_change <- 2**fold_change
  PX <- data.frame(Protein = "NULL", Gene = limma_res$feature, Peptide = "NULL", 
                   Residue.Both = limma_res$feature, p = "NULL", FC = fold_change) %>%
    dplyr::mutate(Residue.Both = sub("^.*-", "", Residue.Both)) %>%
    dplyr::mutate(Residue.Both = gsub("[a-z]", ";", Residue.Both)) %>%
    dplyr::mutate(Residue.Both = gsub(";$", "", Residue.Both),
                  Gene = sub("^(.*)-.*$", "\\1", Gene))
  
  ksea_res <- KSEA.Scores(KSDB, PX, NetworKIN = TRUE, NetworKIN.cutoff = 5) %>%
    dplyr::select(Kinase.Gene, m, FDR, z.score) %>%
    dplyr::rename(pathway = Kinase.Gene, enrichment = z.score,
                  adj_p_val = FDR, set_size = m) %>%
    dplyr::mutate(contrast = contrasts) %>%
    filter(set_size >= 3)
  write.table(x = ksea_res, file = tbl_path, 
              sep = "\t", quote = F)
  
}


KSEA_helper <- function(m, contrasts, coef.str){
  for (contrast in contrasts){
    print(contrast)
    KSEA_helper_single(m, contrast, coef.str)
  }
  
  file_pattern <- paste0("^KSEA_.*_vs_.*_", "NetworKIN_5", ".txt$")
  combined <- data.frame()
  for (file_path in list.files(pattern = file_pattern)){
    xx <- read.table(file_path, sep = "\t")
    combined <- rbind(combined, xx)
  }
  write.table(x = combined, file = paste0("KSEA_", "NetworKIN_5", "_combined.txt"), 
              sep = "\t", quote = F)
}





