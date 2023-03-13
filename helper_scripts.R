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


site_splitter <- function(m){
  f_data <- fData(m) %>%
    mutate(og_site = rownames(.)) %>%
    mutate(site_end = sub("^.*-([A-Za-z0-9]+)$", "\\1", og_site),
           site_base = sub("-[A-Za-z0-9]+$", "", og_site))
  
  f_data_expand <- data.frame()
  
  for (og_site_i in f_data$og_site){
    endings <- str_split(f_data[og_site_i, "site_end"], "[a-z]")[[1]] %>% head(-1)
    f_data_expand <- f_data_expand %>% 
      rbind(data.frame(og_site = og_site_i, split_site = paste(f_data[og_site_i, "site_base"], endings, sep = "-")))
  }
  
  m_long <- exprs(m) %>%
    as.data.frame() %>% mutate(og_site = rownames(.)) %>%
    tidyr::pivot_longer(-og_site, names_to = "Sample", values_to = "value") %>%
    filter(!is.na(value)) %>%
    left_join(f_data_expand, by = "og_site") %>%
    select(split_site, Sample, value)
  
  new_mat <- tidyr::pivot_wider(m_long, values_from = "value", names_from = "Sample", values_fn = mean) %>%
    as.data.frame()
  rownames(new_mat) <- new_mat$split_site
  new_mat <- new_mat[, -1]
  new_mat <- new_mat[, colnames(exprs(m))]
  
  f_data_split <- f_data_expand %>%
    group_by(split_site) %>% summarize(og_sites = list(unlist(og_site))) %>% as.data.frame()
  rownames(f_data_split) <- f_data_split$split_site
  
  m_split <- MSnSet(exprs = as.matrix(new_mat), pData = pData(m), fData = f_data_split[rownames(new_mat), ])
  return(m_split)
}




## Modification of MSnSetUtils function.
plot_pca <- function(eset, phenotype = NULL, shape = NULL, label = NULL, z_score = TRUE,
                     princomp_center = TRUE, show_ellipse = TRUE, components = 1:2, biplot = FALSE,
                     biplot_labels = NULL, standardize = TRUE, save_dfs = NULL,
                     num_features = 6L, show_NA = TRUE,
                     legend_title = phenotype,
                     arrow_args = list(), label_args = list(), ...) {
  
  # Handling coloring by phenotype. Do this first, in case
  # rows are removed when show_NA = FALSE
  if (!is.null(phenotype)) {
    colorBy <- pData(eset)[, phenotype]
    # If not showing missing values, remove those samples
    if (!show_NA) {
      idx <- !is.na(colorBy)
      eset <- eset[, idx]
      colorBy <- colorBy[idx]
    }
  } else {
    show_ellipse <- FALSE
    colorBy <- NULL
  }
  if (!is.null(shape)){
    shapeBy <- pData(eset)[, shape]
    if (!show_NA) {
      idx <- !is.na(shapeBy)
      eset <- eset[, idx]
      shapeBy <- shapeBy[idx]
    }
  } else {
    shapeBy <- NULL
  }
  
  # Check that components are valid
  if (length(components) != 2) {
    stop(sprintf("components must be a vector of length 2, not %d.",
                 length(components)))
  }
  if (!all(components %in% 1:ncol(eset))) {
    stop(sprintf("The values of components must be between 1 and %d.",
                 ncol(eset)))
  }
  
  complete_rows <- complete.cases(exprs(eset))
  
  # Check that there are enough complete rows for PCA
  if (sum(complete_rows) < 2) {
    stop("There are fewer than 2 rows with non-missing data.")
  }
  
  message(sprintf("Subsetting to %d complete rows for PCA.",
                  sum(complete_rows)))
  
  # Subset to complete rows
  eset <- eset[complete_rows, ]
  
  # If z_score, convert to Z-Scores by sample (row when transposed)
  if (z_score) {
    z <- t(scale(exprs(eset), center = TRUE, scale = TRUE))
  } else {
    z <- t(exprs(eset))
  }
  
  ## PCA
  # By default, center = TRUE, scale. = FALSE
  pca_res <- prcomp(z, center = princomp_center)
  
  u <- pca_res$x # Scores
  v <- pca_res$rotation # Eigenvectors
  
  if (standardize) {
    n <- nrow(u)
    lam <- pca_res$sdev * sqrt(n)
    
    # Scale u down and v up. Product is still the same
    u <- t(t(u) / lam)
    v <- t(t(v) * lam)
  }
  
  # Determine ratio between scale of v and u
  u_range <- apply(u[, components], 2, function(x) abs(range(x)))
  v_range <- apply(v[, components], 2, function(x) abs(range(x)))
  
  ratio <- max(v_range / u_range) # ratio for scaling v and secondary axes
  v <- v / ratio # scale v
  
  if (!is.null(save_dfs)){
    assign(save_dfs, list("sample_decomposition" = u %>% as.data.frame(), "feature_decomposition" = v %>% as.data.frame() %>%
                            mutate(feature = rownames(.))), envir = globalenv())
  }
  
  # Data frames for plotting
  df.u <- as.data.frame(u[, components])
  df.v <- as.data.frame(v[, components])
  
  # Percent of variance explained by each PC
  d <- pca_res$sdev # Standard deviations
  var_expl <- round(100 * d ^ 2 / sum(d ^ 2), digits = 2)[components]
  axis_labs <- sprintf("PC%d (%g%%)", #"%sPC%d (%g%%)",
                       # ifelse(obs.scale == 0, "Standardized ", ""),
                       components,
                       var_expl)
  
  # If colorBy is not NULL, add that column to df
  if (!is.null(colorBy)) {
    df.u$colorBy <- colorBy
  }
  if (!is.null(shapeBy)) {
    df.u$shapeBy <- shapeBy
  }
  
  ## Visualization
  # Base plot
  p <- ggplot(data = df.u, mapping = aes(x = df.u[, 1], y = df.u[, 2], color = colorBy, shape = shapeBy)) +
    geom_hline(yintercept = 0, lty = "longdash", color = "darkgrey") +
    geom_vline(xintercept = 0, lty = "longdash", color = "darkgrey") +
    labs(x = axis_labs[1], y = axis_labs[2]) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  # 50% confidence ellipse layer first so they are
  # beneath the layer of points or labels.
  if (show_ellipse & !is.numeric(colorBy)) {
    p <- p +
      stat_ellipse(mapping = aes(fill = colorBy, color = NULL),
                   geom = "polygon", type = "norm",
                   level = 0.5, alpha = 0.1, show.legend = TRUE)
  }
  
  # If label is NULL, add points. Otherwise, add labels
  if (is.null(label)) {
    p <- p +
      geom_point(...)
  } else {
    labels <- pData(eset)[, label]
    p <- p +
      geom_text(mapping = aes(label = labels), ...)
  }
  
  # Set titles for color and fill legend
  p <- p +
    guides(color = guide_legend(title = legend_title),
           fill = guide_legend(title = legend_title))
  
  # If colorBy is numeric, use a colorbar
  if (is.numeric(colorBy)) {
    p <- p +
      guides(color = guide_colorbar(title = legend_title))
  }
  
  ## Biplot
  if (biplot) {
    # Get the indices of the top influential features
    # from each principal component. num_features determines how
    # many to select from each component.
    top_features <- lapply(1:2, function(i) {
      order(abs(df.v)[, i], decreasing = TRUE)[1:num_features]
    })
    top_features <- unique(unlist(top_features))
    
    # Subset loadings to top features and rename columns
    df.v <- df.v[top_features, ]
    colnames(df.v) <- c("xend", "yend")
    df.v$x <- df.v$y <- 0
    
    # If biplot_labels is not provided, default to row names
    if (is.null(biplot_labels)) {
      df.v$labels <- rownames(df.v)
    } else {
      df.v$labels <- fData(eset)[top_features, biplot_labels]
    }
    
    scale_args <- list(expand = expansion(mult = rep(0.1, 2)),
                       sec.axis = sec_axis(~ . * ratio))
    
    # Arguments for geom_segment
    arrow_args <- list(mapping = aes(x = x, y = y, xend = xend, yend = yend),
                       arrow = arrow(length = unit(0.5, "line")),
                       data = df.v, color = "red3") %>%
      # Allow user-supplied args to overwrite defaults
      modifyList(val = arrow_args, keep.null = TRUE)
    
    # Arguments for geom_label_repel
    label_args <- list(mapping = aes(x = xend, y = yend, label = labels),
                       data = df.v,
                       color = arrow_args[["color"]],
                       max.overlaps = Inf,
                       min.segment.length = 0,
                       fill = alpha("white", 0.5)) %>%
      # Allow user-supplied args to overwrite defaults
      modifyList(val = label_args, keep.null = TRUE)
    
    # Add segments with arrows and text labels
    p <- p +
      # Add extra padding around plot area and secondary axes for v units
      do.call(scale_x_continuous, scale_args) +
      do.call(scale_y_continuous, scale_args) +
      do.call(geom_segment, arrow_args) +
      do.call(geom_label_repel, label_args) +
      theme(axis.text.y.right = element_text(color = arrow_args[["color"]]),
            axis.text.x.top = element_text(color = arrow_args[["color"]]),
            axis.ticks.y.right = element_line(color = arrow_args[["color"]]),
            axis.ticks.x.top = element_line(color = arrow_args[["color"]]))
  }
  
  return(p)
}


arrow_plotter <- function(arrow_df, top_n = 10, text_size = 4){
  
  # ppt_pathways <- c("REACTOME NEUTROPHIL DEGRANULATION", 
  #                   "REACTOME EXTRACELLULAR MATRIX ORGANIZATION", 
  #                   "REACTOME HEMOSTASIS", "REACTOME AMYLOID FIBER FORMATION", 
  #                   "REACTOME IRE1ALPHA ACTIVATES CHAPERONES", 
  #                   "REACTOME MITOCHONDRIAL FATTY ACID BETA OXIDATION", 
  #                   "REACTOME COMPLEX I BIOGENESIS", "REACTOME TRANSLATION", 
  #                   "REACTOME ANTIGEN PRESENTATION FOLDING ASSEMBLY AND PEPTIDE LOADING OF CLASS I MHC")
  # 
  # arrow_df <- data.frame(x_c = -4:4, y_c = sample(-4:4))
  # arrow_df <- arrow_df %>%
  #   mutate(label_name = ppt_pathways,
  #          angle = case_when(x_c < 0 ~ atan(y_c/x_c),
  #                            TRUE ~ atan(y_c/x_c)),
  #          angle = angle * 180/pi) 
  # top_n = 10
  # label_offset <- 0.007
  
  arrow_df$arrow_len <- sqrt(arrow_df$x_c**2 + arrow_df$y_c**2)
  # arrow_df$label_x = arrow_df$x_c - (label_offset*arrow_df$arrow_len*arrow_df$y_c)
  # arrow_df$label_y = arrow_df$y_c + (label_offset*arrow_df$arrow_len*arrow_df$x_c)
  
  arrow_df <- arrow_df %>% arrange(-arrow_len) %>% head(top_n)
  
  range <- max(abs(arrow_df$x_c), abs(arrow_df$y_c))
  
  
  p <- ggplot(arrow_df, aes(xend = x_c, yend = y_c)) + 
    geom_segment(x = 0, y = 0, arrow = arrow(angle = 25, length = unit(0.25, "cm"))) + 
    ggrepel::geom_text_repel(arrow_df, mapping = aes(x = x_c, y = y_c, label = label_name), 
                    nudge_x = 0.005, nudge_y = 0.005, box.padding = 0.5, 
                    segment.color = "black", segment.alpha = 1, size = text_size) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + xlim(-range, range) + ylim(-range, range)
  
  return(p)
}





