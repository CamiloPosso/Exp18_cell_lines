library(shiny)
library(dplyr)
library(venn)
library(stringr)
library(pheatmap)
library(ggpubr)
library(ggpolypath)

## plotting functions

plot_function <- function(enrichment_df, database, treatments, plot_type, p_cutoff, n_pathways, exclude_rest = T){
  plot_title <- paste(treatments, collapse = ", ")
  enrichment_df <- enrichment_df %>%
    filter(DB == database) %>%
    mutate(Pathway = gsub("_", " ", Pathway),
           Pathway = str_wrap(Pathway, width = 40))
  other_treatments <- setdiff(c("GVD", "GV", "GD", "G_Early", "G_Late"), treatments)
  n_treat <- length(treatments)
  
  enrichment_df_s <- enrichment_df %>%
    filter(p.adjust < p_cutoff,
           Treatment %in% treatments) %>%
    unique() %>% arrange(-p.adjust) %>% group_by(Pathway) %>%
    mutate(total = n(), positive = (sign(enrichment) == 1)) %>% 
    mutate(plot_group = case_when(total == n_treat & all(positive) ~ "sim",
                                  total == n_treat & all(!positive) ~ "sim",
                                  TRUE ~ "diff"))
  
  ## These are kinases which show matching significant activity in both GVD and
  ## the combo G with the excluded treatment.
  excluded_pathways_diff <- enrichment_df %>%
    filter(Pathway %in% enrichment_df_s$Pathway,
           Treatment %in% other_treatments,
           p.adjust < p_cutoff) %>%
    pull(Pathway) %>% unique()
  
  excluded_pathways_common <- enrichment_df %>%
    filter(Pathway %in% enrichment_df_s$Pathway,
           Treatment %in% other_treatments,
           p.adjust < p_cutoff) %>%
    group_by(Pathway) %>%
    mutate(positive = (sign(enrichment) == 1),
           plot_group = case_when(all(positive) ~ "sim",
                                  all(!positive) ~ "sim",
                                  TRUE ~ "diff")) %>%
    ungroup() %>%
    filter(plot_group == "sim")
  
  exclude_positive <- excluded_pathways_common %>% filter(positive) %>% 
    pull(Pathway) %>% unique() %>% 
    intersect(enrichment_df_s %>% filter(plot_group == "sim", positive) %>% pull(Pathway) %>% unique())
  exclude_negative <- excluded_pathways_common %>% filter(!positive) %>% 
    pull(Pathway) %>% unique() %>%
    intersect(enrichment_df_s %>% filter(plot_group == "sim", !positive) %>% pull(Pathway) %>% unique())
  
  excluded_pathways_common <- c(exclude_positive, exclude_negative)
  
  if (!exclude_rest){
    excluded_pathways_diff <- c()
    excluded_pathways_common <- c()
  } else {
    plot_title <- paste(plot_title, "(excluding others)")
  }
  
  chosen_terms_diff <- enrichment_df_s %>%
    filter(plot_group == "diff",
           !Pathway %in% excluded_pathways_diff) %>%
    group_by(Treatment) %>% top_n(n = n_pathways, wt = -p.adjust) %>%
    pull(Pathway) %>% unique()
  
  chosen_terms_sim <- enrichment_df_s %>%
    filter(plot_group == "sim",
           !Pathway %in% excluded_pathways_common) %>%
    group_by(Treatment) %>% top_n(n = n_pathways, wt = -p.adjust) %>%
    pull(Pathway) %>% unique()
  
  plots <- list()
  
  if (plot_type == "diff"){
    chosen_terms <- chosen_terms_diff
    # plot_path <- paste0("gsea_", t2g_name, "_", chosen_treatment, "_combinations_", chosen_state, "_differences.png")
    plot_title <- paste0(plot_title, " - differences")
  } else {
    chosen_terms <- chosen_terms_sim
    # plot_path <- paste0("gsea_", t2g_name, "_", chosen_treatment, "_combinations_", chosen_state, "_similarties.png")
    plot_title <- paste0(plot_title, " - similarities")
  }
  
  if (length(chosen_terms) >= 2){
    mat_treatments <- c("G_Early", "G_Late", "GD", "GV", "GVD")
    mat <- enrichment_df %>%
      filter(Pathway %in% chosen_terms) %>%
      select(Pathway, enrichment, Treatment) %>%
      tidyr::pivot_wider(values_from = enrichment, names_from = Treatment)
    mat.rownames <- mat$Pathway
    mat <- mat[, -1] %>% as.matrix()
    rownames(mat) <- mat.rownames
    
    mat <- mat[, mat_treatments]
    
    mat_pval <- enrichment_df %>%
      filter(Pathway %in% chosen_terms) %>%
      select(Pathway, p.adjust, Treatment) %>%
      tidyr::pivot_wider(values_from = p.adjust, names_from = Treatment)
    mat_pval <- mat_pval[, -1] %>% as.matrix()
    mat_pval <- mat_pval[, mat_treatments]
    index_pval <- mat_pval < p_cutoff
    mat_pval[mat_pval != -1] <- ""
    mat_pval[index_pval] <- "*"
    
    ann_col <- data.frame(Treatment = sub(" .+$", "", colnames(mat)))
    rownames(ann_col) <- colnames(mat)
    cat_colors <- c("G_Early" = "darkorange2", "G_Late" = "firebrick2",
                    "GD" = "mediumpurple", "GV" = "dodgerblue3", "GVD" = "darkseagreen2")
    for (treatment in  other_treatments){
      cat_colors[[treatment]] <- "white"
    }
    ann_colors <- list(`Treatment` = cat_colors)
    
    color_pal <- rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
    color_pal <- c(colorRampPalette(color_pal[1:4])(15), colorRampPalette(color_pal[5:9])(20))
    breaks_buddy <- c(seq(-2.5, -0.001, length.out = 15),
                      0,
                      seq(0.001, 2.5, length.out = 20))
    
    p <- pheatmap(mat, cluster_cols = F, display_numbers = mat_pval, annotation_col = ann_col,
                                   annotation_colors = ann_colors, main = plot_title, width = 9, height = 12,
                                   fontsize_number = 13, breaks = breaks_buddy, color = color_pal, fontsize_row = 9)
  } else {
    p <- ggplot() + ggtitle("Less than 2 pathways")
  }
  return(p)
}


venn_plot_function <- function(enrichment_df, database, p_cutoff){
  enrichment_df <- enrichment_df %>%
    filter(DB == database)
  
  enrichment_df_s <- enrichment_df %>%
    filter(p.adjust < p_cutoff) %>%
    select(Treatment, Pathway, enrichment)
  rownames(enrichment_df_s) <- 1:nrow(enrichment_df_s)
  
  up_paths <-lapply(unique(enrichment_df$Treatment), function(x) {
                                                       enrichment_df_s %>% 
                                                       filter(Treatment == x,
                                                              enrichment > 0) %>% 
                                                       pull(Pathway)})
  down_paths <-lapply(unique(enrichment_df$Treatment), function(x) {
                                                         enrichment_df_s %>% 
                                                         filter(Treatment == x,
                                                                enrichment < 0) %>% 
                                                         pull(Pathway)})
  names(up_paths) <- unique(enrichment_df$Treatment)
  names(down_paths) <- unique(enrichment_df$Treatment)
  p_up <- venn(up_paths, zcolor = 'style', ggplot = TRUE, plotsize = 45) + ggtitle("")
  p_down <- venn(down_paths, zcolor = 'style', ggplot = TRUE, plotsize = 45) + ggtitle("")
  # p <- ggarrange(p_up, p_down, ncol = 1)
  cowplot::plot_grid(p_up, p_down, ncol = 1,
                     labels = c("Overexpressed vs Parental", "Underexpressed vs Parental"), label_size = 10)
}


#### App #####

# Define UI 
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("Functional enrichment vs Parental"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    
    h6("We use GSEA and KSEA to compare each treatment to the Parental samples. Select the p-value cutoff and the desired treatments to compare."),
    h6("The Venn diagrams below count how many overexpressed and underexpressed pathways are common between the various treatments"),
    
    numericInput("p_cutoff", "Adjusted p-value cutoff:", 0.05, 0, 0.15),
    plotOutput("summaryplot", height = "450", width = "300"),
  
    h4("Select data and database"),
    selectInput("mode", label = "Mode:",
                choices = c("Global" = "global",
                            "Phospho" = "phospho")),
    uiOutput("global_db"),
    uiOutput("phospho_db"),
    
    h4("Select treatments"),
    
    checkboxGroupInput("treatments", label = "Treatment", 
                       choices = c("G_Early", "G_Late", "GD", "GV", "GVD"), selected = c("GD", "GV")),
    numericInput("n_pathways", "Number of Pathways per Treatment", 8, 2, 151),
    h6("Select here to remove pathways active in the unselected treatments from the heatmaps"),
    checkboxInput("exclude_rest", label = "Exclude the other treatments", TRUE),
    downloadButton('save_plot', 'Save plot'),
    selectInput("save_mode", label = "Which plot to save?",
                choices = c("Differences" = "diff",
                            "Similarities" = "sim"))
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(width = 7,
    h3(textOutput("data_type")),
    plotOutput("plot1", height = "600", width = "600"),
    plotOutput("plot2", height = "600", width = "600")
  )
)


# Define server logic to plot stuff
server <- function(input, output, session) {
  print(getwd())
  xx <- read.table("enrichment_combined.txt", sep = "\t")
  
  output$phospho_db <- renderUI({
    req(input$mode == "phospho")
    selectInput("DB_p", label = "Database:",
                choices = c("PSP + NetworKIN 5" = "KSEA NetworKIN 5"))
  })
  
  output$global_db <- renderUI({
    req(input$mode == "global")
    selectInput("DB_g", label = "Database:",
                choices = c("GO Biological Process" = "GSEA GOBP",
                            "KEGG" = "GSEA KEGG",
                            "Hallmark" = "GSEA Hallmark",
                            "Reactome" = "GSEA Reactome"))
  })
  p_cutoff <- reactive({input$p_cutoff})
  n_pathways <- reactive({input$n_pathways})
  data_type <- reactive({input$mode})
  appdata <- reactive({
    if (input$mode == "global"){
      xx %>%
        filter(DB == input$DB_g)
    } else {
      xx %>%
        filter(DB == input$DB_p)
    }
  })
  chosen_DB <- reactive({
    if (input$mode == "global"){
      input$DB_g
    } else {
      input$DB_p
    }
  })
  # output$data_type <- renderText({nrow(appdata())})
  # output$data_type <- renderText({chosen_DB()})
  output$data_type <- renderText({treatments()})
  treatments <- reactive({
    input$treatments
  })
  exclude_rest <- reactive({
    input$exclude_rest
  })
  
  output$plot1 <- renderPlot({
    plot_function(xx, chosen_DB(), treatments(), "sim", p_cutoff(), n_pathways(), exclude_rest())
  })
  output$plot2 <- renderPlot({
    plot_function(xx, chosen_DB(), treatments(), "diff", p_cutoff(), n_pathways(), exclude_rest())
  })
  
  output$summaryplot <- renderPlot({
    p = venn_plot_function(xx, chosen_DB(), p_cutoff())
    print(p)
  })
  
  save_mode <- reactive({
    input$save_mode
  })
  
  output$save_plot <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(input$mode, chosen_DB(), input$p_cutoff, paste0(save_mode(), ".pdf"), sep = "_")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      p <- reactive({
        plot_function(xx, chosen_DB(), treatments(), save_mode(), p_cutoff(), n_pathways(), exclude_rest())
      })
      
      pdf(file = file)
      p()
      dev.off()
    }
  )
}

shinyApp(ui, server)













