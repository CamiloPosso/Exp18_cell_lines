library(shiny)
library(dplyr)
library(venn)
library(stringr)
library(pheatmap)
library(ggpubr)
library(ggpolypath)
library(readxl)

## plotting functions

plot_function <- function(enrichment_df, pathway, mat, cluster_cols, meta_df, ...){
  meta_df <- meta_df %>%
    mutate(State = case_when(State == "Spontaneous" ~ "Late",
                             TRUE ~ State))
  meta_df <- meta_df[colnames(mat), ]
  
  genes <- enrichment_df %>%
    filter(Pathway == pathway) %>%
    pull(Features) %>% unlist() %>% unique()
  # other_treatments <- setdiff(c("GVD", "GV", "GD", "G_Early", "G_Late"), treatments)
  # n_treat <- length(treatments)
  
  mat <- mat[intersect(genes, rownames(mat)), ]
  ann_col <- meta_df %>%
    select(Treatment, State)
  cat_colors1 <- c("darkorange2", "firebrick2", "mediumpurple", "dodgerblue3", "darkolivegreen4")
  cat_colors2 <- c("gray17", "gray42", "gray69")
  names(cat_colors1) <- unique(ann_col$Treatment)
  names(cat_colors2) <- unique(ann_col$State)
  # cat_colors <- c("ASF1A" = "darkorange2", "ASF1B" = "firebrick2",
  #                 "ASF1AB" = "mediumpurple", "ASF1B_IL" = "dodgerblue3")
  ann_colors <- list(`Treatment` = cat_colors1,
                     `State` = cat_colors2)
  
  color_pal <- rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
  color_pal <- c(colorRampPalette(color_pal[1:4])(15), colorRampPalette(color_pal[5:9])(20))
  y_max = (abs(min(mat, na.rm = T)) + abs(max(mat, na.rm = T)))/1.4
  breaks_buddy <- c(seq(-y_max, -0.001, length.out = 15),
                    0,
                    seq(0.001, y_max, length.out = 20))
  
  p <- pheatmap(mat, cluster_cols = cluster_cols, annotation_col = ann_col, annotation_colors = ann_colors, 
                main = pathway, width = 9, height = 12, fontsize_number = 13, breaks = breaks_buddy, 
                color = color_pal, fontsize_row = 9, ...)
}


#### App #####

# Define UI 
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("ASF TLK Pathway Heatmaps"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    
    h4("Select dataset and database"),
    selectInput("dataset", label = "Dataset:",
                choices = c("Exp 18" = "exp18",
                            "Exp 18 + Exp 12 MOLM14 samples" = "exp18_exp12")),
    selectInput("mode", label = "Data type:",
                choices = c("Global" = "global",
                            "Phospho" = "phospho")),
    uiOutput("global_db"),
    uiOutput("phospho_db"),
    
    h5("Click and delete the pathway below to search for a new one."),
    uiOutput("search_bar"),
    numericInput("width", "Width: ", 10),
    numericInput("height", "Height: ", 9),
    downloadButton('save_plot', 'Save plot')
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(width = 7,
            h3(textOutput("testing")),
            plotOutput("plot1", height = "600", width = "600"),
            # plotOutput("plot2", height = "600", width = "600")
  )
)


# Define server logic to plot stuff
server <- function(input, output, session) {
  
  meta_exp18 <- read_excel("pathway_heatmap_data/Exp18_metadata_01-14-2021.xlsx") %>%
    filter(!is.na(Index)) %>%
    mutate(Index = as.character(Index),
           Plex = as.character(Plex),
           comp1_ = paste(State, Treatment, sep = "_and_")) %>% as.data.frame() %>%
    mutate(Sample = paste0("X", Index)) %>%
    select(Sample, Treatment, State)
  rownames(meta_exp18) <- meta_exp18$Sample
  meta_combined <- read.table("pathway_heatmap_data/exp_12_18_combined_metadata.txt", sep = "\t") %>%
    select(Sample, Treatment, State)

  xx <- read.table("functional_data/Exp18_enrichment_combined.txt", sep = "\t")
  yy <- read.table("functional_data/Exp18_Exp12_enrichment_combined.txt", sep = "\t")
  
  global_crosstab_combined <- read.table("pathway_heatmap_data/exp_12_18_combined_global_corrected.txt", sep = "\t")
  phospho_crosstab_combined <- read.table("pathway_heatmap_data/exp_12_18_combined_phospho_corrected.txt", sep = "\t")
  global_crosstab_single <- read.table("pathway_heatmap_data/ptrc_ex18_crosstab_global_gene_corrected.txt")
  phospho_crosstab_single <- read.table("pathway_heatmap_data/ptrc_ex18_crosstab_phospho_siteID_corrected.txt")
  
  pathway_df <- readRDS("pathway_heatmap_data/pathway_features.RDS") %>%
    mutate(Features = str_split(Features, ", "))
  
  enrichment_df <- reactive({
    if (input$dataset == "exp18"){
      if (input$mode == "global"){
        xx %>% filter(DB == input$DB_g)
      } else {
        xx %>% filter(DB == input$DB_p)
      }
    } else {
      if (input$mode == "global"){
        yy %>%
          filter(DB == input$DB_g)
      } else {
        yy %>%
          filter(DB == input$DB_p)
      }
    }
  })
  
  meta_df <- reactive({
    if (input$dataset == "exp18"){
      meta_exp18
    } else {
      meta_combined
    }
  })
  
  mat <- reactive({
    if (input$dataset == "exp18"){
      if (input$mode == "global"){
        global_crosstab_single
      } else {
        phospho_crosstab_single
      }
    } else {
      if (input$mode == "global"){
        global_crosstab_combined
      } else {
        phospho_crosstab_combined
      }
    }
  })
  
  output$phospho_db <- renderUI({
    req(input$mode != "global")
    selectInput("DB_p", label = "Database:",
                choices = c("PSP + Networkin" = "KSEA NetworKIN 5"))
  })
  
  output$global_db <- renderUI({
    req(input$mode == "global")
    selectInput("DB_g", label = "Database:",
                choices = c("Reactome" = "GSEA Reactome",
                            "Hallmark" = "GSEA Hallmark",
                            "KEGG" = "GSEA KEGG",
                            "Gene Ontology BP" = "GSEA GOBP"))
  })
  
  output$search_bar <- renderUI({
    fluidRow(
      selectizeInput(
        inputId = "search_bar", 
        label = "Search Bar",
        multiple = FALSE,
        selected = "REACTOME ACTIVATION OF ATR IN RESPONSE TO REPLICATION STRESS",
        choices = c("Search Bar" = "", enrichment_df() %>% pull(Pathway) %>% unique()),
        options = list(
          create = FALSE,
          placeholder = "Search Me",
          maxItems = '1',
          onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
          onType = I("function (str) {if (str === \"\") {this.close();}}")
        )
      )
    )
  })
  
  output$plot1 <- renderPlot({
    plot_function(pathway_df, pathway(), mat(), TRUE, meta_df())
  })
  
  output$testing <- renderText({
    input$search_bar
  })
  pathway <- reactive({
    input$search_bar
  })
  cluster_cols <- reactive({
    input$cluster_cols
  })
  width <- reactive({
    input$width
  })
  height <- reactive({
    input$height
  })
  
  output$save_plot <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      file = paste(input$dataset, input$mode, paste0(pathway(), ".pdf"), sep = "_")
      print(file)
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      p <- reactive({
        plot_function(pathway_df, pathway(), mat(), TRUE, meta_df())
      })
      pdf(file = file, width = width(), height = height())
      p()
      dev.off()
    }
  )
  
}

shinyApp(ui, server)













