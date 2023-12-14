library(ggplot2)
library(ggrepel)
library(dplyr)

ppt_pathways <- c("REACTOME NEUTROPHIL DEGRANULATION", 
                  "REACTOME EXTRACELLULAR MATRIX ORGANIZATION", 
                  "REACTOME HEMOSTASIS", "REACTOME AMYLOID FIBER FORMATION", 
                  "REACTOME IRE1ALPHA ACTIVATES CHAPERONES", 
                  "REACTOME MITOCHONDRIAL FATTY ACID BETA OXIDATION", 
                  "REACTOME COMPLEX I BIOGENESIS", "REACTOME TRANSLATION", 
                  "REACTOME ANTIGEN PRESENTATION FOLDING ASSEMBLY AND PEPTIDE LOADING OF CLASS I MHC")

arrow_df <- data.frame(x_c = -4:4, y_c = sample(-4:4))
arrow_df <- arrow_df %>%
  mutate(label_name = ppt_pathways,
         angle = case_when(x_c < 0 ~ atan(y_c/x_c),
                           TRUE ~ atan(y_c/x_c)),
         angle = angle * 180/pi) 
top_n = 10
label_offset <- 0.007

arrow_df$arrow_len <- sqrt(arrow_df$x_c**2 + arrow_df$y_c**2)
arrow_df$label_x = arrow_df$x_c - (label_offset*arrow_df$arrow_len*arrow_df$y_c)
arrow_df$label_y = arrow_df$y_c + (label_offset*arrow_df$arrow_len*arrow_df$x_c)

arrow_df <- arrow_df %>% arrange(-arrow_len) %>% head(top_n)


p <- ggplot(arrow_df, aes(xend = x_c, yend = y_c)) + 
  geom_segment(x = 0, y = 0, arrow = arrow(angle = 25, length = unit(0.25, "cm"))) + 
  geom_text_repel(arrow_df, mapping = aes(x = label_x, y = label_y, label = label_name), 
                         nudge_x = 0.5, nudge_y = 0.5, box.padding = 0.5, 
                  segment.color = "white", segment.alpha = 1, size = 4) +
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
        plot.background=element_blank())
p








