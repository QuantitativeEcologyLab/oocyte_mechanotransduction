#This script generates figures S4 and S5 using reactome data

# Written by Michael Noonan
# Last updated: March 25, 2026

#Load in any necessary packages
library(ggplot2)
library(scico)

#Function for plotting Reactome pathway enrichment results
plot_reactome <- function(csv_file, outfile, plot_title) {
  
  reactome_data <- read.csv(csv_file)
  
  #Select top 20 pathways by combined score and clean pathway names
  filtered_r <- reactome_data[order(reactome_data$Combined.Score, decreasing = TRUE), ]
  filtered_r <- head(filtered_r, 20)
  filtered_r$Term <- gsub(" R-HSA-[0-9]+", "", filtered_r$Term)
  
  p <- ggplot(filtered_r,
              aes(x = Combined.Score,
                  y = reorder(Term, Combined.Score),
                  size = Combined.Score,
                  color = -log10(Adjusted.P.value))) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(name = "Combined Score") +
    scale_color_gradientn(colors = scico(100, palette = "batlow"),
                          name = "-log10(Adj. P-value)") +
    labs(x = "Combined Score", y = NULL, title = plot_title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"),
          legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 12, face = "bold"),
          legend.position = "inside",
          legend.position.inside = c(0.85, 0.25),
          axis.line = element_line(linetype = "solid"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          axis.text.y = element_text(size = 8))
  
  ggsave(outfile, plot = p, width = 10, height = 8,
         dpi = 600, bg = "transparent")
}


#-------------------------------------------------------------
#Figure S4 - Pathway-level analysis of oocytes matured under dynamic, static, or 2D conditions
#-------------------------------------------------------------

#Figure S4a - Dynamic vs. 2D
plot_reactome("data/Selected_Reactome_down_2DvsDyn.csv",
              "figures/figure_s4a_reactome_pathway_down_dynamic_vs_2D.png",
              "Up-regulated pathways in dynamic vs. 2D cultures")


#Figure S4b - Dynamic vs. Static
plot_reactome("data/Reactome_2022_table_up_DynvsST.csv",
              "figures/figure_s4b_reactome_pathway_down_dynamic_vs_static.png",
              "Up-regulated pathways in dynamic vs. static cultures")


#-------------------------------------------------------------
#Figure S5 - Pathway-level analysis of blastocysts derived from oocytes matured under dynamic, or 2D conditions
#-------------------------------------------------------------

#Figure S5a - Dynamic cell culture
plot_reactome("data/Reactome_Pathways_2024_unique_Dyn.csv",
              "figures/figure_s5a_reactome_unique_Dynamic.png",
              "Reactome Pathway Enrichment: Unique Dynamic genes")
#Figure S5b - 2D cell culture
plot_reactome("data/Reactome_Pathways_2024_unique_2D.csv",
              "figures/figure_s5b_reactome_unique_2D.png",
              "Reactome Pathway Enrichment: Unique 2D genes")
