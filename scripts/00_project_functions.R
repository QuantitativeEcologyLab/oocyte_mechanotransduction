#This script builds some helper functions used throughout the project

# Written by Michael Noonan
# Last updated: March 25, 2026


#-------------------------------------------------------------
#Function for adding gene names
#-------------------------------------------------------------

add_genename <- function(sig_df, gene_expression_data) {
  df        <- as.data.frame(sig_df)
  df$gene_id <- rownames(df)
  df <- merge(df, gene_expression_data[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
  df$gene_name <- ifelse(is.na(df$gene_name), df$gene_id, df$gene_name)
  df
}

#-------------------------------------------------------------
#Function for generating the volcano plots for figure 2
#-------------------------------------------------------------

create_volcano_plot <- function(results_df, comparison_label, output_filename,
                                alpha_threshold = 0.05, fold_change_cutoff = 1.0, top_n = 30) {
  
  #Replace missing gene names with gene IDs
  results_df$gene_name <- ifelse(is.na(results_df$gene_name), results_df$gene_id, results_df$gene_name)
  
  #Classify significance
  results_df$Significance <- ifelse(!is.na(results_df$padj) & results_df$padj < alpha_threshold, 
                                    "Significant", "Not Significant")
  
  #Get top N significant genes by adjusted p-value
  sig_df   <- results_df[results_df$Significance == "Significant", ]
  sig_df   <- sig_df[order(sig_df$padj), ]
  top_genes <- head(sig_df$gene_name, top_n)
  
  p <- 
    ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Significance), size = 2, alpha = 0.6) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red", name = "")) +
    scale_x_continuous(limits = c(-30, 30), oob = scales::squish) +
    geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linetype = "solid"),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.position = "none",
          legend.position.inside = c(0.8,0.95),
          legend.title = element_blank(),
          legend.text = element_text(size=8, family = "sans", face = "bold"),
          legend.background = element_rect(fill = "transparent"),
          legend.key.size = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          plot.title = element_text(size = 12, family = "sans", face = "bold"),
          axis.title.y = element_text(size=12, family = "sans", face = "bold"),
          axis.title.x = element_text(size=12, family = "sans", face = "bold"),
          axis.text.y = element_text(size=10, family = "sans"),
          axis.text.x  = element_text(size=10, family = "sans"),
          strip.background=element_blank(),
          plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
    labs(title = paste(comparison_label),
         x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
    ggrepel::geom_text_repel(data = subset(results_df, gene_name %in% top_genes), aes(label = gene_name),
                             size = 3, box.padding = 0.5, point.padding = 0.5,
                             max.overlaps = Inf, segment.color = "grey50") 
  
  ggsave(p, width = 8, height = 6,
         units = "in",
         dpi = 600,
         bg = "transparent",
         file = output_filename)
}

#-------------------------------------------------------------
#Function for importing ROAST files
#-------------------------------------------------------------

read_roast <- function(file) {
  df           <- read.csv(file, row.names = 1)
  df$PropSelected <- ifelse(df$Direction == "Up", df$PropUp, -df$PropDown)
  df[, "PropSelected", drop = FALSE]
}


#-------------------------------------------------------------
#Function for computing standard errors
#-------------------------------------------------------------

std <- function(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))

