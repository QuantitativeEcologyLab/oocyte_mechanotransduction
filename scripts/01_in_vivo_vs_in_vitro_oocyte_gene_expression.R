# This script compares oocyte gene expression in vivo vs in vitro.
# This includes:
# - DEG analysis (limma) on bovine MII oocytes
# - ROAST analysis on mechanotransduction-specific genes 
# - Code for generating figures 1 and S1 (volcano plot, heatmaps, bar plots per gene group)

# Written by Michael Noonan
# last updated: March 25, 2026

#Load in any necessary packages
library(limma)
library(ggplot2)
library(pheatmap)
library(scico)

#Import the relevant datasets and project functions
source("scripts/00_data_import.R")
source("scripts/00_project_functions.R")



#-------------------------------------------------------------
# DEG analysis (limma) on bovine MII oocytes
#-------------------------------------------------------------

#Subset to gene expression for bovine-MII
cow_MII_IDs <- which(meta_data$Species == "cow" & meta_data$Stage == "MII")
cow_MII <- expression_data[, cow_MII_IDs]
IDs_cowMII <- meta_data[cow_MII_IDs, ]
group <- factor(IDs_cowMII$Collection)

#define the contrasts for the limma analysis
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(vitro - vivo, levels = design)

#Identify the differentially expressed genes via limma
#Note: p-values were adjusted for multiple comparisons
fit <- lmFit(cow_MII, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust = "fdr", number = Inf)
results$logP <- -log10(results$P.Value)
results$Significance <- ifelse(
 results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Significant", "Not Significant"
)

#Save the output
write.csv(results, "results/DEG_results_in_vivo_vs_in_vitro_bovine_MII_oocytes.csv", row.names = TRUE)



# Define key thresholds
logFC_thresh <- 1
pval_thresh <- 0.05

#Identify genes that meet the threshold criteria
DEGs <- results[abs(results$logFC) >= logFC_thresh & results$adj.P.Val < pval_thresh, ]
upregulated <- DEGs[DEGs$logFC > 0, ]
downregulated <- DEGs[DEGs$logFC < 0, ]

#Return key results
cat("Total genes tested:", nrow(results), "\n")
cat("Total DEGs:", nrow(DEGs), "\n")
cat("Upregulated:", nrow(upregulated), "\n")
cat("Downregulated:", nrow(downregulated), "\n")



#-------------------------------------------------------------
# ROAST analysis on mechanotransduction-specific genes for MII phase
#-------------------------------------------------------------

#Subset to bovine MII samples
MII_IDs <- which(meta_data$Species == "cow" & meta_data$Stage == "MII")
cow_MII <- expression_data[, MII_IDs]

#Build design with collection as the only term
collection_MII <- as.numeric(as.factor(meta_data$Collection[MII_IDs]))
design_MII <- model.matrix(~collection_MII)

#Conduct pathway analyses
RES_MII <- list()
for (i in seq_along(mechanotransduction_gene_sets)) {
 gene_index <- which(gene_names %in% mechanotransduction_gene_sets[[i]])
 ROAST_MII <- mroast(cow_MII, gene_index, design_MII, contrast = 2)
 ROAST_MII$Stage <- "MII"
 ROAST_MII$Collection <- "vitro - vivo"
 ROAST_MII$Pathway <- names(mechanotransduction_gene_sets)[i]
 ROAST_MII$Species <- "Bovine"
 RES_MII[[length(RES_MII) + 1]] <- ROAST_MII
}

roast_results_MII <- na.omit(do.call(rbind, RES_MII))

#Get signed proportions
roast_results_MII$SignedProp[roast_results_MII$Direction == "Up"] <- roast_results_MII$PropUp[roast_results_MII$Direction == "Up"]
roast_results_MII$SignedProp[roast_results_MII$Direction == "Down"] <- -roast_results_MII$PropDown[roast_results_MII$Direction == "Down"]

#Save results
write.csv(roast_results_MII, "results/ROAST_Bovine_metaphaseII_vivo_vs_vitro.csv", row.names = FALSE)


#-------------------------------------------------------------
# ROAST analysis on mechanotransduction pathways in bovine blastocysts
#-------------------------------------------------------------

#Subset to bovine blastocyst samples
BL_IDs <- which(meta_data$Species == "cow" & meta_data$Stage == "BL")
cow_BL <- expression_data[, BL_IDs]
IDs_BL <- meta_data[BL_IDs, ]

#Build design with collection as the only term
collection_BL <- as.numeric(as.factor(IDs_BL$Collection))
design_BL <- model.matrix(~collection_BL)

#Conduct pathway analyses
RES_BL <- list()
for (i in seq_along(mechanotransduction_gene_sets)) {
 gene_index <- which(gene_names %in% mechanotransduction_gene_sets[[i]])
 ROAST_BL <- mroast(cow_BL, gene_index, design_BL, contrast = 2)
 ROAST_BL$Stage <- "BL"
 ROAST_BL$Collection <- "vitro - vivo"
 ROAST_BL$Pathway <- names(mechanotransduction_gene_sets)[i]
 ROAST_BL$Species <- "Bovine"
 RES_BL[[length(RES_BL) + 1]] <- ROAST_BL
}

roast_results_BL <- na.omit(do.call(rbind, RES_BL))

#Get signed proportions
roast_results_BL$SignedProp <- ifelse(
 roast_results_BL$Direction == "Up",
 roast_results_BL$PropUp,
 -roast_results_BL$PropDown
)

write.csv(roast_results_BL, "results/ROAST_bovine_blastocyst_vivo_vs_vitro.csv", row.names = FALSE)

#-------------------------------------------------------------
# Figure 1
#-------------------------------------------------------------

#---------------
#Figure 1a - Volcano plot of differentially expressed genes
#---------------

fig_1_a <- 
 ggplot(results, aes(x = logFC, y = logP, color = Significance)) +
 geom_point(alpha = 0.6, pch = 16) +
 scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey"), name = "") +
 geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
 geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
 theme_bw() +
 theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 panel.border = element_blank(),
 axis.line = element_line(linetype = "solid"),
 panel.background = element_rect(fill = "transparent"),
 plot.background = element_rect(fill = "transparent", color = NA),
 legend.position = "inside",
 legend.position.inside = c(0.8,0.95),
 legend.title = element_text(size=8, family = "sans", face = "bold"),
 legend.text = element_text(size=8, family = "sans", face = "bold"),
 legend.background = element_rect(fill = "transparent"),
 legend.key.size = unit(0.2, 'cm'),
 legend.spacing.y = unit(0.1, 'cm'),
 plot.title = element_text(hjust = .01, vjust = -4, size = 12, family = "sans", face = "bold"),
 axis.title.y = element_text(size=12, family = "sans", face = "bold"),
 axis.title.x = element_text(size=12, family = "sans", face = "bold"),
 axis.text.y = element_text(size=10, family = "sans"),
 axis.text.x = element_text(size=10, family = "sans"),
 strip.background=element_blank(),
 plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
 labs(x = "log2 Fold Change", y = "-log10(P-value)") +
 
 ggrepel::geom_text_repel(
 data = subset(results, adj.P.Val < 0.05 & abs(logFC) > 2),
 aes(label = rownames(subset(results, adj.P.Val < 0.05 & abs(logFC) > 2))),
 size = 4, max.overlaps = 10, show.legend = F, segment.color = NA
 )

ggsave(fig_1_a, width = 8, height = 6,
 units = "in",
 dpi = 600,
 bg = "transparent",
 file = "figures/figure_1a_volcano_plot_MII_cow.png")




#---------------
#Figure 1b - heatmap of top 50 genes
#---------------

levels(group) <- c("in vitro", "in vivo")
ann_col <- data.frame(Origin = group)
rownames(ann_col) <- colnames(cow_MII)
top_genes <- rownames(results)[1:50]

pheatmap(cow_MII[top_genes, ],
 annotation_col = ann_col,
 annotation_colors = list(Origin = c("in vivo" = "seagreen4", "in vitro" = "dodgerblue4")),
 #scale = "row",
 cluster_cols = FALSE,
 show_rownames = TRUE,
 show_colnames = FALSE,
 color = scico(1000, palette = 'lipari'),
 breaks = seq(0, 19, length.out = 1001),
 filename = "figures/figure_1b_heatmap_top50_genes_MII_cow.png",
 width = 8, height = 10) 



#---------------
#Figure 1c - heatmap of mechanotransduction-related genes
#---------------
GENE_LIST <- unique(na.omit(unlist(mechanotransduction_gene_sets)))
mechanotransduction_genes <- intersect(GENE_LIST, rownames(cow_MII))

pheatmap(cow_MII[mechanotransduction_genes, ],
 annotation_col = ann_col,
 annotation_colors = list(Origin = c("in vivo" = "seagreen4", "in vitro" = "dodgerblue4")),
 #scale = "row",
 cluster_cols = FALSE,
 show_rownames = TRUE,
 show_colnames = FALSE,
 color = scico(1000, palette = 'lipari'),
 breaks = seq(0, 19, length.out = 1001),
 filename = "figures/figure_1c_heatmap_of_mechanotransduction_genes_MII_cow.png",
 width = 8, height = 14)



#-------------------------------------------------------------
# Figure S1 - Mechanotransduction receptor expression in vitro vs in vivo
#-------------------------------------------------------------

# Obtain expression data for the relevant genes
data_vv <- data.frame(expression_data[,meta_data$Species == "cow" & meta_data$Stage == "MII"])
colnames(data_vv) <- make.unique(meta_data[meta_data$Species == "cow" & meta_data$Stage == "MII","Collection"])
data_vv$GeneName <- rownames(data_vv)


#Convert from wide to long format for plotting
#Identify the columns starting with vivo or vitro
cols_to_pivot <- grep("^(vivo|vitro)", names(data_vv), value = TRUE)

long_data_vv <- reshape(
 data_vv,
 varying = cols_to_pivot,
 v.names = "NormalizedCounts",
 timevar = "SampleID",
 times = cols_to_pivot,
 direction = "long"
)

#Add Group column (using grepl logic)
long_data_vv$Group <- ifelse(grepl("^vivo", long_data_vv$SampleID), "vivo",
 ifelse(grepl("^vitro", long_data_vv$SampleID), "vitro", NA))

#Filter relevant genes
target_genes <- unique(na.omit(unlist(mechanotransduction_gene_sets)))
long_data_vv <- subset(long_data_vv, GeneName %in% target_genes)


group_colors <- c("vivo" = "seagreen4", "vitro" = "dodgerblue4")

# Pathway-level bar plots
for (grp_name in names(gene_groups_receptors)) {
 gene_list <- gene_groups_receptors[[grp_name]]
 group_data <- long_data_vv[long_data_vv$GeneName %in% gene_list,]
 
 summary_list <- by(group_data, list(group_data$Group, group_data$GeneName), function(x) {
 data.frame(Group = x$Group[1], 
 GeneName = x$GeneName[1],
 mean_count = mean(x$NormalizedCounts, na.rm = TRUE),
 sd_count = sd(x$NormalizedCounts, na.rm = TRUE),
 se_count = std(x$NormalizedCounts))
 })
 summary_data <- do.call(rbind, summary_list)
 
 p <- 
 ggplot() +
 geom_bar(data = summary_data,
 aes(x = GeneName, y = mean_count, fill = Group),
 stat = "identity", width = 0.6, position = position_dodge(0.8), alpha = 0.7) +
 geom_errorbar(data = summary_data,
 aes(x = GeneName, ymin = mean_count - se_count*1.96,
 ymax = mean_count + se_count*1.96, group = Group),
 width = 0.2, position = position_dodge(0.8)) +
 geom_jitter(data = group_data,
 aes(x = GeneName, y = NormalizedCounts, color = Group),
 position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
 size = 2, alpha = 0.9) +
 scale_fill_manual(values = group_colors, name = "Origin") +
 scale_color_manual(values = group_colors, name = "Origin") +
 labs(title = paste("Normalized Counts:", grp_name), x = "Gene", y = "Normalized Counts") +
 theme_bw() +
 theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 panel.background = element_rect(fill = "transparent"),
 plot.background = element_rect(fill = "transparent", color = NA),
 legend.title = element_text(size=12, family = "sans", face = "bold"),
 legend.text = element_text(size=12, family = "sans", face = "bold"),
 legend.background = element_rect(fill = "transparent"),
 legend.spacing.y = unit(0.1, 'cm'),
 axis.title.y = element_text(size=12, family = "sans", face = "bold"),
 axis.title.x = element_text(size=12, family = "sans", face = "bold"),
 axis.text.y = element_text(size=10, family = "sans"),
 axis.text.x = element_text(size=10, family = "sans",angle = 45, hjust = 1),
 strip.background=element_blank(),
 plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
 scale_x_discrete(expand = c(0,0)) +
 scale_y_continuous(expand = c(.01,0))
 
 ggsave(p, file = paste0("figures/figure_s1_panels/",grp_name, "_BarPlot_vivo_vitro.png"), width = 8, height = 6, dpi = 300)
}




#-------------------------------------------------------------
# Figure S2 - Mechanotransduction pathway expression in vitro vs in vivo
#-------------------------------------------------------------

#Pathway-level bar plots
for (grp_name in names(gene_groups_pathways)) {
 gene_list <- gene_groups_pathways[[grp_name]]
 group_data <- long_data_vv[long_data_vv$GeneName %in% gene_list,]
 
 summary_list <- by(group_data, list(group_data$Group, group_data$GeneName), function(x) {
 data.frame(Group = x$Group[1], 
 GeneName = x$GeneName[1],
 mean_count = mean(x$NormalizedCounts, na.rm = TRUE),
 sd_count = sd(x$NormalizedCounts, na.rm = TRUE),
 se_count = std(x$NormalizedCounts))
 })
 summary_data <- do.call(rbind, summary_list)
 
 p <- 
 ggplot() +
 geom_bar(data = summary_data,
 aes(x = GeneName, y = mean_count, fill = Group),
 stat = "identity", width = 0.6, position = position_dodge(0.8), alpha = 0.7) +
 geom_errorbar(data = summary_data,
 aes(x = GeneName, ymin = mean_count - se_count*1.96,
 ymax = mean_count + se_count*1.96, group = Group),
 width = 0.2, position = position_dodge(0.8)) +
 geom_jitter(data = group_data,
 aes(x = GeneName, y = NormalizedCounts, color = Group),
 position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
 size = 2, alpha = 0.9) +
 scale_fill_manual(values = group_colors, name = "Origin") +
 scale_color_manual(values = group_colors, name = "Origin") +
 labs(title = paste("Normalized Counts:", grp_name), x = "Gene", y = "Normalized Counts") +
 theme_bw() +
 theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 panel.background = element_rect(fill = "transparent"),
 plot.background = element_rect(fill = "transparent", color = NA),
 legend.title = element_text(size=12, family = "sans", face = "bold"),
 legend.text = element_text(size=12, family = "sans", face = "bold"),
 legend.background = element_rect(fill = "transparent"),
 legend.spacing.y = unit(0.1, 'cm'),
 axis.title.y = element_text(size=12, family = "sans", face = "bold"),
 axis.title.x = element_text(size=12, family = "sans", face = "bold"),
 axis.text.y = element_text(size=10, family = "sans"),
 axis.text.x = element_text(size=10, family = "sans",angle = 45, hjust = 1),
 strip.background=element_blank(),
 plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
 scale_x_discrete(expand = c(0,0)) +
 scale_y_continuous(expand = c(.01,0))
 
 ggsave(p, file = paste0("figures/figure_s2_panels/",grp_name, "_BarPlot_vivo_vitro.png"), width = 8, height = 6, dpi = 300)
}


#---------------
# Supplementary Figure S7 - ROAST heatmap of mechanotransduction pathways
#---------------

#Build signed proportion matrix sorted by value
roast_BL_sorted <- roast_results_BL[order(roast_results_BL$SignedProp, decreasing = TRUE), ]
roast_mat_BL <- matrix(roast_BL_sorted$SignedProp, ncol = 1)
rownames(roast_mat_BL) <- roast_BL_sorted$Pathway
colnames(roast_mat_BL) <- "In vivo vs In vitro"

#Build number annotation matrix with significance symbols
num_mat_BL <- matrix(
 paste0(sprintf("%.2f", roast_mat_BL),
 ifelse(roast_BL_sorted$FDR < 0.05, " *",
 ifelse(roast_BL_sorted$FDR < 0.09, " .", ""))),
 ncol = 1,
 dimnames = list(rownames(roast_mat_BL), colnames(roast_mat_BL))
)

pheatmap(roast_mat_BL,
 color = scico(1000, palette = "vik", alpha = 0.8),
 breaks = seq(-1, 1, length.out = 1001),
 cluster_rows = FALSE,
 cluster_cols = FALSE,
 display_numbers = num_mat_BL,
 number_color = "black",
 border_color = NA,
 fontsize_row = 10,
 fontsize_col = 11,
 fontsize_number = 12,
 filename = "figures/figure_s7_ROAST_heatmap_BL_vivo_vs_vitro.png",
 width = 4.5, height = 9)
