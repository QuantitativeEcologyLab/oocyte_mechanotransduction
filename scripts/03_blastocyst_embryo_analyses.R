# This script compares blastocyst development rates and transcriptomic profiles
# of individual blastocysts derived from oocytes matured under 2D or dynamic
# culture conditions. This includes:
# - Linear mixed location-scale model comparing blastocyst and hatching rates (Figure 4a)
# - DEG analysis (limma/voom) on embryo gene expression between culture conditions
# - ROAST gene set analysis on mechanotransduction pathways in embryos
# - PCA of embryo gene expression (Figure 4b)
# - Heatmap of top 50 differentially expressed mechanotransduction genes (Figure 4c)
# - ROAST heatmap of mechanotransduction pathway enrichment in embryos (Figure 4d)

# Written by Michael Noonan
# Last updated: March 25, 2026

#Load in any necessary packages
library(nlme)
library(limma)
library(edgeR)
library(scico)

#Import the relevant datasets and project functions
source("scripts/00_data_import.R")
source("scripts/00_project_functions.R")


#-------------------------------------------------------------
# Comparison of IVF-rates across culture conditions
#-------------------------------------------------------------

#First get some summary statistics on the IVF data
summary_stats <- do.call(rbind, lapply(split(ivf_data, ivf_data$Group), function(x) {
 data.frame(
 Group = x$Group[1],
 n = nrow(x),
 Blastocyst_mean = mean(x$Blastocyst, na.rm = TRUE),
 Blastocyst_sd = sd(x$Blastocyst, na.rm = TRUE),
 Hatched_mean = mean(x$Hatched_blastocyst, na.rm = TRUE),
 Hatched_sd = sd(x$Hatched_blastocyst, na.rm = TRUE)
 )
}))
summary_stats

#Linear mixed, location-scale model on the blastocyst rate
#between 2D and dynamic cell cultures
mod_blast_var <- lme(Blastocyst ~ Group,
 random = ~1 | Replicate,
 weights = varIdent(form = ~1 | Group),
 data = ivf_data)

summary(mod_blast_var)
intervals(mod_blast_var)

#-------------------------------------------------------------
# Differentially expressed genes in embryos from 2D vs. dynamic cell cultures
#-------------------------------------------------------------

#Build design matrix
levels(group_embryo) <- c("twoD", "Dynamic")
design <- model.matrix(~0 + group_embryo)
colnames(design) <- levels(group_embryo)
contrast_matrix <- makeContrasts(Dynamic - twoD, levels = design)

#Filter low expression (CPM >= 1 in at least n_min samples)
n_min <- min(table(group_embryo))
keep <- rowSums(cpm(counts_mat_embryo) >= 1) >= n_min
counts_mat_f <- counts_mat_embryo[keep, ]

cat("Genes before filter:", nrow(counts_mat_embryo), "after:", nrow(counts_mat_f), "\n")

#Define thresholds
logFC_thresh <- 1
fdr_thresh <- 0.05

#TMM normalisation and voom transformation
dge_emb <- DGEList(counts = counts_mat_f)
dge_emb <- calcNormFactors(dge_emb, method = "TMM")
v_emb <- voom(dge_emb, design, plot = FALSE)

#Fit linear model and extract DE results
fit_emb <- lmFit(v_emb, design)
fit_emb <- contrasts.fit(fit_emb, contrast_matrix)
fit_emb <- eBayes(fit_emb, trend = TRUE)
res_emb <- topTable(fit_emb, number = Inf, adjust.method = "fdr")
res_emb$logP <- -log10(pmax(res_emb$P.Value, .Machine$double.xmin))

write.csv(res_emb, "results/embryo_differentially_expressed_genes_Dyn_vs_2D.csv", row.names = TRUE)

DEGs_emb <- res_emb[abs(res_emb$logFC) >= logFC_thresh & res_emb$adj.P.Val < fdr_thresh, ]
up_in_Dyn <- DEGs_emb[DEGs_emb$logFC > 0, ]
down_in_Dyn <- DEGs_emb[DEGs_emb$logFC < 0, ]

cat("Total genes tested:", nrow(res_emb), "\n",
 "DEGs:", nrow(DEGs_emb), "\n",
 "Up in Dyn:", nrow(up_in_Dyn), "| Down in Dyn:", nrow(down_in_Dyn), "\n")


#-------------------------------------------------------------
#ROAST gene set analysis on mechanotransduction pathways in embryos from 2D vs. dynamic cell cultures
#-------------------------------------------------------------

#Save original row names before uppercasing -- needed for heatmap gene matching later
original_rownames_emb <- rownames(v_emb$E)

#Uppercase row names for gene matching
rownames(v_emb$E) <- toupper(rownames(v_emb$E))

#Map gene sets to indices in the expression matrix
gene_sets_upper <- lapply(mechanotransduction_gene_sets, toupper)
gene_sets_idx_emb <- ids2indices(gene_sets_upper,
 rownames(v_emb$E),
 remove.empty = TRUE)

#Run ROAST
roast_emb <- mroast(v_emb$E,
 index = gene_sets_idx_emb,
 design = design,
 contrast = contrast_matrix)

#Add signed proportion and enrichment direction columns
roast_emb$SignedProp <- ifelse(roast_emb$Direction == "Up", roast_emb$PropUp, -roast_emb$PropDown)
roast_emb$EnrichedIn <- ifelse(roast_emb$Direction == "Up", "Dynamic", "2D")
roast_emb$Pathway <- rownames(roast_emb)

write.csv(roast_emb, "results/embryo_ROAST_Dynamic_vs_2D.csv", row.names = FALSE)


#-------------------------------------------------------------
# Figure 4 - Transcriptomic profiling of individual blastocysts derived from oocytes
#-------------------------------------------------------------

#---------------
#Figure 4a - IVF success rate boxplots
#---------------

FIG_4A <-
 ggplot(ivf_data, aes(x = Group, y = Blastocyst, fill = Group)) +
 geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5, colour = "black") +
 geom_jitter(aes(colour = Group), size = 2.5, width = 0.1, height = 0) +
 stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "black") +
 scale_fill_manual(values = c(Dyn = "#818231", `2D` = "#001959")) +
 scale_colour_manual(values = c(Dyn = "#818231", `2D` = "#001959")) +
 theme_bw() +
 theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 panel.border = element_blank(),
 panel.background = element_rect(fill = "transparent"),
 plot.background = element_rect(fill = "transparent", color = NA),
 legend.title = element_text(size = 12, family = "sans", face = "bold", hjust = 0.5),
 legend.text = element_text(size = 12, family = "sans", face = "bold"),
 legend.background = element_rect(fill = "transparent"),
 axis.line = element_line(linetype = "solid"),
 plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
 axis.title.x = element_text(size = 12, face = "bold"),
 axis.title.y = element_text(size = 12, face = "bold"),
 axis.text = element_text(size = 10),
 axis.text.y = element_text(size = 10),
 axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
 legend.position = "none") +
 labs(y = "Blastocyst rate (% of total oocytes)", x = "") +
 scale_y_continuous(limits = c(0, 30), expand = c(0, 0))

ggsave(FIG_4A,
 file = "figures/figure_4a_blastocyst_rates.png",
 width = 7, height = 4, dpi = 300)


#---------------
#Figure 4b - PCA on embryo gene expression between cultures
#---------------

#Conduct PCA on voom-normalised log-CPM values across all expressed genes
pca_emb <- prcomp(t(v_emb$E), scale. = FALSE)
percentVar <- round(100 * (pca_emb$sdev^2) / sum(pca_emb$sdev^2), 1)

#Build data frame based on the results
pca_df_emb <- data.frame(
 PC1 = pca_emb$x[, 1],
 PC2 = pca_emb$x[, 2],
 Group = group_embryo,
 Sample = colnames(v_emb$E)
)

#Define axis labels
DIM_1 <- paste0("PC1 (", percentVar[1], "% variance)")
DIM_2 <- paste0("PC2 (", percentVar[2], "% variance)")


FIG_4B <-
 ggplot(pca_df_emb, aes(x = PC1, y = PC2, color = Group)) +
 geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
 geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1) +
 geom_point(size = 3) +
 scale_color_manual(values = c(Dynamic = "#818231", twoD = "#001959"),
 labels = c("Dynamic", "2D")) +
 labs(x = DIM_1, y = DIM_2, color = "") +
 theme_bw() +
 theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 panel.border = element_rect(colour = "black", linewidth = 1),
 panel.background = element_rect(fill = "transparent"),
 plot.background = element_rect(fill = "transparent", color = NA),
 legend.background = element_rect(fill = "transparent"),
 legend.text = element_text(size = 10, face = "bold"),
 legend.position = "inside",
 legend.position.inside = c(0.85, 0.9),
 axis.title = element_text(size = 12, face = "bold"),
 axis.text = element_text(size = 10))

ggsave(FIG_4B,
 file = "figures/figure_4b_PCA_embryo.png",
 width = 9, height = 4.7,
 dpi = 600, bg = "transparent")


#---------------
#Figure 4c - heatmap of top 50 DEGs in embryo between cultures
#---------------

#Select top 50 mechanotransduction genes by adjusted p-value then fold change
all_genes_emb <- toupper(unique(unlist(mechanotransduction_gene_sets)))
res_emb_mechano <- res_emb[rownames(res_emb) %in% all_genes_emb, ]
top50_mechano <- rownames(res_emb_mechano)[
 order(res_emb_mechano$adj.P.Val, -abs(res_emb_mechano$logFC), na.last = NA)
][1:min(50, nrow(res_emb_mechano))]

#Build annotation data frame
levels(group_embryo) <- c("2D", "Dynamic")
ann_col_emb <- data.frame(Culture = group_embryo)
rownames(ann_col_emb) <- colnames(v_emb$E)

pheatmap(v_emb$E[top50_mechano, ],
 annotation_col = ann_col_emb,
 annotation_colors = list(Culture = c(Dynamic = "#818231", `2D` = "#001959")),
 cluster_cols = FALSE,
 cluster_rows = TRUE,
 show_rownames = TRUE,
 show_colnames = FALSE,
 color = scico(1000, palette = "lipari"),
 #breaks = seq(0, 16, length.out = 1001),
 filename = "figures/figure_4c_heatmap_top50_mechano_embryo.png",
 width = 8, height = 10)


#---------------
# Figure 4d - ROAST heatmap of mechanotransduction pathways in embryos
#---------------

#Build signed proportion matrix sorted by value
roast_emb_sorted <- roast_emb[order(roast_emb$SignedProp, decreasing = TRUE), ]
roast_mat <- matrix(roast_emb_sorted$SignedProp, ncol = 1)
rownames(roast_mat) <- rownames(roast_emb_sorted)
colnames(roast_mat) <- "Dynamic vs 2D"

#Build number annotation matrix with significance symbols
num_mat <- matrix(
 paste0(sprintf("%.2f", roast_mat),
 ifelse(roast_emb_sorted$FDR < 0.05, " *",
 ifelse(roast_emb_sorted$FDR < 0.09, " .", ""))),
 ncol = 1,
 dimnames = list(rownames(roast_mat), colnames(roast_mat))
)

pheatmap(roast_mat,
 color = scico(1000, palette = "vik", alpha = 0.8),
 breaks = seq(-1, 1, length.out = 1001),
 cluster_rows = FALSE,
 cluster_cols = FALSE,
 display_numbers = num_mat,
 number_color = "black",
 border_color = NA,
 fontsize_row = 10,
 fontsize_col = 11,
 filename = "figures/figure_4d_ROAST_heatmap_embryo.png",
 width = 4.5, height = 9)



#---------------
# Figure S8 - ROAST heatmap of mechanotransduction pathways in embryos
#---------------

#Build signed proportion matrix sorted by value
roast_emb_sorted <- roast_emb[order(roast_emb$SignedProp, decreasing = TRUE), ]
roast_mat <- matrix(roast_emb_sorted$SignedProp, ncol = 1)
rownames(roast_mat) <- rownames(roast_emb_sorted)
colnames(roast_mat) <- "Dynamic vs 2D"

#Build number annotation matrix with significance symbols
num_mat <- matrix(
  paste0(sprintf("%.2f", roast_mat),
         ifelse(roast_emb_sorted$PValue < 0.05, " *",
                ifelse(roast_emb_sorted$PValue < 0.09, " .", ""))),
  ncol = 1,
  dimnames = list(rownames(roast_mat), colnames(roast_mat))
)

pheatmap(roast_mat,
         color = scico(1000, palette = "vik", alpha = 0.8),
         breaks = seq(-1, 1, length.out = 1001),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = num_mat,
         number_color = "black",
         border_color = NA,
         fontsize_row = 10,
         fontsize_col = 11,
         filename = "figures/figure_S8_ROAST_heatmap_embryo_uncorrected_significance.png",
         width = 4.5, height = 9)
