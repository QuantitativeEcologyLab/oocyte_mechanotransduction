# This script compares oocyte gene expression across three in vitro culture
# conditions: 2D, dynamic, and static. This includes:
# - DESeq2 differential expression analysis across culture conditions
# - ROAST gene set analysis on mechanotransduction pathways
# - Venn diagram of expressed genes across conditions (Figure 2a)
# - Heatmap of top 100 DEGs from 2D vs Dynamic contrast (Figure 2b)
# - PCA of gene expression across culture types (Figure 2c)
# - Volcano plots for each pairwise contrast (Figures 2d-f)
# - Heatmap of mechanotransduction-related gene expression (Figure 3a-b)

# Written by Michael Noonan
# Last updated: March 25, 2026

#Load in any necessary packages
library(ggvenn)
library(scico)
library(DESeq2)
library(plotly)
library(ggplot2)
library(edgeR)
library(pheatmap)


#Import the relevant datasets and project functions
source("scripts/00_data_import.R")
source("scripts/00_project_functions.R")

#-------------------------------------------------------------
#DESeq2 differential expression analysis across culture conditions
#-------------------------------------------------------------

#Drop gene metadata columns (gene_id, gene_name, gene_length)
count_data <- oocyte_expression[,-c(1:3)]

#Assign gene IDs as row names for DESeq2
rownames(count_data) <- oocyte_expression$gene_id

#Define condition factor with four replicates per group
conditions <- factor(c(rep("twoD", 4), rep("dynamic", 4), rep("static", 4)))

#Build sample metadata data frame for DESeq2 analysis
col_data <- data.frame(row.names = colnames(count_data), condition = conditions)

#Build DESeqDataSet with "condition" as the model term
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data,
 design = ~condition)
dds <- dds[rowSums(counts(dds)) > 1, ]

#Run the analysis
dds <- DESeq(dds)

res_2D_vs_Dyn <- results(dds, contrast = c("condition", "twoD", "dynamic"))
res_2D_vs_St <- results(dds, contrast = c("condition", "twoD", "static"))
res_Dyn_vs_St <- results(dds, contrast = c("condition", "dynamic", "static"))

#Save the results
write.csv(res_2D_vs_Dyn, file = "results/oocyte_DEG_Dynvs2D.csv", row.names = TRUE)
write.csv(res_2D_vs_St, file = "results/oocyte_DEG_Staticvs2D.csv", row.names = TRUE)
write.csv(res_Dyn_vs_St, file = "results/oocyte_DEG_DynvsStatic.csv", row.names = TRUE)

#Setting as 0.06 to identify genes with marginal significance
alpha <- 0.06
sig_2D_vs_Dyn <- res_2D_vs_Dyn[which(res_2D_vs_Dyn$padj < alpha), ]
sig_2D_vs_St <- res_2D_vs_St[ which(res_2D_vs_St$padj < alpha), ]
sig_Dyn_vs_St <- res_Dyn_vs_St[which(res_Dyn_vs_St$padj < alpha), ]

cat("Sig genes 2D vs Dyn:", nrow(sig_2D_vs_Dyn),
 "| 2D vs St:", nrow(sig_2D_vs_St),
 "| Dyn vs St:", nrow(sig_Dyn_vs_St), "\n")

#Merge gene names
sig_2D_vs_Dyn_df <- add_genename(sig_2D_vs_Dyn, oocyte_expression)
sig_2D_vs_St_df <- add_genename(sig_2D_vs_St, oocyte_expression)
sig_Dyn_vs_St_df <- add_genename(sig_Dyn_vs_St, oocyte_expression)


#-------------------------------------------------------------
#ROAST gene set analysis on mechanotransduction pathways
#-------------------------------------------------------------

#Build count matrix from oocyte_expression
all_genes <- unique(unlist(mechanotransduction_gene_sets))
all_gene_data <- oocyte_expression[oocyte_expression$gene_name %in% all_genes, ]

#Handle duplicate gene names by keeping highest mean expression
sample_cols <- grepl("^twoD_|^dynamic_|^static_", names(oocyte_expression))
sample_ids <- names(oocyte_expression)[sample_cols]
all_gene_data$mean_expr <- rowMeans(all_gene_data[, sample_ids])
all_gene_data <- all_gene_data[order(all_gene_data$gene_name, -all_gene_data$mean_expr), ]
all_gene_data <- all_gene_data[!duplicated(all_gene_data$gene_name), ]
all_gene_data$mean_expr <- NULL

#Build count matrix
counts_mat <- as.matrix(all_gene_data[, sample_ids])
rownames(counts_mat) <- toupper(all_gene_data$gene_name)

#Build group factor matching column order
group_roast <- factor(
 ifelse(grepl("^twoD_", sample_ids), "2D",
 ifelse(grepl("^dynamic_", sample_ids), "Dynamic",
 ifelse(grepl("^static_", sample_ids), "Static", NA)))
)

design_roast <- model.matrix(~0 + group_roast)
colnames(design_roast) <- levels(group_roast)

#Normalise and apply voom transformation
dge_roast <- DGEList(counts = counts_mat)
dge_roast <- calcNormFactors(dge_roast)
v_roast <- voom(dge_roast, design_roast)

gene_sets_upper <- lapply(mechanotransduction_gene_sets, toupper)

#ROAST function
run_oocyte_roast <- function(group_a, group_b, v, design, mechanotransduction_gene_sets, outfile) {
 cv <- numeric(ncol(design))
 cv[which(colnames(design) == group_a)] <- 1
 cv[which(colnames(design) == group_b)] <- -1
 res <- mroast(v$E, index = mechanotransduction_gene_sets, design = design, contrast = cv)
 write.csv(res, outfile, row.names = TRUE)
 res
}

#Run function across the three relevant comparisons
roast_Dynvs2D <- run_oocyte_roast("Dynamic", "2D", v_roast, design_roast, gene_sets_upper, "results/ROAST_Dynvs2D.csv")
roast_DynvsSt <- run_oocyte_roast("Dynamic", "Static", v_roast, design_roast, gene_sets_upper, "results/ROAST_DynvsSt.csv")
roast_2DvsSt <- run_oocyte_roast("2D", "Static", v_roast, design_roast, gene_sets_upper, "results/ROAST_2DvsSt.csv")


#-------------------------------------------------------------
#Figure 2 - Transcriptomic profiling of oocytes matured under different in vitro conditions
#-------------------------------------------------------------

#---------------
#Figure 2a - venn diagram of expressed genes
#---------------

gene_ids <- oocyte_expression$gene_id
group_2D <- oocyte_expression[, grepl("^twoD_", names(oocyte_expression))]
group_Dyn <- oocyte_expression[, grepl("^dynamic_", names(oocyte_expression))]
group_St <- oocyte_expression[, grepl("^static_", names(oocyte_expression))]

#Subset to only those genes found in >= two samples
expressed_2D <- gene_ids[rowSums(group_2D > 0) >= 2]
expressed_Dyn <- gene_ids[rowSums(group_Dyn > 0) >= 2]
expressed_St <- gene_ids[rowSums(group_St > 0) >= 2]

#Unique genes
unique_2D <- setdiff(expressed_2D, union(expressed_Dyn, expressed_St))
unique_Dyn <- setdiff(expressed_Dyn, union(expressed_2D, expressed_St))
unique_St <- setdiff(expressed_St, union(expressed_2D, expressed_Dyn))



venn_list <- list("2D" = expressed_2D, Dynamic = expressed_Dyn, Static = expressed_St)
ggvenn_plot <- 
 ggvenn(venn_list,
 fill_color = scico(3, palette = "batlow"),
 fill_alpha = 0.6,
 set_name_size = 5) +
 theme(plot.title = element_text(face = "bold"))
ggvenn_plot$layers[[3]]$aes_params$fontface <- "bold"
ggvenn_plot$layers[[4]]$data$text <- gsub("(\\d)(\\d{3})\n", "\\1,\\2\n", ggvenn_plot$layers[[4]]$data$text)

ggvenn_plot

ggsave("figures/figure_2a_venn_diagram.png", plot = ggvenn_plot, width = 8, height = 6, dpi = 600)


#---------------
#Figure 2b - heatmap of the top 100 most differentially expressed genes
#---------------

#Some data carpentry
norm_counts <- counts(dds, normalized = TRUE)
top100_genes <- rownames(sig_2D_vs_Dyn[order(sig_2D_vs_Dyn$padj), ])[1:100]
top100_counts <- norm_counts[top100_genes, ]
log_top100 <- log2(top100_counts + 1)

#Convert factor names
cond_factor <- factor(conditions)
levels(cond_factor) <- c("Dynamic", "Static", "2D")
annotation_col <- data.frame(Culture = cond_factor)
rownames(annotation_col) <- colnames(log_top100)

#Vector of colour names
COLS <- setNames(scico(3, palette = "batlow"), c("2D", "Dynamic", "Static"))


pheatmap(log_top100,
 cluster_rows = TRUE,
 annotation_col = annotation_col,
 annotation_colors = list(Culture = c(Dynamic = "#818231", Static = "#F9CCF9", `2D` = "#001959")),
 cluster_cols = FALSE,
 show_rownames = FALSE,
 show_colnames = FALSE,
 color = scico(1000, palette = 'lipari'),
 breaks = seq(0, 16, length.out = 1001),
 filename = "figures/figure_2b_heatmap_top100_genes_2D_vs_Dyn.png",
 width = 8, height = 10)


#---------------
#Figure 2c - PCA of gene expression across the 3 culture types
#---------------

#Apply a variance-stabilising transformation retaining experimental variation
vsd <- vst(dds, blind = FALSE)

#Conduct PCA on voom-normalised log-CPM values across all expressed genes
pca <- prcomp(t(assay(vsd)), scale. = FALSE)
percentVar <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)

#Build PCA data frame and join sample metadata
pca_df <- merge(
 cbind(as.data.frame(pca$x[, 1:2]), Sample = rownames(pca$x)),
 as.data.frame(colData(dds)),
 by.x = "Sample", by.y = "row.names"
)

#Convert factor names
levels(pca_df$condition) <- c("Dynamic", "Static", "2D")

#Define axis labels
DIM_1 <- paste0("PC1 (", percentVar[1], "% variance)")
DIM_2 <- paste0("PC2 (", percentVar[2], "% variance)")


FIG_2C <-
 ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
 geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
 geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1) +
 geom_point(size = 3) +
 scale_color_manual(values = c(Dynamic = "#818231", Static = "#F9CCF9", `2D` = "#001959")) +
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

ggsave(FIG_2C,
 file = "figures/figure_2c_PCA_culture_types.png",
 width = 6, height = 5,
 dpi = 600, bg = "transparent")


#---------------
#Figure 2d-f - Volcano plots of culture contrasts
#---------------


res_2D_vs_Dyn_df <- add_genename(res_2D_vs_Dyn, oocyte_expression)
res_2D_vs_St_df <- add_genename(res_2D_vs_St, oocyte_expression)
res_Dyn_vs_St_df <- add_genename(res_Dyn_vs_St, oocyte_expression)


create_volcano_plot(res_2D_vs_Dyn_df, "2D vs Dynamic", "figures/figure_2d_volcano_plot_2D_vs_Dyn.png")
create_volcano_plot(res_2D_vs_St_df, "2D vs Static", "figures/figure_2e_volcano_plot_2D_vs_St.png")
create_volcano_plot(res_Dyn_vs_St_df, "Dynamic vs Static", "figures/figure_2f_volcano_plot_Dyn_vs_St.png")



#-------------------------------------------------------------
#Figure 3 - Transcriptomic profiling of oocytes matured under different in vitro conditions
#-------------------------------------------------------------

#---------------
#Figure 3a
#---------------

#Identify any mechanotransduction genes present in the expression data
all_genes <- unique(unlist(mechanotransduction_gene_sets))
all_gene_data <- oocyte_expression[oocyte_expression$gene_name %in% all_genes, ]

#Extract sample columns
sample_cols <- grepl("^twoD_|^dynamic_|^static_", names(oocyte_expression))
sample_ids <- names(oocyte_expression)[sample_cols]

#Subset to sample columns plus gene_name
counts_sub <- all_gene_data[, c("gene_name", sample_ids)]

#Log_2 transform counts (helps to visualise)
counts_sub[, sample_ids] <- log2(counts_sub[, sample_ids] + 1)

counts_sub$mean_expr <- rowMeans(counts_sub[, sample_ids])
counts_sub <- counts_sub[order(counts_sub$gene_name, -counts_sub$mean_expr), ]
counts_sub <- counts_sub[!duplicated(counts_sub$gene_name), ]
counts_sub$mean_expr <- NULL

#Build expression matrix
expr_matrix <- as.data.frame(counts_sub[, sample_ids])
rownames(expr_matrix) <- counts_sub$gene_name

#Build annotation data frame
group <- ifelse(grepl("^twoD_", sample_ids), "2D",
 ifelse(grepl("^dynamic_", sample_ids), "Dynamic",
 ifelse(grepl("^static_", sample_ids), "Static", NA)))

annotation_col_heat <- data.frame(Culture = factor(group, levels = c("2D", "Dynamic", "Static")))
rownames(annotation_col_heat) <- sample_ids


pheatmap(as.matrix(expr_matrix),
 annotation_col = annotation_col_heat,
 annotation_colors = list(Culture = COLS),
 color = scico(1000, palette = 'lipari'),
 breaks = seq(0, 16, length.out = 1001),
 fontsize_row = 6,
 fontsize_col = 10,
 show_rownames = TRUE,
 cluster_rows = TRUE, 
 cluster_cols = FALSE,
 show_colnames = FALSE,
 filename = "figures/figure_3a_heatmap_mechanotransduction_Genes.png",
 width = 8, height = 14)


#---------------
#Figure 3b
#---------------


# Read and format ROAST results
r1 <- read_roast("results/ROAST_Dynvs2D.csv")
r2 <- read_roast("results/ROAST_2DvsSt.csv")
r3 <- read_roast("results/ROAST_DynvsSt.csv")

colnames(r1) <- "Prop_Dyn_vs_2D"
colnames(r2) <- "Prop_2D_vs_St"
colnames(r3) <- "Prop_Dyn_vs_St"

#Read FDR.Mixed values for significance annotation
fdr1 <- read.csv("results/ROAST_Dynvs2D.csv", row.names = 1)["FDR"]
fdr2 <- read.csv("results/ROAST_2DvsSt.csv", row.names = 1)["FDR"]
fdr3 <- read.csv("results/ROAST_DynvsSt.csv", row.names = 1)["FDR"]

colnames(fdr1) <- "FDR_Dyn_vs_2D"
colnames(fdr2) <- "FDR_2D_vs_St"
colnames(fdr3) <- "FDR_Dyn_vs_St"
#Merge proportion and FDR matrices
merged_roast <- merge(merge(r1, r2, by = "row.names"), r3, by.x = "Row.names", by.y = "row.names")
rownames(merged_roast) <- merged_roast$Row.names
merged_roast <- merged_roast[, -1]

merged_fdr <- merge(merge(fdr1, fdr2, by = "row.names"), fdr3, by.x = "Row.names", by.y = "row.names")
rownames(merged_fdr) <- merged_fdr$Row.names
merged_fdr <- merged_fdr[, -1]

heatmap_mat <- as.matrix(merged_roast)
fdr_mat <- as.matrix(merged_fdr[rownames(heatmap_mat), ])

#Remove any rows with NAs explicitly before building annotation matrix
keep_rows <- complete.cases(heatmap_mat)
heatmap_mat <- heatmap_mat[keep_rows, ]
fdr_mat <- fdr_mat[rownames(heatmap_mat), ]

#Build significance annotation matrix
num_mat_3b <- matrix(
 paste0(sprintf("%.2f", heatmap_mat),
 ifelse(fdr_mat < 0.05, " *",
 ifelse(fdr_mat < 0.09, " .", ""))),
 nrow = nrow(heatmap_mat),
 dimnames = dimnames(heatmap_mat)
)

#Generate the heatmap
pheatmap(heatmap_mat,
 color = scico(1000, palette = "vik", alpha = 0.8),
 breaks = seq(-1, 1, length.out = 1001),
 cluster_rows = FALSE,
 cluster_cols = FALSE,
 fontsize_row = 8,
 fontsize_col = 10,
 border_color = NA,
 display_numbers = num_mat_3b,
 number_color = "black",
 filename = "figures/figure_3b_ROAST_Heatmap.png")




#-------------------------------------------------------------
#Figure S6 - Heatmap with uncorrected p-values for discussion
#-------------------------------------------------------------


# Read and format ROAST results
r1 <- read_roast("results/ROAST_Dynvs2D.csv")
r2 <- read_roast("results/ROAST_2DvsSt.csv")
r3 <- read_roast("results/ROAST_DynvsSt.csv")

colnames(r1) <- "Prop_Dyn_vs_2D"
colnames(r2) <- "Prop_2D_vs_St"
colnames(r3) <- "Prop_Dyn_vs_St"

#Read FDR.Mixed values for significance annotation
fdr1 <- read.csv("results/ROAST_Dynvs2D.csv", row.names = 1)["PValue"]
fdr2 <- read.csv("results/ROAST_2DvsSt.csv", row.names = 1)["PValue"]
fdr3 <- read.csv("results/ROAST_DynvsSt.csv", row.names = 1)["PValue"]

colnames(fdr1) <- "FDR_Dyn_vs_2D"
colnames(fdr2) <- "FDR_2D_vs_St"
colnames(fdr3) <- "FDR_Dyn_vs_St"
#Merge proportion and FDR matrices
merged_roast <- merge(merge(r1, r2, by = "row.names"), r3, by.x = "Row.names", by.y = "row.names")
rownames(merged_roast) <- merged_roast$Row.names
merged_roast <- merged_roast[, -1]

merged_fdr <- merge(merge(fdr1, fdr2, by = "row.names"), fdr3, by.x = "Row.names", by.y = "row.names")
rownames(merged_fdr) <- merged_fdr$Row.names
merged_fdr <- merged_fdr[, -1]

heatmap_mat <- as.matrix(merged_roast)
fdr_mat <- as.matrix(merged_fdr[rownames(heatmap_mat), ])

#Remove any rows with NAs explicitly before building annotation matrix
keep_rows <- complete.cases(heatmap_mat)
heatmap_mat <- heatmap_mat[keep_rows, ]
fdr_mat <- fdr_mat[rownames(heatmap_mat), ]

#Build significance annotation matrix
num_mat_3b <- matrix(
  paste0(sprintf("%.2f", heatmap_mat),
         ifelse(fdr_mat < 0.05, " *",
                ifelse(fdr_mat < 0.09, " .", ""))),
  nrow = nrow(heatmap_mat),
  dimnames = dimnames(heatmap_mat)
)

#Generate the heatmap
pheatmap(heatmap_mat,
         color = scico(1000, palette = "vik", alpha = 0.8),
         breaks = seq(-1, 1, length.out = 1001),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 8,
         fontsize_col = 10,
         border_color = NA,
         display_numbers = num_mat_3b,
         number_color = "black",
         filename = "figures/figure_S6_ROAST_heatmap_oocyte_uncorrected_significance.png")



