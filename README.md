# Oocyte Mechanotransduction

This repository contains the data and scripts required to reproduce all analyses and figures in the paper "Oocyte mechanotransduction is altered during *in vitro* maturation and can be enhanced through millifluidic mechanical stimulation" by Franko et al.

## Dependencies

All analyses were conducted in R version 4.3.3 (2024-02-29) on macOS 15.7.4. The following packages are required:
```r
install.packages(c(
  "nlme",        # v3.1-168 - linear mixed models
  "pheatmap",    # v1.0.13  - heatmaps
  "ggplot2",     # v3.5.2   - figures
  "ggvenn",      # v0.1.10  - Venn diagrams
  "ggrepel",     # v0.9.6   - label repulsion
  "plotly",      # v4.11.0  - interactive figures
  "scico",       # v1.5.0   - colour palettes
  "htmlwidgets", # v1.6.4   - saving interactive figures
  "ellipse"      # v0.5.0   - confidence ellipses for PCA
))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "limma",   # v3.58.1  - differential expression analysis
  "edgeR",   # v4.0.16  - count data normalisation
  "DESeq2"   # v1.42.1  - differential expression analysis
))
```

## Repository Structure
```
data/
    sample_metadata.csv                               # Sample metadata for published gene expression data
    published_gene_expression_data.csv                # PQN-normalised gene expression data (Milazzotto et al. 2022)
    oocyte_gene_expression_data.csv                   # Experimental oocyte RNA-seq count data
    embryo_gene_expression_data.csv                   # Experimental embryo RNA-seq count data
    IVF_rates.csv                                     # Blastocyst and hatching rate data
    Selected_Reactome_down_2DvsDyn.csv                # Reactome ORA results: downregulated in 2D vs Dynamic
    Reactome_2022_table_up_DynvsST.csv                # Reactome ORA results: upregulated in Dynamic vs Static
    Reactome_Pathways_2024_unique_Dyn.csv             # Reactome ORA results: unique Dynamic embryo genes
    Reactome_Pathways_2024_unique_2D.csv              # Reactome ORA results: unique 2D embryo genes

results/                                              # Output directory for analysis results (CSV files)

figures/                                              # Output directory for all figures
    figure_s1_panels/                                 # Individual panels for Supplementary Figure 1
    figure_s2_panels/                                 # Individual panels for Supplementary Figure 2

scripts/
    00_data_import.R
    00_project_functions.R
    01_in_vivo_vs_in_vitro_oocyte_gene_expression.R
    02_culture_type_analyses.R
    03_blastocyst_embryo_analyses.R
    04_supplementary_figures_4_and_5.R
```

## Data Availability

The experimental oocyte and embryo gene expression data, IVF rate data, and Reactome enrichment results required to reproduce all analyses are available in this repository.

The published gene expression data used in script 01 are from Milazzotto et al. (2022) https://doi.org/10.1016/j.isci.2022.103904 and are included here with permission. Raw RNA-seq data are available from the Gene Expression Omnibus (GEO).

## Workflow

Scripts should be run in the following order. Scripts 01 onward source `00_data_import.R` and `00_project_functions.R` automatically and do not need to be run separately.

| Script | Description |
|--------|-------------|
| `00_data_import.R` | Imports all datasets and defines gene sets used across all analyses |
| `00_project_functions.R` | Helper functions used across all scripts |
| `01_in_vivo_vs_in_vitro_oocyte_gene_expression.R` | DEG and ROAST analyses comparing in vivo vs in vitro bovine MII oocytes and blastocysts; generates Figures 1, S1, S2, and S6 |
| `02_culture_type_analyses.R` | DESeq2 DEG and ROAST analyses across 2D, dynamic, and static culture conditions; generates Figures 2, and 3|
| `03_blastocyst_embryo_analyses.R` | IVF rate modelling and DEG/ROAST analyses of embryo gene expression; generates Figure 4 |
| `04_supplementary_figures_4_and_5.R` | Reactome pathway enrichment visualizations; generates Figures S4 and S5 |

## Citation

Franko et al. (under review). Oocyte mechanotransduction is altered during *in vitro* maturation and can be enhanced through millifluidic mechanical stimulation.

## License

MIT