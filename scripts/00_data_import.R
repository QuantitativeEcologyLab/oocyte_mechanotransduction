# This script imports the various datasets and runs some
# data carpentry steps to get them into the correct format for any subsequent analyses.
# It also defines important global variables to be used across all analyses

# Written by Michael Noonan
# last updated: March 24, 2026


#-------------------------------------------------------------
# Data import
#-------------------------------------------------------------

#---------------
#Import sample metadata for the published gene expression data
#---------------

meta_data <- read.csv("data/sample_metadata.csv")


#---------------
#Import the gene expression data
#---------------

# Note data are from Milazzotto et al. 2022 https://doi.org/10.1016/j.isci.2022.103904
# They have been scaled via PQN normalisation for details see linked reference
expression_data <- read.csv("data/published_gene_expression_data.csv")

#Some data carpentry to get the genes into the correct format for downstream analyses
gene_names <- expression_data[-c(13133:13134), 1]
expression_data <- expression_data[1:13132, -1]
expression_data <- apply(expression_data, 2, as.numeric)
rownames(expression_data) <- gene_names


#---------------
#Import the oocyte gene expression data from the experimental work
#---------------

raw_oocyte_expression <- read.csv("data/oocyte_gene_expression_data.csv")


#Filter out lowly expressed genes
oocyte_expression <- raw_oocyte_expression[
  rowSums(raw_oocyte_expression[, grepl("^twoD_", names(raw_oocyte_expression))] > 0) >= 2 |
    rowSums(raw_oocyte_expression[, grepl("^dynamic_", names(raw_oocyte_expression))] > 0) >= 2 |
    rowSums(raw_oocyte_expression[, grepl("^static_", names(raw_oocyte_expression))] > 0) >= 2,
]

#Remove raw data from memory
rm(raw_oocyte_expression)
#cat("Genes before filtering:", nrow(raw_oocyte_expression), "\n")
#cat("Genes after filtering:", nrow(raw_oocyte_expression), "\n")


#---------------
#Import the embryo gene expression data from the experimental work
#---------------

# Import embryo count data
embryo_data <- read.csv("data/embryo_gene_expression_data.csv")

# Build count matrix
sample_cols_embryo <- grepl("^twoD_|^dynamic_", names(embryo_data))
counts_mat_embryo <- as.matrix(embryo_data[, sample_cols_embryo])
mode(counts_mat_embryo) <- "numeric"
rownames(counts_mat_embryo) <- embryo_data$gene_name

# Define group factor
group_embryo <- factor(
  ifelse(grepl("^twoD_", colnames(counts_mat_embryo)), "2D",
         ifelse(grepl("^dynamic_", colnames(counts_mat_embryo)), "Dynamic", NA)),
  levels = c("2D", "Dynamic")
)


#---------------
#Import the IVF data
#---------------

ivf_data <- read.csv("data/IVF_rates.csv")


#-------------------------------------------------------------
# Define Gene Sets
#-------------------------------------------------------------

#This section builds lists of important gene sets used across all analyses.
#Gene sets were identified via:
# - Reactome pathways: https://reactome.org/
# - Di et al. 2023 https://doi.org/10.1038/s41392-023-01501-9
# - Vining and Mooney 2017 https://doi.org/10.1038/nrm.2017.108

#Full set of mechanotransduction genes
mechanotransduction_gene_sets <- list(
  Integrins = c("ITGA1","ITGA2","ITGA3","ITGA5","ITGA6","ITGA7","ITGA8",
                "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5"),
  I_FAK_ERK = c("PTK2","SRC","MAPK1","MAPK3","PXN","GRB2","TLN1","VCL",
                "ITGA1","ITGA2","ITGA3","ITGA5","ITGA6","ITGA7","ITGA8",
                "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5"),
  I_RhoA = c("RHOA","ROCK1","ROCK2","MYL12B","ACTB","LIMK1","LIMK2",
             "CFL1","CFL2","ITGA1","ITGA2","ITGA3","ITGA5","ITGA6",
             "ITGA7","ITGA8","ITGB1","ITGB2","ITGB3","ITGB4","ITGB5"),
  I_ILK = c("ILK","PARVA","TNS1","PXN","TLN1",
            "ITGA1","ITGA2","ITGA3","ITGA5","ITGA6","ITGA7","ITGA8",
            "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5"),
  P_PI3K_Akt = c("PIK3CA","PIK3CB","PIK3CD","AKT1","AKT2","AKT3","MTOR",
                 "PIEZO1","PIEZO2"),
  P_Calpain = c("CAPN1","CAPN2","PIEZO1","PIEZO2"),
  P_MMP = c("MMP1","MMP2","MMP3","MMP7","MMP9","MMP14","MMP15","MMP16",
            "PIEZO1","PIEZO2"),
  P_RhoA = c("RHOA","ROCK1","ROCK2","MYL12B","LIMK1","LIMK2","CFL1","CFL2",
             "PIEZO1","PIEZO2"),
  G_cAMP_PKA = c("GNAS","PRKACA","PRKACB","ADCY1","ADCY2","ADCY3","ADCY4","ADCY5",
                 "ADGRL1","ADGRL2","ADGRL3"),
  G_Ras_MAPK = c("HRAS","KRAS","NRAS","RAF1","MAPK1","MAPK3",
                 "ADGRL1","ADGRL2","ADGRL3"),
  G_RhoA = c("RHOA","ROCK1","ROCK2","LIMK1","LIMK2",
             "ADGRL1","ADGRL2","ADGRL3"),
  T_TRP = c("TRPA1","TRPC1","TRPC2","TRPC3","TRPC6","TRPV1","TRPV4",
            "TMC1","TMC2"),
  C_YAP_TAZ = c("YAP1","WWTR1","TAZ","TEAD1","TEAD2","TEAD3","TEAD4",
                "CDH1","CDH2","CDH11"),
  C_RhoA = c("RHOA","ROCK1","ROCK2","CDH1","CDH2","CDH11"),
  S_PI3K_Akt = c("PIK3CA","PIK3CB","AKT1","AKT2","AKT3","MTOR",
                 "TRPC1","TRPC3","TRPC6","TRPV4"),
  S_RhoA = c("RHOA","ROCK1","ROCK2","TRPC1","TRPC3","TRPC6","TRPV4"),
  TLR_NFkB = c("NFKB1","RELA","IKBKB","TLR1","TLR2","TLR4","TLR5","TLR6","TLR9"),
  TLR_MAPK = c("MAPK1","MAPK3","MAPK14","TLR1","TLR2","TLR4","TLR5","TLR6","TLR9"),
  YAP_TAZ = c("YAP1","WWTR1","TEAD1","TEAD2","TEAD3","TEAD4"),
  I_FAX_Paxillin = c("PTK2","SRC","MAPK1","MAPK3","PXN",
                     "ITGA1","ITGA2","ITGA3","ITGA5","ITGA6","ITGA7","ITGA8",
                     "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5"),
  ATC_Actin = c("ACTA1","ACTB","ACTG1","MYH9","MYL12A","MYL12B"),
  ATC_Tubulin = c("TUBA1A","TUBA4A","TUBB","TUBB3","TUBG1"),
  ATC_Actin_associated = c("VCL","TLN1","FLNA","PALLD","ZYX"),
  ATC_RhoA = c("RHOA","ROCK1","ROCK2","LIMK1","LIMK2","CFL1","CFL2"),
  Hippo = c("MST1","MST2","LATS1","LATS2","MOB1A","MOB1B"),
  PIEZO = c("PIEZO1","PIEZO2"),
  GPCRs = c("ADGRL1","ADGRL2","ADGRL3"),
  TMCs = c("TMC1","TMC2"),
  Cadherins = c("CDH1","CDH2","CDH11"),
  SACs = c("TRPC1","TRPC3","TRPC6","TRPV4"),
  TLRs = c("TLR1","TLR2","TLR4","TLR5","TLR6","TLR9"),
  Ion_channels = c("KCNK","ASIC")
)

# Gene sets for mechanotransduction receptors and cytoskeleton components
# used for generating supplementary figure 1
gene_groups_receptors <- list(
  Integrins = c("ITGA1","ITGA2","ITGA3","ITGA5","ITGA6","ITGA7","ITGA8",
                "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5"),
  ATC_Actin = c("ACTA1","ACTB","ACTG1","MYH9","MYL12A","MYL12B"),
  ATC_Tubulin = c("TUBA1A","TUBA4A","TUBB","TUBB3","TUBG1"),
  ATC_Actin_associated = c("VCL","TLN1","FLNA","PALLD","ZYX"),
  PIEZO = c("PIEZO1","PIEZO2"),
  GPCRs = c("ADGRL1","ADGRL2","ADGRL3"),
  TMCs = c("TMC1","TMC2"),
  Cadherins = c("CDH1","CDH2","CDH11"),
  SACs = c("TRPC1","TRPC3","TRPC6","TRPV4"),
  TLRs = c("TLR1","TLR2","TLR4","TLR5","TLR6","TLR9")#,
  #Ion_channels = c("KCNK","ASIC")
)


# Gene sets for mechanotransduction pathways
# used for generating supplementary figure 2
gene_groups_pathways <- list(
  FAK_ERK = c("PTK2","SRC","MAPK1","MAPK3","PXN","GRB2","TLN1","VCL"),
  RhoA = c("RHOA","ROCK1","ROCK2","MYL12B","ACTB","LIMK1","LIMK2","CFL1","CFL2"),
  ILK = c("ILK","PARVA","TNS1","PXN","TLN1"),
  PI3K_Akt = c("PIK3CA","PIK3CB","PIK3CD","AKT1","AKT2","AKT3","MTOR"),
  Calpain = c("CAPN1","CAPN2"),
  MMP = c("MMP1","MMP2","MMP3","MMP7","MMP9","MMP14","MMP15","MMP16"),
  cAMP_PKA = c("GNAS","PRKACA","PRKACB"),
  Ras_MAPK = c("HRAS","KRAS","NRAS","RAF1","MAPK1","MAPK3"),
  YAP_TAZ = c("YAP1","WWTR1","TAZ","TEAD1","TEAD2","TEAD3","TEAD4"),
  NFkB = c("NFKB1","RELA","IKBKB"),
  MAPK = c("MAPK1","MAPK3","MAPK14"),
  FAX_Paxillin = c("PTK2","SRC","MAPK1","MAPK3","PXN")
)


#Secondary mechanotransduction pathways
gene_groups_emb <- list(
  Integrins = c("ITGA1","ITGA2","ITGA3","ITGA5","ITGA6","ITGA7","ITGA8",
                "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5"),
  ATC_Actin = c("ACTA1","ACTB","ACTG1","MYH9","MYL12A","MYL12B"),
  ATC_Tubulin = c("TUBA1A","TUBA4A","TUBB","TUBB3","TUBG1"),
  ATC_Actin_associated = c("VCL","TLN1","FLNA","PALLD","ZYX"),
  PIEZO = c("PIEZO1","PIEZO2"),
  GPCRs = c("ADGRL1","ADGRL2","ADGRL3"),
  TMCs = c("TMC1", "TMC2", "TMC3", "TMC4", "TMC5",
           "TMC6", "TMC7", "TMC8"),
  Cadherins = c("CDH1","CDH2","CDH11"),
  SACs = c("TRPC1","TRPC3","TRPC6","TRPV4"),
  TLRs = c("TLR1","TLR2","TLR4","TLR5","TLR6","TLR9"),
  Ion_channels_family = c("KCNK1", "KCNK2", "KCNK3", "KCNK4", "KCNK5",
                          "KCNK6", "KCNK9", "KCNK10", "KCNK12", "KCNK15",
                          "KCNK16", "KCNK17", "KCNK18",
                          "ASIC1", "ASIC2", "ASIC3", "ASIC4"),
  FAK_ERK = c("PTK2","SRC","MAPK1","MAPK3","PXN","GRB2","TLN1","VCL"),
  RhoA = c("RHOA","ROCK1","ROCK2","MYL12B","ACTB","LIMK1","LIMK2","CFL1","CFL2"),
  ILK = c("ILK","PARVA","TNS1","PXN","TLN1"),
  PI3K_Akt = c("PIK3CA","PIK3CB","PIK3CD","AKT1","AKT2","AKT3","MTOR"),
  Calpain = c("CAPN1","CAPN2"),
  MMP = c("MMP1","MMP2","MMP3","MMP7","MMP9","MMP14","MMP15","MMP16"),
  cAMP_PKA = c("GNAS","PRKACA","PRKACB"),
  Ras_MAPK = c("HRAS","KRAS","NRAS","RAF1","MAPK1","MAPK3"),
  YAP_TAZ = c("YAP1","WWTR1","TAZ","TEAD1","TEAD2","TEAD3","TEAD4"),
  NFkB = c("NFKB1","RELA","IKBKB"),
  MAPK = c("MAPK1","MAPK3","MAPK14"),
  FAK_Paxillin = c("PTK2","SRC","MAPK1","MAPK3","PXN")
)
