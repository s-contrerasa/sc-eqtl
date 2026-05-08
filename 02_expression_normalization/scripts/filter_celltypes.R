# Cell type filtering and annotation
#
# This script filters the normalized Seurat object to the 22 canonical immune
# cell types present in the OneK1K dataset, then assigns each cell to:
#   - one of 14 fine-grained eQTL cell types (used for eQTL mapping)
#   - one of 6 main immune lineages (used for pseudobulk aggregation)
#
# Cell type labels come from the predicted.celltype.l2 field, which contains
# Azimuth reference-based annotations.
#
# Filtering happens AFTER normalization so that SCTransform uses all QC-passed
# cells for variance estimation before subsetting to immune lineages of interest.

## Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## Paths -------------------------------------------------------------------
source("config/paths.env")
input_file <- file.path(Sys.getenv("NORM_DIR"), "sctransform_norm.RDS")
output_dir <- Sys.getenv("NORM_DIR")

## Load normalized object --------------------------------------------------
cat("Loading normalized Seurat object...\n")
data <- readRDS(input_file)
cat("Total cells before filtering:", ncol(data), "\n")

## Step 1 — Filter to 22 canonical immune cell types -----------------------
# These are the cell types present in the OneK1K dataset that have sufficient
# cell numbers for downstream eQTL analysis.
celltypes_keep <- c(
  # T cells — CD4
  "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating",
  # T cells — CD8
  "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating",
  # NK cells
  "NK", "NK_CD56bright", "NK Proliferating",
  # B cells
  "B naive", "B intermediate", "B memory", "Plasmablast",
  # Monocytes
  "CD14 Mono", "CD16 Mono",
  # Dendritic cells
  "cDC1", "cDC2", "ASDC", "pDC"
)

data <- subset(data, subset = predicted.celltype.l2 %in% celltypes_keep)
cat("Cells after cell type filtering:", ncol(data),
    "(", round(100 * ncol(data) / ncol(readRDS(input_file)), 1), "% retained )\n")

## Step 2 — Annotate 14 eQTL cell types ------------------------------------
# Fine-grained groupings used for cell-type-specific eQTL mapping.
# Matches Table S1 in Yazar et al. 2022.
data@meta.data <- data@meta.data %>%
  mutate(eqtl_celltype = case_when(
    predicted.celltype.l2 %in% c("CD4 Naive", "CD4 TCM")        ~ "CD4_NC",
    predicted.celltype.l2 %in% c("CD4 TEM", "CD4 CTL")          ~ "CD4_ET",
    predicted.celltype.l2 == "CD4 Proliferating"                 ~ "CD4_SOX4",
    predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TCM")        ~ "CD8_NC",
    predicted.celltype.l2 == "CD8 TEM"                           ~ "CD8_ET",
    predicted.celltype.l2 == "CD8 Proliferating"                 ~ "CD8_S100B",
    predicted.celltype.l2 %in% c("NK", "NK_CD56bright")         ~ "NK",
    predicted.celltype.l2 == "NK Proliferating"                  ~ "NK_R",
    predicted.celltype.l2 %in% c("B naive", "B intermediate")   ~ "B_IN",
    predicted.celltype.l2 == "B memory"                          ~ "B_Mem",
    predicted.celltype.l2 == "Plasmablast"                       ~ "Plasma",
    predicted.celltype.l2 == "CD14 Mono"                         ~ "Monocyte_C",
    predicted.celltype.l2 == "CD16 Mono"                         ~ "Monocyte_NC",
    predicted.celltype.l2 %in% c("cDC1", "cDC2", "ASDC", "pDC") ~ "DC"
  ))

## Step 3 — Annotate 6 main immune lineages --------------------------------
# Coarser groupings used for pseudobulk aggregation and PEER factor computation.
data@meta.data <- data@meta.data %>%
  mutate(main_lineage = case_when(
    predicted.celltype.l2 %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM",
                                  "CD4 CTL", "CD4 Proliferating")       ~ "CD4",
    predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TCM", "CD8 TEM",
                                  "CD8 Proliferating")                   ~ "CD8",
    predicted.celltype.l2 %in% c("NK", "NK_CD56bright", "NK Proliferating") ~ "NK",
    predicted.celltype.l2 %in% c("B naive", "B intermediate",
                                  "B memory", "Plasmablast")             ~ "Bcells",
    predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono")              ~ "Monocyte",
    predicted.celltype.l2 %in% c("cDC1", "cDC2", "ASDC", "pDC")        ~ "DC"
  ))

## Verify no unassigned cells remain ---------------------------------------
n_unassigned <- sum(is.na(data@meta.data$eqtl_celltype))
if (n_unassigned > 0) warning(n_unassigned, " cells could not be assigned an eQTL cell type.")

cat("\nCell counts by main lineage:\n")
print(table(data@meta.data$main_lineage))

cat("\nCell counts by eQTL cell type:\n")
print(table(data@meta.data$eqtl_celltype))

## Save output -------------------------------------------------------------
cat("\nSaving filtered and annotated object...\n")
saveRDS(data, file.path(output_dir, "sctransform_norm_filtered.RDS"))
cat("Done.\n")
