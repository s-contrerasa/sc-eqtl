# SCTransform normalization of QC-filtered single-cell RNA-seq data (Yazar et al. 2022)
# Regresses out pool effects. percent.mt is excluded because mitochondrial genes
# are absent from this dataset's expression matrix.

## Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
})

## Paths -------------------------------------------------------------------
source("config/paths.env")
input_file <- file.path(Sys.getenv("QC_DIR"),   "QC.RDS")
output_dir <- Sys.getenv("NORM_DIR")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## Load data ---------------------------------------------------------------
cat("Loading QC-filtered Seurat object...\n")
data <- readRDS(input_file)
stopifnot("pool_number" %in% colnames(data@meta.data))
cat("Cells loaded:", ncol(data), "\n")

## SCTransform normalization -----------------------------------------------
cat("Running SCTransform...\n")
options(future.globals.maxSize = 1024^3 * 100)

data <- SCTransform(
  object          = data,
  vars.to.regress = "pool_number",
  conserve.memory = TRUE,
  verbose         = TRUE
)

## Save output -------------------------------------------------------------
cat("Saving normalized object...\n")
saveRDS(data, file.path(output_dir, "sctransform_norm.RDS"))
cat("Done.\n")
sessionInfo()
