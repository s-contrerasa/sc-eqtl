# Quality control for OneK1K single-cell RNA-seq data
# Pool-aware filtering using per-pool normalized distributions (Yazar et al. 2022)

## Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(bestNormalize)
})

## Paths -------------------------------------------------------------------
source("config/paths.env")
input_file <- EXPR_INPUT
output_dir <- QC_DIR
fig_dir    <- file.path(output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Load data ---------------------------------------------------------------
cat("Loading pre-QC Seurat object...\n")
data <- readRDS(input_file)
stopifnot("pool_number" %in% colnames(data@meta.data))
cat("Cells before QC:", ncol(data), "\n")

## QC metrics --------------------------------------------------------------
# percent.mt is excluded: mitochondrial genes are absent from the expression
# matrix in this dataset, so the metric cannot be computed.
meta.data <- data@meta.data
metrics <- c("nCount_RNA", "nFeature_RNA")

## Helper functions --------------------------------------------------------

getCutoffs <- function(object, group, metric, lowSD = 3, highSD = 2) {
  split_obj <- object %>% split(.[[group]])

  valid <- map_lgl(split_obj, ~ {
    x <- pull(.x, !!sym(metric))
    sum(!is.na(x)) >= 10 && sd(x, na.rm = TRUE) > 0
  })
  split_obj <- split_obj[valid]
  if (length(split_obj) == 0) stop(paste("No valid pools for metric:", metric))

  res      <- split_obj %>% map(~ pull(.x, !!sym(metric))) %>% map(orderNorm)
  lower_z  <- res %>% map(~ mean(.$x.t) - sd(.$x.t) * lowSD)
  higher_z <- res %>% map(~ mean(.$x.t) + sd(.$x.t) * highSD)
  lower    <- map2(res, lower_z,  predict, inverse = TRUE)
  higher   <- map2(res, higher_z, predict, inverse = TRUE)

  tibble(
    pool_number = names(res),
    metric      = metric,
    lower       = unlist(lower),
    higher      = unlist(higher)
  )
}

splitCutoffs <- function(df) split(df, df$pool_number)

findOutliers <- function(cutoffs, metric, group, meta) {
  map_dfr(names(cutoffs), function(p) {
    pool_meta <- meta %>%
      rownames_to_column("barcode") %>%
      filter(.data[[group]] == p)

    bind_rows(
      pool_meta %>% filter(.data[[metric]] < cutoffs[[p]]$lower)  %>% mutate(type = "lower"),
      pool_meta %>% filter(.data[[metric]] > cutoffs[[p]]$higher) %>% mutate(type = "higher")
    )
  })
}

## Compute cutoffs ---------------------------------------------------------
cat("Computing per-pool QC cutoffs...\n")
cutoffs <- map_dfr(metrics, ~ getCutoffs(meta.data, "pool_number", .))

## Identify and remove outliers --------------------------------------------
cutoffs_by_pool <- splitCutoffs(cutoffs)

outliers <- map_dfr(
  metrics,
  ~ findOutliers(cutoffs_by_pool, ., "pool_number", meta.data),
  .id = "origin"
) %>% distinct(barcode, .keep_all = TRUE)

cat("Outlier cells identified:", nrow(outliers), "\n")

data_qc <- subset(data, cells = setdiff(colnames(data), outliers$barcode))
cat("Cells after QC:", ncol(data_qc), "\n")

## Latent variable (pool batch grouping) -----------------------------------
b1 <- levels(data_qc@meta.data$pool_number)[1:33]
data_qc$latent <- factor(ifelse(data_qc@meta.data$pool_number %in% b1, "b1", "b2"))

## QC figures --------------------------------------------------------------
plot_qc_metric <- function(meta_pre, meta_post, metric, label) {
  df <- bind_rows(
    meta_pre  %>% transmute(value = .data[[metric]], stage = "Pre-QC"),
    meta_post %>% transmute(value = .data[[metric]], stage = "Post-QC")
  ) %>% mutate(stage = factor(stage, levels = c("Pre-QC", "Post-QC")))

  ggplot(df, aes(x = stage, y = value, fill = stage)) +
    geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
    scale_fill_manual(values = c("Pre-QC" = "#d9534f", "Post-QC" = "#5cb85c")) +
    labs(title = label, x = NULL, y = label) +
    theme_classic() +
    theme(legend.position = "none")
}

p1 <- plot_qc_metric(meta.data, data_qc@meta.data, "nCount_RNA",   "UMI count per cell")
p2 <- plot_qc_metric(meta.data, data_qc@meta.data, "nFeature_RNA", "Genes detected per cell")

p_combined <- p1 + p2 + patchwork::plot_annotation(
  title = "OneK1K — QC filtering summary"
)

print(p_combined)
ggsave(file.path(fig_dir, "qc_pre_post.pdf"), p_combined, width = 8, height = 5)
ggsave(file.path(fig_dir, "qc_pre_post.png"), p_combined, width = 8, height = 5, dpi = 300)
cat("QC figure saved to", fig_dir, "\n")

## Save outputs ------------------------------------------------------------
cat("Saving outputs...\n")
saveRDS(data_qc,           file.path(output_dir, "QC.RDS"))
saveRDS(outliers,          file.path(output_dir, "outliers.RDS"))
saveRDS(cutoffs,           file.path(output_dir, "cutoffs.RDS"))
saveRDS(data_qc@meta.data, file.path(output_dir, "QC_metadata.RDS"))
saveRDS(
  tibble(before = ncol(data), after = ncol(data_qc), removed = ncol(data) - ncol(data_qc)),
  file.path(output_dir, "QC_summary.RDS")
)
cat("Done.\n")
