# QC visualization and summary
# Run this after scrna_qc.R to inspect QC results interactively.

## Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

## Paths -------------------------------------------------------------------
source("config/paths.env")
qc_dir  <- Sys.getenv("QC_DIR")
fig_dir <- file.path(qc_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Load QC outputs ---------------------------------------------------------
cat("Loading QC outputs...\n")
data_qc    <- readRDS(file.path(qc_dir, "QC.RDS"))
cutoffs    <- readRDS(file.path(qc_dir, "cutoffs.RDS"))
outliers   <- readRDS(file.path(qc_dir, "outliers.RDS"))
qc_summary <- readRDS(file.path(qc_dir, "QC_summary.RDS"))
meta_pre   <- readRDS(file.path(qc_dir, "preQC_singlets_metadata.RDS"))
meta_post  <- data_qc@meta.data

## 1. Summary table --------------------------------------------------------
cat("\n--- QC Summary ---\n")
print(qc_summary)

pool_summary <- meta_pre %>%
  rownames_to_column("barcode") %>%
  group_by(pool_number) %>%
  summarise(cells_before = n(), .groups = "drop") %>%
  left_join(
    meta_post %>%
      rownames_to_column("barcode") %>%
      group_by(pool_number) %>%
      summarise(cells_after = n(), .groups = "drop"),
    by = "pool_number"
  ) %>%
  mutate(
    cells_after   = replace_na(cells_after, 0),
    cells_removed = cells_before - cells_after,
    pct_removed   = round(100 * cells_removed / cells_before, 1)
  )

cat("\n--- Per-pool QC summary ---\n")
print(pool_summary, n = Inf)

## 2. Pre vs post violin plots ---------------------------------------------
plot_metric <- function(pre, post, metric, label) {
  bind_rows(
    pre  %>% transmute(value = .data[[metric]], stage = "Pre-QC"),
    post %>% transmute(value = .data[[metric]], stage = "Post-QC")
  ) %>%
    mutate(stage = factor(stage, levels = c("Pre-QC", "Post-QC"))) %>%
    ggplot(aes(x = stage, y = value, fill = stage)) +
    geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_manual(values = c("Pre-QC" = "#d9534f", "Post-QC" = "#5cb85c")) +
    labs(x = NULL, y = label) +
    theme_classic() +
    theme(legend.position = "none")
}

p_counts   <- plot_metric(meta_pre, meta_post, "nCount_RNA",   "UMI count per cell")
p_features <- plot_metric(meta_pre, meta_post, "nFeature_RNA", "Genes detected per cell")

p_violin <- (p_counts | p_features) +
  plot_annotation(title = "Pre vs Post QC — cell-level metrics")

print(p_violin)
ggsave(file.path(fig_dir, "qc_violin_pre_post.pdf"), p_violin, width = 8, height = 5)
ggsave(file.path(fig_dir, "qc_violin_pre_post.png"), p_violin, width = 8, height = 5, dpi = 300)

## 3. Per-pool cutoffs -----------------------------------------------------
p_cutoffs <- cutoffs %>%
  pivot_longer(c(lower, higher), names_to = "threshold", values_to = "value") %>%
  ggplot(aes(x = as.numeric(pool_number), y = value,
             color = threshold, linetype = threshold)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = c(lower = "#2171b5", higher = "#cb181d")) +
  labs(
    title = "Per-pool QC thresholds",
    x     = "Pool number",
    y     = "Threshold value",
    color = "Threshold", linetype = "Threshold"
  ) +
  theme_classic()

print(p_cutoffs)
ggsave(file.path(fig_dir, "qc_cutoffs_per_pool.pdf"), p_cutoffs, width = 10, height = 5)
ggsave(file.path(fig_dir, "qc_cutoffs_per_pool.png"), p_cutoffs, width = 10, height = 5, dpi = 300)

## 4. Cells removed per pool -----------------------------------------------
p_removed <- pool_summary %>%
  ggplot(aes(x = reorder(pool_number, -pct_removed), y = pct_removed)) +
  geom_col(fill = "#d9534f", alpha = 0.8) +
  labs(
    title = "Percentage of cells removed per pool",
    x     = "Pool",
    y     = "% cells removed"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))

print(p_removed)
ggsave(file.path(fig_dir, "qc_cells_removed_per_pool.pdf"), p_removed, width = 12, height = 5)
ggsave(file.path(fig_dir, "qc_cells_removed_per_pool.png"), p_removed, width = 12, height = 5, dpi = 300)

cat("\nAll figures saved to", fig_dir, "\n")
