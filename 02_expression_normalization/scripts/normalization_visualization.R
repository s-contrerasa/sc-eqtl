# Normalization visualization
# Run interactively after sctransform_normalize.R and filter_celltypes.R.
# Produces variable feature plots and cell type composition summaries.

## Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

## Paths -------------------------------------------------------------------
source("config/paths.env")
norm_dir <- Sys.getenv("NORM_DIR")
fig_dir  <- file.path(norm_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Load filtered object ----------------------------------------------------
cat("Loading filtered Seurat object...\n")
data <- readRDS(file.path(norm_dir, "sctransform_norm_filtered.RDS"))

## 1. Variable features plot -----------------------------------------------
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)

p_vf_base  <- VariableFeaturePlot(data)
p_vf_label <- LabelPoints(plot = p_vf_base, points = top10, repel = TRUE)
p_vf <- p_vf_base + p_vf_label +
  plot_annotation(title = "Top 2000 variable features after SCTransform")

print(p_vf)
ggsave(file.path(fig_dir, "variable_features.pdf"), p_vf, width = 12, height = 5)
ggsave(file.path(fig_dir, "variable_features.png"), p_vf, width = 12, height = 5, dpi = 300)

## 2. Cell counts by main lineage ------------------------------------------
lineage_counts <- data@meta.data %>%
  count(main_lineage) %>%
  mutate(main_lineage = reorder(main_lineage, n))

p_lineage <- ggplot(lineage_counts, aes(x = main_lineage, y = n, fill = main_lineage)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = scales::comma(n)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma) +
  labs(title = "Cell counts by main immune lineage", x = NULL, y = "Number of cells") +
  theme_classic()

print(p_lineage)
ggsave(file.path(fig_dir, "cell_counts_main_lineage.pdf"), p_lineage, width = 7, height = 5)
ggsave(file.path(fig_dir, "cell_counts_main_lineage.png"), p_lineage, width = 7, height = 5, dpi = 300)

## 3. Cell counts by eQTL cell type ----------------------------------------
eqtl_counts <- data@meta.data %>%
  count(eqtl_celltype) %>%
  mutate(eqtl_celltype = reorder(eqtl_celltype, n))

p_eqtl <- ggplot(eqtl_counts, aes(x = eqtl_celltype, y = n, fill = eqtl_celltype)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = scales::comma(n)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma) +
  labs(title = "Cell counts by eQTL cell type (14 groups)", x = NULL, y = "Number of cells") +
  theme_classic()

print(p_eqtl)
ggsave(file.path(fig_dir, "cell_counts_eqtl_celltypes.pdf"), p_eqtl, width = 7, height = 6)
ggsave(file.path(fig_dir, "cell_counts_eqtl_celltypes.png"), p_eqtl, width = 7, height = 6, dpi = 300)

## 4. Cell type composition per donor (main lineage) -----------------------
donor_composition <- data@meta.data %>%
  count(donor_id, main_lineage) %>%
  group_by(donor_id) %>%
  mutate(pct = n / sum(n) * 100)

p_composition <- ggplot(donor_composition,
                        aes(x = donor_id, y = pct, fill = main_lineage)) +
  geom_col(width = 1) +
  labs(
    title = "Cell type composition per donor",
    x     = "Donor",
    y     = "% of cells",
    fill  = "Lineage"
  ) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

print(p_composition)
ggsave(file.path(fig_dir, "donor_celltype_composition.pdf"), p_composition, width = 12, height = 5)
ggsave(file.path(fig_dir, "donor_celltype_composition.png"), p_composition, width = 12, height = 5, dpi = 300)

cat("\nAll figures saved to", fig_dir, "\n")
