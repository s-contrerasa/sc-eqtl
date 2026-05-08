# 00 — Input Data

Input data are derived from Yazar et al. (2022) and obtained from the GEO repository (accession GSE196830).
File paths are system-specific and defined locally in `config/paths.env`.

---

## Gene Expression (Phenotypes)

| File | Format | Description |
|------|--------|-------------|
| `QC.RDS` | Seurat RDS | Quality-controlled single-cell RNA-seq object |

**Dataset:** OneK1K Phase 1 (Yazar et al. 2022, *Science*)
**Cells:** ~1.27 million PBMCs
**Donors:** ~981 individuals
**Pools:** 75 sequencing pools
**Cell types:** 22 immune cell subtypes grouped into 6 main lineages (CD4, CD8, NK, B cells, Monocyte, DC)

---

## Genotypes

| Files | Format | Description |
|-------|--------|-------------|
| `chr_1.bed/bim/fam` – `chr_22.bed/bim/fam` | PLINK binary | Genotype data split by chromosome (autosomes only) |

**Reference genome:** GRCh38

---

## Source Publication

Yazar S, et al. (2022). Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease. *Science*, 376(6589).
DOI: [10.1126/science.abf3041](https://doi.org/10.1126/science.abf3041)

Raw data available at: GEO accession [GSE196830](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196830)
