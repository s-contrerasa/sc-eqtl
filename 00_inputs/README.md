# 00 — Input Data

This pipeline starts from two sets of input files provided by the Gao Lab (Purdue University).
These files are stored on the Negishi HPC cluster and are not included in this repository.

---

## Gene Expression (Phenotypes)

**Location:** `/depot/gao824/boran/onek1k/gene_expression/`

| File | Format | Description |
|------|--------|-------------|
| `QC.RDS` | Seurat RDS | Quality-controlled single-cell RNA-seq object. ~613K singlets across 75 sequencing pools and ~981 donors. |
| `onek1k.h5ad` | AnnData HDF5 | Raw expression data in Python-compatible format. |
| `onek1k.h5seurat` | H5Seurat | Raw expression data in Seurat-compatible HDF5 format. |

**Dataset:** OneK1K Phase 1 (Yazar et al. 2022, *Science*)
**Cells:** ~1.27 million PBMCs
**Donors:** ~981 individuals
**Pools:** 75 sequencing pools
**Cell types:** 22 immune cell subtypes grouped into 6 main lineages (CD4, CD8, NK, B cells, Monocyte, DC)

---

## Genotypes

**Location:** `/depot/gao824/boran/onek1k/plink_file/`

| Files | Format | Description |
|-------|--------|-------------|
| `chr_1.bed/bim/fam` – `chr_22.bed/bim/fam` | PLINK binary | Genotype data split by chromosome (autosomes only) |
| `combined_ids.txt` | Text | Sample IDs matching expression donors |

**Reference genome:** GRCh38
**Variants:** ~8M SNPs across 22 autosomes

---

## Source Publication

Yazar S, et al. (2022). Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease. *Science*, 376(6589).
DOI: [10.1126/science.abf3041](https://doi.org/10.1126/science.abf3041)

Raw data available at: GEO accession [GSE196830](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196830)
