# DESeq2 Differential Expression Analysis (GSE162698) (M1 vs. M0)

This repository contains an end-to-end **DESeq2 RNA-seq differential expression** workflow in R using the public GEO dataset **GSE162698**. The pipeline loads the provided raw counts matrix, creates a DESeq2 object, performs QC and PCA, runs differential expression (**M1 vs M0**), generates a volcano plot, and performs GO + pathway enrichment via **g:Profiler**.

---

## Data sources

- **Raw counts (GEO Series raw counts matrix):** `GSE162698_raw_counts_GRCh38.p13_NCBI.tsv.gz`  
  Download page: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162698

- **Sample metadata (curated):** `GSE162698_metadata.csv`  
  This metadata table was **manually curated by the author** for this project by reading sample annotations on GEO (GSE/GSM pages) and standardizing them into a clean analysis-ready format (see details below).

---

## Metadata description (GSE162698_metadata.csv)

This file contains the per-sample annotations required to run the DESeq2 design and generate figures.

**Columns**
- `Sample`: GEO sample accession (e.g., `GSM...`). Must match the **column names of the counts matrix**.
- `Polarization`: condition label used in the DESeq2 design (e.g., `M0`, `M1`, `M2 (IL10)`, `M2 (IL4)`, `TAM-like`).
- `Donor`: donor identifier used for visualization (e.g., shapes in PCA plots).
- `Sample_Name`: human-readable sample label (for plotting/figure captions).

---

## Repo structre

- `data/`
  - `GSE162698_raw_counts_GRCh38.p13_NCBI.tsv`
  - `GSE162698_metadata.csv`
- `scripts/`
  - `01_load_data_create_dds.R`
  - `02_pca_analysis.R`
  - `03_normalized_counts_table.R`
  - `04_differential_expression_deseq2.R`
  - `05_volcano_plot.R`
  - `06_go_enrichment_gprofiler.R`
  - `07_pathway_enrichment_gprofiler.R`

---

## Pipeline overview

### 1) Load data + create DESeq2 object
**Script:** `scripts/01_load_data_create_dds.R`

- Loads raw counts table
- Loads curated metadata (`GSE162698_metadata.csv`)
- Validates sample order (`stopifnot(all(rownames(metadata) == colnames(data)))`)
- Builds `dds` with design: `~ Polarization`

### 2) PCA analysis + PC loadings
**Script:** `scripts/02_pca_analysis.R`

- VST transform (`vst(dds, blind = TRUE)`)
- PCA based on top 500 variable genes
- PC1 vs PC2 plot + custom legend panel (`PC1_vs_PC2.png`)
- PC3 vs PC4 plot + custom legend panel (`PC3_vs_PC4.png`)
- Extracts and annotates PC loadings → `PC_loadings.csv`
- Top loadings plot → `PC_loadings.png`

### 3) Normalized counts table
**Script:** `scripts/03_normalized_counts_table.R`

- Size factor estimation and normalized counts
- Gene ID detection (ENSEMBL / ENTREZID / SYMBOL)
- Annotation via `org.Hs.eg.db`
- Export: `normalized_counts_annotated.csv`

### 4) Differential expression (DESeq2)
**Script:** `scripts/04_differential_expression_deseq2.R`

- Runs `DESeq(dds)`
- Dispersion QC plot (`plotDispEsts`)
- Extracts results: `contrast = c("Polarization", "M1", "M0")`
- LFC threshold + shrinkage (`ashr`)
- Annotates results via `org.Hs.eg.db`
- Export: `DESeq2_results_annotated.csv`

### 5) Volcano plot
**Script:** `scripts/05_volcano_plot.R`

- Categorizes genes by direction/significance (`padj <= 0.01`, `|log2FC| >= 1`)
- Labels top 10 up + top 10 down genes
- Export: `volcano_plot.png`

### 6) GO enrichment (g:Profiler)
**Script:** `scripts/06_go_enrichment_gprofiler.R`

- GO enrichment: BP/MF/CC
- Custom background genes
- Export: `gprofiler_GO.csv`
- Publication figure (dot plot + top-10 table): `GO_enrichment_with_table.png`
- Also generates an interactive `gostplot()` in-session

### 7) Pathway enrichment (g:Profiler)
**Script:** `scripts/07_pathway_enrichment_gprofiler.R`

- Pathway enrichment: KEGG / Reactome / WikiPathways
- Export: `gprofiler_pathways.csv`
- Publication figure (dot plot + top-10 table): `Pathway_Enrichments_with_table.png`
- Extracts genes for a selected pathway (`REAC:R-HSA-909733`) and exports DE stats:
  - `REAC_R-HSA-909733_genes_DESeq2.csv`

---

## Outputs

Typical outputs written by the scripts:

- `PC1_vs_PC2.png`
- `PC3_vs_PC4.png`
- `PC_loadings.csv`
- `PC_loadings.png`
- `normalized_counts_annotated.csv`
- `DESeq2_results_annotated.csv`
- `volcano_plot.png`
- `gprofiler_GO.csv`
- `GO_enrichment_with_table.png`
- `gprofiler_pathways.csv`
- `Pathway_Enrichments_with_table.png`
- `REAC_R-HSA-909733_genes_DESeq2.csv`

---

## Requirements

This workflow is written in **R** (version 4.5.2) and optimized for running on **Windows 11** OS RStudio (version 2025.09.2 Build 418)

---

## How to run

1. Download this repository
2. Place the scripts, data, and metadata files in the working directory.
3. Open R / RStudio.
4. Run scripts in order:

```r
source("01_load_data_create_dds.R")
source("02_pca_analysis.R")
source("03_normalized_counts_table.R")
source("04_differential_expression_deseq2.R")
source("05_volcano_plot.R")
source("06_go_enrichment_gprofiler.R")
source("07_pathway_enrichment_gprofiler.R")


                                              
