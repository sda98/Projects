# DESeq2 Differential Expression Analysis (GSE162698) — M1 vs M0

This repository contains an end-to-end **RNA-seq differential expression** workflow in **R** using **DESeq2** and the public GEO dataset **GSE162698**. The pipeline loads the raw counts matrix, builds a DESeq2 dataset, generates and annotates **normalized counts** performs **PCA**, runs differential expression (**M1 vs M0**), generates a **volcano plot**, and performs **GO + pathway enrichment** using **g:Profiler (gprofiler2)**. The outputs are saved in the **results** folder created in the working directory.

---

## Data sources

- **Raw counts (GEO Series raw counts matrix):** `GSE162698_raw_counts_GRCh38.p13_NCBI.tsv.gz` 
  Download page: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162698 (Also preloaded here in data/GSE162698_raw_counts_GRCh38.p13_NCBI.tsv.gz)

- **Sample metadata (curated):** `GSE162698_metadata.csv`
  This metadata table was **manually curated by the author** by reviewing GEO sample annotations (GSE/GSM pages) and standardizing them into an analysis-ready format.

---

## Metadata description (`GSE162698_metadata.csv`)

This file contains per-sample annotations used to define the DESeq2 design and to generate figures. The file is loaded here in data/GSE162698_metadata.csv

**Columns**
- `Sample`: GEO sample accession (e.g., `GSM...`). Must match the **column names of the counts matrix**.
- `Polarization`: condition label used in the DESeq2 design (e.g., `M0`, `M1`, `M2 (IL10)`, `M2 (IL4)`, `TAM-like`).
- `Donor`: donor identifier (used for visualization, e.g., point shapes in PCA plots).
- `Sample_Name`: human-readable sample label (optional; useful for plotting/figure captions).

---

## Repository structure

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
- Validates alignment between metadata and count matrix  
  (`stopifnot(all(rownames(metadata) == colnames(data)))`)
- Creates `dds` with design: `~ Polarization`

### 2) PCA analysis + PC loadings
**Script:** `scripts/02_pca_analysis.R`

- Variance stabilizing transformation: `vst(dds, blind = TRUE)`
- PCA on the top 500 most variable genes
- PCA plots:
  - `PC1_vs_PC2.png`
  - `PC3_vs_PC4.png`
- Extracts and annotates PC loadings:
  - `PC_loadings.csv`
  - `PC_loadings.png` (top absolute loadings per PC)

### 3) Normalized counts table
**Script:** `scripts/03_normalized_counts_table.R`

- Estimates size factors and computes DESeq2 normalized counts
- Detects gene ID type (ENSEMBL / ENTREZID / SYMBOL)
- Annotates genes via `org.Hs.eg.db`
- Exports: `normalized_counts_annotated.csv`

### 4) Differential expression (DESeq2)
**Script:** `scripts/04_differential_expression_deseq2.R`

- Runs differential expression: `dds_results <- DESeq(dds)`
- QC: dispersion plot (`plotDispEsts`) → `Dispersion_plot.png`
- Extracts DE results for: `contrast = c("Polarization", "M1", "M0")`
- Generates MA plots:
  - `MA_plot_unshrunk.png`
  - `MA_plot_shrunk.png` (LFC shrinkage via `ashr`)
- Annotates results via `org.Hs.eg.db`
- Exports: `DESeq2_results_annotated.csv`

### 5) Volcano plot
**Script:** `scripts/05_volcano_plot.R`

- Classifies genes by direction/significance (default: `padj <= 0.01`, `|log2FC| >= 1`)
- Labels top 10 upregulated + top 10 downregulated genes
- Exports: `volcano_plot.png`

### 6) GO enrichment (g:Profiler)
**Script:** `scripts/06_go_enrichment_gprofiler.R`

- GO enrichment for:
  - `GO:BP`, `GO:MF`, `GO:CC`
- Uses significant DE genes (Padj < 0.01, |Log2FC| > 1) as query and all tested genes as background
- Exports GO enrichment table: `gprofiler_GO.csv`
- Exports figure (dot plot + top-10 GO terms table): `GO_enrichment_with_table.png`
- Also produces an interactive `gostplot()` when run in-session

### 7) Pathway enrichment (g:Profiler)
**Script:** `scripts/07_pathway_enrichment_gprofiler.R`

- Pathway enrichment for:
  - `KEGG`, `REAC` (Reactome), `WP` (WikiPathways)
- Uses significant DE genes (Padj < 0.01, |Log2FC| > 1) as query and all tested genes as background
-  Exports pathways enrichment table: `gprofiler_pathways.csv`
- Exports figure (dot plot + top-10 pathways table): `Pathways_enrichments_with_table.png`
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
- `Dispersion_plot.png`
- `MA_plot_unshrunk.png`
- `MA_plot_shrunk.png`
- `volcano_plot.png`
- `gprofiler_GO.csv`
- `GO_enrichment_with_table.png`
- `gprofiler_pathways.csv`
- `Pathways_enrichments_with_table.png`
- `REAC_R-HSA-909733_genes_DESeq2.csv`

---

## Requirements

- **R** (tested with **R 4.5.2**)  
- Tested on **Windows 11** using **RStudio 2025.09.2 Build 418**

---

## How to run

1. Download or clone this repository.
2. Set your working directory to /Projects/Bioinformatics/Transcriptomics of Macrophages
3. Open R / RStudio.
4. Run scripts in order:

```r
source("scripts/01_load_data_create_dds.R")
source("scripts/02_pca_analysis.R")
source("scripts/03_normalized_counts_table.R")
source("scripts/04_differential_expression_deseq2.R")
source("scripts/05_volcano_plot.R")
source("scripts/06_go_enrichment_gprofiler.R")
source("scripts/07_pathway_enrichment_gprofiler.R")
                                      
