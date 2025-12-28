# ============================================================
# 01_load_data_create_dds.R
# Purpose:
#   - Load GEO raw counts matrix
#   - Load sample metadata
#   - Create DESeq2 dataset object (dds)
#
# Outputs (created in environment):
#   - results folder in the working directory
#   - data     (count matrix)
#   - metadata (metadata)
#   - dds      (DESeqDataSet)
# ============================================================


# ---- Sources ----
# Data source: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162698
#   - Series RNA-seq raw counts matrix: GSE162698_raw_counts_GRCh38.p13_NCBI.tsv.gz
# Metadata source: Built by the author (Sergey Dadoyan), see the README file for details


# ---- Libraries ----
library(DESeq2)
library(dplyr)
library(readr)
library(tibble)

# ---- Creating results folder ----

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ---- Load counts matrix ----
data <- read_tsv("data/GSE162698_raw_counts_GRCh38.p13_NCBI.tsv") %>%
  mutate(GeneID = as.character(GeneID)) %>%
  column_to_rownames("GeneID")


# ---- Load sample metadata ----
metadata <- read.csv("data/GSE162698_metadata.csv") %>%
  column_to_rownames("Sample")


# ---- Validate alignment between metadata and count matrix ----
stopifnot(all(rownames(metadata) == colnames(data)))


# ---- Create DESeq2 object ----
metadata$Donor <- factor(metadata$Donor)
metadata$Polarization <- factor(metadata$Polarization)

dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = metadata,
  design = ~ Polarization 
)

