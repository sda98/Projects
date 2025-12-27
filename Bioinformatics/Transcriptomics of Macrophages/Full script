# ============================================================
# 01_load_data_create_dds.R
# Purpose: Load count matrix + metadata and create DESeq2 dds
# ============================================================

# Data source: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162698 (Series RNA-seq raw counts matrix: GSE162698_raw_counts_GRCh38.p13_NCBI.tsv.gz)
# Metadata source: 

# Uploading packages

library(DESeq2)
library(dplyr)
library(readr)
library(tibble)

# Setting directions

setwd("C:/your/working/directory")
getwd()

# Loading data

data <- read_tsv("GSE162698_raw_counts_GRCh38.p13_NCBI.tsv") %>%
  mutate(GeneID = as.character(GeneID)) %>%
  column_to_rownames("GeneID")

# Loading metadata

metadata <- read.csv("Metadata.csv") %>%
  column_to_rownames("Sample")

# Matching data and metadata

stopifnot(all(rownames(metadata) == colnames(data)))

# Creating DESeq2 object

metadata$Donor <- factor(metadata$Donor)
metadata$Polarization <- factor(metadata$Polarization)

dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = metadata,
  design = ~ Polarization 
)

