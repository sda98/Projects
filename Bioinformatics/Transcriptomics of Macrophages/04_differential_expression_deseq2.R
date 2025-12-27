# ============================================================
# 04_differential_expression_deseq2.R
# Purpose:
#   - Run DESeq2 differential expression
#   - Visualize dispersion estimates
#   - Extract results for a specified contrast
#   - Apply LFC shrinkage (ashr) + MA plots
#   - Annotate results with org.Hs.eg.db
#   - Export annotated DE results table
#
# Inputs (expected to already exist in environment):
#   - dds (DESeqDataSet)
#
# Output files:
#   - DESeq2_results_annotated.csv
# ============================================================


# ---- Libraries ----
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readr)
library(ashr)


# ---- Working directory ----
setwd("C:/your/working/directory")
getwd()


# ============================================================
# Run DESeq2
# ============================================================

# Differential Expression Analysis

#DESeq2 

dds_results <- DESeq(dds)


# ============================================================
# QC: dispersion model fit
# ============================================================

#Disperison model

plotDispEsts(dds_results)


# ============================================================
# Extract DE results + MA plots + LFC shrinkage
# ============================================================

# Extracting the results

DESeq_results <- results(dds_results, contrast = c("Polarization", "M1", "M0"), alpha = 0.01, lfcThreshold = 1)

#MA plot and Log Shrinkage

plotMA(DESeq_results, ylim = c(-15,15))
DESeq_results <- lfcShrink(dds_results, contrast = c("Polarization", "M1", "M0"), res=DESeq_results, type = 'ashr')
plotMA(DESeq_results, ylim = c(-15,15))


# ============================================================
# Annotate DESeq2 results (org.Hs.eg.db)
# ============================================================

## Annotation

## Start from DESeq_results
res_tbl <- as.data.frame(DESeq_results) %>%
  rownames_to_column("gene_id")

## Detect ID type from result rownames
gene_ids_raw <- res_tbl$gene_id
is_ensembl <- all(grepl("^ENSG\\d+", gene_ids_raw))
is_numeric <- all(grepl("^\\d+$", gene_ids_raw))
keytype <- if (is_ensembl) "ENSEMBL" else if (is_numeric) "ENTREZID" else "SYMBOL"

## Build gene_key (strip Ensembl version if present)
gene_key <- if (keytype == "ENSEMBL") sub("\\..*$", "", gene_ids_raw) else gene_ids_raw
key_col <- keytype

## Annotation lookup
anno_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(gene_key),
  keytype = keytype,
  columns = unique(c(key_col, "ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
) %>% distinct()

anno <- anno_raw %>%
  group_by(.data[[key_col]]) %>%
  summarise(
    ENTREZID = paste(unique(na.omit(ENTREZID)), collapse = ";"),
    ENSEMBL  = paste(unique(na.omit(ENSEMBL)),  collapse = ";"),
    SYMBOL   = paste(unique(na.omit(SYMBOL)),   collapse = ";"),
    GENENAME = paste(unique(na.omit(GENENAME)), collapse = ";"),
    .groups = "drop"
  ) %>%
  rename(gene_key = !!key_col)

## Join annotations onto DESeq results
deseq_tbl <- res_tbl %>%
  mutate(gene_key = gene_key) %>%
  left_join(anno, by = "gene_key") %>%
  mutate(SYMBOL_clean = if_else(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, gene_id)) %>%
  dplyr::select(
    gene_id, ENSEMBL, SYMBOL, GENENAME, SYMBOL_clean,
    baseMean, log2FoldChange, lfcSE, pvalue, padj
  )

nrow(deseq_tbl) == nrow(res_tbl)

readr::write_csv(deseq_tbl, "DESeq2_results_annotated.csv")
