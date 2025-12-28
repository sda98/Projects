# ============================================================
# 03_normalized_counts_table.R
# Purpose:
#   - Estimate size factors and compute DESeq2 normalized counts
#   - Detect gene ID type (ENSEMBL / ENTREZID / SYMBOL)
#   - Annotate genes using org.Hs.eg.db
#   - Export annotated normalized counts table
#
# Inputs (expected to already exist in environment):
#   - dds (Created in 01_load_data_create_dds.R)
#
# Output files:
#   - normalized_counts_annotated.csv (Anntated table of normailized counts)
# ============================================================


# ---- Libraries ----

library(DESeq2)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readr)


# ---- Compute DESeq2 normalized counts ----

## Normalized counts table

dds_sf <- estimateSizeFactors(dds)
sizeFactors(dds_sf) 
study_normalized_counts <- counts(dds_sf, normalized = TRUE)


# ---- Detect gene identifier type ----

## Gene IDs from normalized count matrix

gene_ids_raw <- rownames(study_normalized_counts)

## Detect ID type

is_ensembl <- all(grepl("^ENSG\\d+", gene_ids_raw))
is_numeric <- all(grepl("^\\d+$", gene_ids_raw))
keytype <- if (is_ensembl) "ENSEMBL" else if (is_numeric) "ENTREZID" else "SYMBOL"

## Strip Ensembl version if needed

gene_key <- if (keytype == "ENSEMBL") sub("\\..*$", "", gene_ids_raw) else gene_ids_raw
key_col <- keytype


# ---- Pull annotations (org.Hs.eg.db) ----

## Pull annotations

anno_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(gene_key),
  keytype = keytype,
  columns = unique(c(key_col, "ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
) %>%
  distinct()

## Collapse to 1 row per gene_key to prevent join duplication

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


# ---- Build annotated normalized-counts table ----

## Build final normalized-counts table

norm_counts_tbl <- as.data.frame(study_normalized_counts) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_key = gene_key) %>%
  left_join(anno, by = "gene_key") %>%
  mutate(
    SYMBOL_clean = if_else(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, gene_id)
  ) %>%
  dplyr::select(gene_id, ENSEMBL, SYMBOL, GENENAME, SYMBOL_clean, starts_with("GSM"))


# ---- Sanity check + export ----

nrow(study_normalized_counts) == nrow(norm_counts_tbl)

readr::write_csv(norm_counts_tbl, "results/normalized_counts_annotated.csv")
