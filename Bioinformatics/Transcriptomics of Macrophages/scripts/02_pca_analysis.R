# ============================================================
# 02_pca_analysis.R
# Purpose:
#   - VST transform
#   - Create PCA plots (PC1 vs PC2, PC3 vs PC4) with custom legend panel
#   - Extract and annotate PC loadings + plot top absolute loadings for 4 principal components
#
# Inputs (expected to already exist in environment):
#   - dds      (Created in 01_load_data_create_dds.R)
#   - metadata (Created in 01_load_data_create_dds.R)
#
# Outputs (written to working directory):
#   - PC1_vs_PC2.png (PC1 vs PC2 plot)
#   - PC3_vs_PC4.png (PC3 vs PC4 plot)
#   - PC_loadings.csv (PC loadings table)
#   - PC_loadings.png (Top 10 PC absolute loadings per PC bar plot)
# ============================================================


# ---- Libraries ----
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(matrixStats)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidytext)
library(readr)
library(ragg)
library(grid)

# ---- VST transform ----

## Variance Stabilization of counts

vsd <- vst(dds, blind = TRUE)


# ---- Principal Component Analysis ----

## Preparing data

ntop <- 500
rv <- matrixStats::rowVars(assay(vsd))
top_idx <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

pca <- prcomp(t(assay(vsd)[top_idx, ]), center = TRUE, scale. = FALSE)

pca_data <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4]
)

pca_data <- cbind(pca_data, metadata[rownames(pca$x), , drop = FALSE])

pc_variance <- (pca$sdev^2) / sum(pca$sdev^2) * 100


## Plot settings (colors / shapes) 

pol_cols <- c("M0"="red",
              "M1"="blue",
              "M2 (IL10)"="purple",
              "M2 (IL4)"="darkgreen",
              "TAM-like"= "gold")


don_shapes <- c(21, 22, 24)

pca_data$Polarization <- factor(pca_data$Polarization, levels = names(pol_cols))

pca_data$Donor <- factor(pca_data$Donor)  


# ============================================================
# PCA: PC1 vs PC2
# ============================================================


## Creating the main plot

p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Donor, fill = Polarization),
             stroke = 0.5, size = 5.5, color = "black",
             alpha = 0.7) +
  scale_fill_manual(values = pol_cols) +
  scale_shape_manual(values = don_shapes) +
  scale_y_continuous(limits = c(-50, 125), breaks = seq(-50, 125, by = 25)) +
  scale_x_continuous(limits = c(-50, 125), breaks = seq(-50, 125, by = 25)) +
  labs(
    x = paste0("PC1: ", formatC(pc_variance[1], format = "f", digits = 1), "%"),
    y = paste0("PC2: ", formatC(pc_variance[2], format = "f", digits = 1), "%")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.title.x.bottom = element_text(size = 22, face = "bold",
                                       margin = margin(t = 15, unit = "pt")),
    legend.position = "none"
  ) +
  coord_fixed(ratio = 1)

## Creating the legend

donor_levels <- levels(pca_data$Donor)
donor_nums <- setNames(as.character(seq_along(donor_levels)), donor_levels)

legend_df <- expand.grid(
  Donor = donor_levels,
  Polarization = levels(pca_data$Polarization),
  KEEP.OUT.ATTRS = FALSE
)

legend_panel <- ggplot(legend_df, aes(x = Donor, y = Polarization, shape = Donor, fill = Polarization)) +
  geom_point(size = 6, stroke = 0.5, color = "black", alpha = 0.7) + 
  scale_fill_manual(values = pol_cols) + 
  scale_shape_manual(values = don_shapes) + 
  scale_x_discrete(position = "top", labels = donor_nums) +
  labs(x = "Donor", y = NULL) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.text.x.top = element_text(size = 21, color = "black",
                                       margin = margin(b = -2, unit = "pt")), 
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x.top = element_text(size = 18, color = "black", face = "bold",
                                        margin = margin(b = 8, unit = "pt")),
        legend.position = "none")


## Final plot

PC1_vs_PC2 <- (p | legend_panel + plot_layout(heights = c(1, 2))) +
  plot_layout(widths = c(6, 1))
PC1_vs_PC2

ggsave("results/PC1_vs_PC2.png", plot = PC1_vs_PC2,
       width = 8, height = 7, units = "in", dpi = 600, device = ragg::agg_png)


# ============================================================
# PCA: PC3 vs PC4
# ============================================================


## Creating the main plot

p_2 <- ggplot(pca_data, aes(x = PC3, y = PC4)) +
  geom_point(aes(shape = Donor, fill = Polarization),
             stroke = 0.5, size = 5.5, color = "black",
             alpha = 0.7) +
  scale_fill_manual(values = pol_cols) +
  scale_shape_manual(values = don_shapes) +
  scale_y_continuous(limits = c(-30, 50), breaks = seq(-30, 50, by = 20)) +
  scale_x_continuous(limits = c(-30, 50), breaks = seq(-30, 50, by = 20)) +
  labs(
    x = paste0("PC3: ", formatC(pc_variance[3], format = "f", digits = 1), "%"),
    y = paste0("PC4: ", formatC(pc_variance[4], format = "f", digits = 1), "%")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.title.x.bottom = element_text(size = 22, face = "bold",
                                       margin = margin(t = 15, unit = "pt")),
    legend.position = "none"
  ) +
  coord_fixed(ratio = 1)


## Final plot

PC3_vs_PC4 <- (p_2 | legend_panel + plot_layout(heights = c(1, 2))) +
  plot_layout(widths = c(6, 1))
PC3_vs_PC4

ggsave("results/PC3_vs_PC4.png", plot = PC3_vs_PC4,
       width = 8, height = 7, units = "in", dpi = 600, device = ragg::agg_png)


# ============================================================
# PC loadings: extraction + annotation + export
# ============================================================

## Extracting PC loadings:

### Loadings matrix (genes x PCs) 

pca_loadings <- pca$rotation

gene_ids_raw <- rownames(pca_loadings)


###  Detect ID type in rownames 

is_ensembl <- all(grepl("^ENSG\\d+", gene_ids_raw))
is_numeric <- all(grepl("^\\d+$", gene_ids_raw))

keytype <- if (is_ensembl) "ENSEMBL" else if (is_numeric) "ENTREZID" else "SYMBOL"

### Strip Ensembl version suffix 

gene_key <- if (keytype == "ENSEMBL") sub("\\..*$", "", gene_ids_raw) else gene_ids_raw

### Annotation table 

key_col <- keytype

anno_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(gene_key),
  keytype = keytype,
  columns = unique(c(key_col, "ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
) %>%
  distinct()

### Collapse 1-to-many mappings so join won't duplicate rows

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

### Build wide table for 4 PCs

pcs <- intersect(c("PC1","PC2","PC3","PC4"), colnames(pca_loadings))
stopifnot(length(pcs) > 0)

load_wide <- as.data.frame(pca_loadings[, pcs, drop = FALSE]) %>%
  mutate(
    gene_id  = gene_ids_raw,
    gene_key = gene_key
  ) %>%
  left_join(anno, by = "gene_key") %>%
  mutate(
    SYMBOL_clean = dplyr::if_else(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, gene_id)
  ) %>%
  dplyr::select(gene_id, ENSEMBL, SYMBOL, GENENAME, SYMBOL_clean, starts_with("PC"))

### Saving loadings table

readr::write_csv(load_wide, "results/PC_loadings.csv")


# ============================================================
# PC loadings: top genes plot
# ============================================================

## TOP PC loadings bar plot

pcs_loadings <- intersect(c("PC1","PC2","PC3","PC4"), names(load_wide))
stopifnot(length(pcs_loadings) > 0)

plot_df <- load_wide %>%
  dplyr::select(SYMBOL_clean, all_of(pcs_loadings)) %>%
  pivot_longer(all_of(pcs_loadings), names_to = "PC", values_to = "loading") %>%
  mutate(abs_loading = abs(loading)) %>%
  group_by(PC) %>%
  slice_max(abs_loading, n = 10, with_ties = TRUE) %>%
  ungroup() %>%
  mutate(SYMBOL = reorder_within(SYMBOL_clean, abs_loading, PC))   

loadings_plot <- ggplot(plot_df, aes(SYMBOL, abs_loading)) +
  geom_col(fill = "#4E79A7", color = "black", linewidth = 0.35, width = 0.75) +
  coord_flip() +
  facet_wrap(~PC, scales = "free_y", axes = "all_x", axis.labels = "all_x") +
  scale_x_reordered() +
  labs(
    title = "Top 10 genes by absolute loading per PC",
    subtitle = "Ranked within each PC",
    x = "Genes",
    y = "Absolute Loading"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black", size = 35),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 29),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size = 27),
    axis.title.x.bottom = element_text(face = "bold", color = "black", size = 28,
                                       margin = margin(t = 8, unit = "pt")),
    axis.title.y = element_text(face = "bold", color = "black", size = 35),
    axis.text.x = element_text(color = "black", size = 18),
    axis.text.y = element_text(color = "black", face = "italic", size = 18),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.4),
    panel.spacing = unit(0.9, "lines")
  )


ggsave("results/PC_loadings.png", plot = loadings_plot,
       width = 12, height = 10, units = "in", dpi = 600, device = ragg::agg_png)
