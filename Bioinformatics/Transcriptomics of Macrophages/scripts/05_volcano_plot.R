# ============================================================
# 05_volcano_plot.R
# Purpose:
#   - Build volcano plot for DESeq2 results (M1 vs M0) with labeled top 10 up- and top 10 downregulated genes
#
# Inputs (expected to already exist in environment):
#   - deseq_tbl (Annotated DESeq2 output table created in 04_differential_expression_deseq2.R)
#
# Output files:
#   - volcano_plot.png
# ============================================================


# ---- Libraries ----

library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(ragg)


# ============================================================
# Prepare volcano plot data
# ============================================================

## Ensuring log2FC and Padj are numeric

deseq_volcano <- as.data.frame(deseq_tbl)
deseq_volcano$log2FoldChange <- as.numeric(deseq_tbl$log2FoldChange)
deseq_volcano$padj <- as.numeric(deseq_volcano$padj)

## Differential expression labels (M1 vs M0)

deseq_volcano$diffexpressed <- "No significant change"
deseq_volcano$diffexpressed[deseq_volcano$log2FoldChange >=  1 & deseq_volcano$padj <= 0.01] <- "Upregulated in M1"
deseq_volcano$diffexpressed[deseq_volcano$log2FoldChange <= -1 & deseq_volcano$padj <= 0.01] <- "Downregulated in M1"

deseq_volcano$diffexpressed <- factor(
  deseq_volcano$diffexpressed,
  levels = c("Upregulated in M1", "No significant change", "Downregulated in M1")
)


## Pick TOP 10 up- and TOP 10 downregulated (among significant ones)

top_up <- deseq_volcano %>%
  filter(!is.na(padj), padj <= 0.01, log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top_down <- deseq_volcano %>%
  filter(!is.na(padj), padj <= 0.01, log2FoldChange <= -1) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top_labels <- c(top_up$SYMBOL_clean, top_down$SYMBOL_clean)

deseq_volcano$delabel <- ifelse(deseq_volcano$SYMBOL_clean %in% top_labels, deseq_volcano$SYMBOL_clean, NA)

deseq_volcano <- deseq_volcano %>% filter(!is.na(padj))

label_df <- subset(deseq_volcano, !is.na(delabel))

# ============================================================
# Volcano plot
# ============================================================

volcano_plot <- ggplot(deseq_volcano, aes(x = log2FoldChange, y = -log10(padj), label = delabel)) +
  geom_point(size = 2.5, shape = 21, color = "black", fill = NA) +
  geom_point(data = subset(deseq_volcano, diffexpressed == "Upregulated in M1"),
             size = 2.5, shape = 21, color = "black", fill = "#FF3333") +
  geom_point(data = subset(deseq_volcano, diffexpressed == "Downregulated in M1"),
             size = 2.5, shape = 21, color = "black", fill = "#3399FF") +
  geom_point(data = subset(deseq_volcano, diffexpressed == "No significant change"),
             size = 2.5, shape = 21, color = "black", fill = "grey50") +
  theme_classic() +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 1) +
  geom_hline(yintercept = 2, linetype = "dotted", color = "black", linewidth = 1) +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(
      label = delabel,
      fill  = ifelse(log2FoldChange > 0, "up", "down")
    ),
    color = "black",
    box.padding = 0.7,
    point.padding = 0,
    segment.color = "black",
    segment.size = 0.4,
    max.overlaps = Inf,
    max.iter = 5000,
    size = 6,
    fontface = "bold.italic",
    direction = "y",
    nudge_x = ifelse(label_df$log2FoldChange > 0,
                     18.5 - label_df$log2FoldChange,
                     -18.5 - label_df$log2FoldChange),
    ylim = c(2, 100),
    segment.curvature = 0,
    segment.angle = 90
  ) +
  scale_fill_manual(
    values = c(
      "up"   = scales::alpha("lightpink", 0.7),
      "down" = scales::alpha("lightblue", 0.7)
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(-19, 19), breaks = seq(-15, 15, by = 3)) +
  scale_y_continuous(limits = c(0, 63), breaks = seq(0, 60, by = 10)) +
  labs(
    x = expression("Log"[2]*"FC"),
    y = expression("-log"[10]*"(P"[adj]*")"),
    title = "Differential expression in M1 with respect to M0") +
  theme(
    axis.title = element_text(size = 24, face = "bold", color ="black"),
    axis.text = element_text(size = 28, color = "black"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black", size = 32, 
                              margin = margin(b = 40, unit = "pt"))
  )

ggsave("results/volcano_plot.png", plot = volcano_plot,
       width = 12, height = 10, units = "in", dpi = 600, device = ragg::agg_png)
