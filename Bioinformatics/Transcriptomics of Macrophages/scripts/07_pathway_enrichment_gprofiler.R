# ============================================================
# 07_pathway_enrichment_gprofiler.R
# Purpose:
#   - Pathway enrichment using g:Profiler (KEGG, Reactome, WikiPathways)
#   - Export enrichment results table
#   - Generate publication-style pathway dot plot + top-10 table panel
#   - Extract genes for a selected pathway and export DE stats
#
# Inputs (expected to already exist in environment):
#   - sig_genes (character vector of significant gene symbols)
#   - bg_genes  (character vector of background gene symbols)
#   - deseq_tbl (data.frame with SYMBOL_clean, log2FoldChange, padj, etc.)
#
# Output files:
#   - gprofiler_pathways.csv
#   - Pathway_Enrichments_with_table.png
#   - REAC_R-HSA-909733_genes_DESeq2.csv
# ============================================================


# ---- Libraries ----
library(dplyr)
library(tidyr)
library(stringr)
library(gprofiler2)
library(ggplot2)
library(ggrepel)
library(scales)
library(patchwork)
library(gridExtra)
library(gtable)
library(grid)
library(readr)


# ============================================================
# Run g:Profiler pathway enrichment (KEGG / Reactome / WP)
# ============================================================

## Pathway analysis

gp_pathway <- gost(
  query = sig_genes,
  organism = "hsapiens",
  correction_method = "g_SCS",
  user_threshold = 0.01,
  sources = c("KEGG", "REAC", "WP"),
  custom_bg = bg_genes,
  evcodes = TRUE,
)

head(gp_pathway$result)
readr::write_csv(gp_pathway$result, "gprofiler_pathways.csv")

gostplot(gp_pathway, capped = FALSE, interactive = TRUE)


# ============================================================
# Static pathway plot (gostplot) + top terms table
# ============================================================

Path <- gostplot(gp_pathway, capped = FALSE, interactive = FALSE)

top_table <- gp_pathway$result %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank, source, term_id, term_name, p_value) %>%
  transmute(
    Rank,
    Source = source,
    `Term ID` = term_id,
    `Term name` = str_trunc(term_name, 70),
    `P (adj)` = formatC(p_value, format = "e", digits = 2)
  )

top_terms <- top_table$`Term ID`


# ============================================================
# Customize Path plot styling (outlined points + highlight top 10)
# ============================================================

# first point layer index
pt_idx <- which(vapply(GO$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]

# mapping for that layer
m <- GO$layers[[pt_idx]]$mapping

# keep original colors as fill (if color was mapped)
if (is.null(m$fill)) {
  if (!is.null(m$colour)) m$fill <- m$colour
  if (!is.null(m$color))  m$fill <- m$color
}

# remove colour mapping -> black outline
m$colour <- NULL
m$color  <- NULL

# remove size mapping -> constant size
m$size <- NULL

Path$layers[[pt_idx]]$mapping <- m

# outlined points with constant size
Path$layers[[pt_idx]]$aes_params$shape  <- 21
Path$layers[[pt_idx]]$aes_params$colour <- "black"
Path$layers[[pt_idx]]$aes_params$stroke <- 0.5
Path$layers[[pt_idx]]$aes_params$size   <- 5.5   # <- pick your constant point size

# palette for GO sources
my_pal <- c("KEGG" = "#dd4477", "REAC" = "#3366cc", "WP" = "#0099c6")
Path <- Path + scale_fill_manual(values = my_pal) + guides(size = "none")
# Base alpha for everything
base_alpha <- 0.5

# Add a rank column to GO$data for top 10
Path$data <- Path$data %>%
  left_join(
    gp_pathway$result %>%
      arrange(p_value) %>%
      slice_head(n = 10) %>%
      mutate(Rank = row_number()) %>%
      dplyr::select(term_id, Rank),
    by = "term_id"
  ) %>%
  mutate(alpha_pt = ifelse(!is.na(Rank), 1, base_alpha))

# Map alpha per point (and turn off alpha legend)
Path <- Path +
  aes(alpha = alpha_pt) +
  scale_alpha_identity(guide = "none")

# Add rank labels for the top 10
label_df <- Path$data %>% filter(!is.na(Rank))

# remove the fixed alpha so mapping can work
Path$layers[[pt_idx]]$aes_params$alpha <- NULL

# map alpha to alpha_pt on the point layer
Path$layers[[pt_idx]]$mapping$alpha <- rlang::sym("alpha_pt")

# use the alpha values as-is (no legend)
Path <- Path + scale_alpha_identity(guide = "none")

label_df <- Path$data %>% dplyr::filter(!is.na(Rank))


# ============================================================
# Publication-style Pathway plot (rank labels + theme)
# ============================================================

p_pub <- Path +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(x = order, y = logpval, label = Rank),
    inherit.aes = FALSE,
    size = 7.5,
    fontface = "bold",
    color = "black",
    fill = "grey95",               # background fill
    label.size = 0.6,              # border thickness
    label.r = unit(5, "pt"),       # corner roundness
    box.padding = 0.9,
    point.padding = 0.82,
    segment.color = "black",
    segment.size = 0.8,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 3
  ) +
  labs(title = "Enriched pathways", y = expression(-log[10](p[adj]))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 33, color = "black", hjust = 0, margin = margin(b = 6)),
    axis.title.y = element_text(face = "bold", size = 22, color = "black", margin = margin(r = 6)),
    axis.text.y  = element_text(size = 22, color = "black"),
    axis.line    = element_line(linewidth = 0.7, color = "black"),
    axis.ticks   = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(3.5, "pt"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 18, color = "black"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 8, 6, 6)
  )

p_pub


# ============================================================
# Table panel (top 10 terms) + alignment to plot
# ============================================================

# slimmer table (this is what makes patchwork alignment actually visible)

table_df <- gp_pathway$result %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number()) %>%
  transmute(
    Rank,
    Source = source,
    `Term ID`   = term_id,
    `Term name` = stringr::str_wrap(term_name, width = 45),
    `P (adj)`      = formatC(p_value, format = "e", digits = 2)
  )


tab_grob <- gridExtra::tableGrob(
  table_df,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 15,
    padding = unit(c(4, 4), "mm"),
    colhead = list(
      fg_params = list(fontface = "bold", col = "black"),
      bg_params = list(fill = "grey95", col = NA)
    ),
    core = list(
      fg_params = list(col = "black", hjust = 0, x = 0.02),
      bg_params = list(fill = c("grey95", "grey93"), col = NA)
    )
  )
)


g <- ggplotGrob(p_pub)

# find panel(s) robustly (works for "panel-1-1", facets, etc.)
panel_rows <- which(grepl("^panel", g$layout$name))
if (length(panel_rows) == 0) stop("Could not find a panel grob in p_pub")

# for facetted plots, take the leftmost panel column
left_panel_col <- min(g$layout$l[panel_rows])

# exact left gutter width up to the left edge of the panel
left_offset <- if (left_panel_col > 1) sum(g$widths[seq_len(left_panel_col - 1)]) else unit(0, "pt")

tab_grob_aligned <- gtable::gtable_add_cols(tab_grob, left_offset, pos = 0)

tab_pw <- patchwork::wrap_elements(full = tab_grob_aligned) +
  theme(plot.margin = margin(0, 7, 0, 0))
# ----------------------------------------------------------------------


# ============================================================
# Assemble final figure (plot + table) + export
# ============================================================

final_Path_fig <- p_pub /
  tab_pw +
  plot_layout(heights = c(2.7, 1.85))

final_Path_fig

ggsave("Pathway_Enrichments_with_table.png", final_Path_fig,
       width = 10, height = 10.5, units = "in", dpi = 600)


# ============================================================
# Selected pathway: extract member genes + export DE stats
# ============================================================

# Selected pathway genes

# sanity check: see what the overlap column is called
names(gp_pathway$result)

# pick the overlap column name robustly
overlap_col <- intersect(c("intersection", "intersections"), names(gp_pathway$result))[1]
stopifnot(!is.na(overlap_col))

genes_pathway <- gp_pathway$result %>%
  filter(term_id == "REAC:R-HSA-909733") %>%
  slice(1) %>%                              # just in case duplicates
  pull(!!sym(overlap_col)) %>%
  str_split(",") %>%
  .[[1]] %>%
  str_trim() %>%
  .[. != ""] %>%
  unique()

length(genes_pathway)
head(genes_pathway, 30)


deseq_long <- deseq_tbl %>%
  transmute(
    symbol = SYMBOL_clean,
    log2FoldChange,
    padj
  ) %>%
  separate_rows(symbol, sep = ";") %>%
  mutate(symbol = str_trim(symbol)) %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct(symbol, .keep_all = TRUE)

pathway_gene_stats <- deseq_long %>%
  filter(symbol %in% genes_pathway) %>%
  arrange(padj)

pathway_gene_stats

readr::write_csv(pathway_gene_stats, "REAC_R-HSA-909733_genes_DESeq2.csv")
