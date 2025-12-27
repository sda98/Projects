# ============================================================
# 06_go_enrichment_gprofiler.R
# Purpose:
#   - Perform GO enrichment (BP/MF/CC) using g:Profiler (gprofiler2)
#   - Use significant DE genes as query + all tested genes as background
#   - Export enrichment table
#   - Generate publication-style GO dot plot + top-10 table panel
#
# Inputs (expected to already exist in environment):
#   - deseq_tbl (data.frame with padj, log2FoldChange, SYMBOL_clean, etc.)
#
# Output files:
#   - gprofiler_GO.csv
#   - GO_enrichment_with_table.png
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
# Build query and background gene lists
# ============================================================

## Gene Ontology Analysis

sig_genes <- deseq_tbl %>%
  filter(!is.na(padj), padj <= 0.01, abs(log2FoldChange) >= 1) %>%
  transmute(symbol = SYMBOL_clean) %>%
  separate_rows(symbol, sep = ";") %>%
  mutate(symbol = str_trim(symbol)) %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct() %>%
  pull(symbol)

length(sig_genes)
head(sig_genes, 20)

bg_genes <- deseq_tbl %>%
  filter(!is.na(padj)) %>%
  transmute(symbol = SYMBOL_clean) %>%
  separate_rows(symbol, sep = ";") %>%
  mutate(symbol = str_trim(symbol)) %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct() %>%
  pull(symbol)


# ============================================================
# Run g:Profiler (GO:BP/GO:MF/GO:CC)
# ============================================================

gp_GO <- gost(
  query = sig_genes,
  organism = "hsapiens",
  correction_method = "g_SCS",
  user_threshold = 0.001,
  sources = c("GO:BP","GO:MF","GO:CC"),
  custom_bg = bg_genes,
  evcodes = TRUE,
)

head(gp_GO$result)
readr::write_csv(gp_GO$result, "gprofiler_GO.csv")


# ============================================================
# g:Profiler plots (interactive + static)
# ============================================================

gostplot(gp_GO, capped = FALSE, interactive = TRUE) # Interactive plot


GO <- gostplot(gp_GO, capped = FALSE, interactive = FALSE)


# ============================================================
# Top terms table (top 10 by p-value)
# ============================================================

top_table <- gp_GO$result %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank, source, term_id, term_name, p_value) %>%
  transmute(
    Rank,
    Source = source,
    `Term ID` = term_id,
    `Term name` = str_trunc(term_name, 70),
    `Adj. p` = formatC(p_value, format = "e", digits = 2)
  )

top_terms <- top_table$`Term ID`


# ============================================================
# Customize GO plot styling (outlined points + highlight top 10)
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

GO$layers[[pt_idx]]$mapping <- m

# outlined points with constant size
GO$layers[[pt_idx]]$aes_params$shape  <- 21
GO$layers[[pt_idx]]$aes_params$colour <- "black"
GO$layers[[pt_idx]]$aes_params$stroke <- 0.5
GO$layers[[pt_idx]]$aes_params$size   <- 5.5   # <- pick your constant point size

# palette for GO sources
my_pal <- c("GO:BP"= "#ff9900", "GO:MF" = "#dc3912", "GO:CC" = "#109618")
GO <- GO + scale_fill_manual(values = my_pal) + guides(size = "none")
# Base alpha for everything
base_alpha <- 0.5

# Add a rank column to GO$data for top 10
GO$data <- GO$data %>%
  left_join(
    gp_GO$result %>%
      arrange(p_value) %>%
      slice_head(n = 10) %>%
      mutate(Rank = row_number()) %>%
      dplyr::select(term_id, Rank),
    by = "term_id"
  ) %>%
  mutate(alpha_pt = ifelse(!is.na(Rank), 1, base_alpha))

# Map alpha per point (and turn off alpha legend)
GO <- GO +
  aes(alpha = alpha_pt) +
  scale_alpha_identity(guide = "none")

# Add rank labels for the top 10
label_df <- GO$data %>% filter(!is.na(Rank))

# remove the fixed alpha so mapping can work
GO$layers[[pt_idx]]$aes_params$alpha <- NULL

# map alpha to alpha_pt on the point layer
GO$layers[[pt_idx]]$mapping$alpha <- rlang::sym("alpha_pt")

# use the alpha values as-is (no legend)
GO <- GO + scale_alpha_identity(guide = "none")

label_df <- GO$data %>% dplyr::filter(!is.na(Rank))


# ============================================================
# Publication-style GO plot (rank labels + theme)
# ============================================================

p_pub <- GO +
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
    box.padding = 1.1,
    point.padding = 0.82,
    segment.color = "black",
    segment.size = 0.8,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 3
  ) +
  labs(title = "GO enrichment", y = expression(-log[10](p[adj]))) +
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

table_df <- gp_GO$result %>%
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

final_GO_fig <- p_pub /
  plot_spacer() /
  tab_pw +
  plot_layout(heights = c(2.7, 0, 1.85))

final_GO_fig <- p_pub /
  tab_pw +
  plot_layout(heights = c(2.7, 1.85))

final_GO_fig


ggsave("GO_enrichment_with_table.png", final_GO_fig,
       width = 9.25, height = 10.5, units = "in", dpi = 600)
