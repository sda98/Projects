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

## Significant genes

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

## Background genes

bg_genes <- deseq_tbl %>%
  filter(!is.na(padj)) %>%
  transmute(symbol = SYMBOL_clean) %>%
  separate_rows(symbol, sep = ";") %>%
  mutate(symbol = str_trim(symbol)) %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct() %>%
  pull(symbol)

length(bg_genes)
head(bg_genes, 20)


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
# Customize GO plot styling (outlined points + highlight top 10)
# ============================================================

## Find the index of the first point layer (the GO term dots) in the ggplot object to edit the point layer

pt_idx <- which(vapply(GO$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]

## Accessing the mappings of this layer

m <- GO$layers[[pt_idx]]$mapping

## Keeping original colors as fills (if color was mapped by gProfiler)

if (is.null(m$fill)) {
  if (!is.null(m$colour)) m$fill <- m$colour
  if (!is.null(m$color))  m$fill <- m$color
}

## Removing colour mapping to edit them later

m$colour <- NULL
m$color  <- NULL

## Removjng size mapping to later assign the constant size

m$size <- NULL

## Reassigning the mapping to m object

GO$layers[[pt_idx]]$mapping <- m

## Editing the point layer aesthetics

GO$layers[[pt_idx]]$aes_params$shape  <- 21
GO$layers[[pt_idx]]$aes_params$colour <- "black"
GO$layers[[pt_idx]]$aes_params$stroke <- 0.5
GO$layers[[pt_idx]]$aes_params$size   <- 5.5   

## Palette for point aesthetics

my_pal <- c("GO:BP"= "#ff9900", "GO:MF" = "#dc3912", "GO:CC" = "#109618")
GO <- GO + scale_fill_manual(values = my_pal) + guides(size = "none")

## Base alpha for everything

base_alpha <- 0.5

## Adding a rank column to GO$data for top 10, and add desired alpha values in the GO$data

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


## Removing the fixed alpha so mapping can work

GO$layers[[pt_idx]]$aes_params$alpha <- NULL

## Mapping alpha to alpha_pt column on the point layer

GO$layers[[pt_idx]]$mapping$alpha <- rlang::sym("alpha_pt")

## Using the alpha values as-is without scaling and without legend

GO <- GO + scale_alpha_identity(guide = "none")

label_df <- GO$data %>% dplyr::filter(!is.na(Rank))

## Removing colour aesthetics from GO

GO$scales$scales <- Filter(
  function(s) !any(s$aesthetics %in% c("colour", "color")),
  GO$scales$scales
)

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
    fill = "grey95",               
    label.size = 0.6,             
    label.r = unit(5, "pt"),       
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
  theme_classic() +
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
# Table panel (top 10 terms) + Table as a plot
# ============================================================

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

# ============================================================
# Assembling final figure (plot + table) 
# ============================================================

final_GO_fig <- p_pub /
  tab_grob +
  plot_layout(heights = c(2.7, 1.85))

final_GO_fig

ggsave("GO_enrichment_with_table.png", final_GO_fig,
       width = 9.25, height = 10.5, units = "in", dpi = 600)
