# ============================================================
# 07_pathway_enrichment_gprofiler.R
# Purpose:
#   - Pathway enrichment using g:Profiler (KEGG, Reactome, WikiPathways)
#   - Export enrichment results table
#   - Generate publication-style pathway dot plot + top-10 table panel
#   - Extract genes for a selected pathway and export DE stats
#
# Inputs (expected to already exist in environment):
#   - sig_genes (character vector of significant gene symbols, created in 06_go_enrichment_gprofiler.R)
#   - bg_genes  (character vector of background gene symbols, created in 06_go_enrichment_gprofiler.R)
#   - deseq_tbl (Annotated DESeq2 output table created in 04_differential_expression_deseq2.R)
#
# Output files:
#   - gprofiler_pathways.csv (Pathways enrichment table, g:SCS-based Padj < 0.001)
#   - Pathways_enrichments_with_table.png (Pathways dot plot with top-10 GO terms table panel)
#   - REAC_R-HSA-909733_genes_DESeq2.csv (table of genes enriched in the pathway of choice, can be changed)
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

gp_pathway <- gost(
  query = sig_genes,
  organism = "hsapiens",
  correction_method = "g_SCS",
  user_threshold = 0.001,             # Can customize the Padj cutoff for Pathways terms if needed
  sources = c("KEGG", "REAC", "WP"),
  custom_bg = bg_genes,
  evcodes = TRUE,
)

head(gp_pathway$result)
readr::write_csv(gp_pathway$result, "gprofiler_pathways.csv")


# ============================================================
# g:Profiler plots (interactive + static)
# ============================================================

gostplot(gp_pathway, capped = FALSE, interactive = TRUE) # Interactive plot


Pathway <- gostplot(gp_pathway, capped = FALSE, interactive = FALSE)


# ============================================================
# Customize Pathways plot styling (outlined points + highlight top 10)
# ============================================================

## Find the index of the first point layer (the Pathways term dots) in the ggplot object to edit the point layer

pt_idx <- which(vapply(Pathway$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]

## Accessing the mappings of this layer

m <- Pathway$layers[[pt_idx]]$mapping

## Keeping original colors as fills (if color was mapped by gProfiler)

if (is.null(m$fill)) {
  if (!is.null(m$colour)) m$fill <- m$colour
  if (!is.null(m$color))  m$fill <- m$color
}

## Removing colour mapping to edit them later

m$colour <- NULL
m$color  <- NULL

## Removing size mapping to later assign the constant size

m$size <- NULL

## Reassigning the mapping to m object

Pathway$layers[[pt_idx]]$mapping <- m

## Editing the point layer aesthetics

Pathway$layers[[pt_idx]]$aes_params$shape  <- 21
Pathway$layers[[pt_idx]]$aes_params$colour <- "black"
Pathway$layers[[pt_idx]]$aes_params$stroke <- 0.5
Pathway$layers[[pt_idx]]$aes_params$size   <- 5.5   

## Palette for point aesthetics

my_pal <- c("KEGG" = "#dd4477", "REAC" = "#3366cc", "WP" = "#0099c6")
Pathway <- Pathway + scale_fill_manual(values = my_pal) + guides(size = "none")

## Base alpha for everything

base_alpha <- 0.5

## Adding a rank column to Pathway$data for top 10, and add desired alpha values in the GO$data

Pathway$data <- Pathway$data %>%
  left_join(
    gp_pathway$result %>%
      arrange(p_value) %>%
      slice_head(n = 10) %>%
      mutate(Rank = row_number()) %>%
      dplyr::select(term_id, Rank),
    by = "term_id"
  ) %>%
  mutate(alpha_pt = ifelse(!is.na(Rank), 1, base_alpha))


## Removing the fixed alpha so mapping can work

Pathway$layers[[pt_idx]]$aes_params$alpha <- NULL

## Mapping alpha to alpha_pt column on the point layer

Pathway$layers[[pt_idx]]$mapping$alpha <- rlang::sym("alpha_pt")

## Using the alpha values as-is without scaling and without legend

Pathway <- Pathway + scale_alpha_identity(guide = "none")

label_df <- Pathway$data %>% dplyr::filter(!is.na(Rank))

## Removing colour aesthetics from Pathway

Pathway$scales$scales <- Filter(
  function(s) !any(s$aesthetics %in% c("colour", "color")),
  Pathway$scales$scales
)

# ============================================================
# Publication-style Pathway plot (rank labels + theme)
# ============================================================

p_pub_pathway <- Pathway +
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
  labs(title = "Pathway enrichment", y = expression(-log[10](p[adj]))) +
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

p_pub_pathway


# ============================================================
# Table panel (top 10 terms) + Table as a plot
# ============================================================

table_df_path <- gp_pathway$result %>%
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


tab_grob_path <- gridExtra::tableGrob(
  table_df_path,
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

final_path_fig <- p_pub_pathway /
  tab_grob_path +
  plot_layout(heights = c(2.7, 1.85))

final_path_fig

ggsave("Pathways_enrichments_with_table.png", final_path_fig,
       width = 9.25, height = 10.5, units = "in", dpi = 600)


# ============================================================
# Selected pathway: extract member genes + export DE stats
# ============================================================

## Sanity check: see what the overlap column is called

names(gp_pathway$result)

## Picking the overlap column names to get gene names

overlap_col <- intersect(c("intersection", "intersections"), names(gp_pathway$result))[1]
stopifnot(!is.na(overlap_col))

## Getting the gene names from the term of interest

genes_pathway <- gp_pathway$result %>%
  filter(term_id == "REAC:R-HSA-909733") %>% # Can repace with another term of intest
  slice(1) %>%                              
  pull(!!sym(overlap_col)) %>%
  str_split(",") %>%
  .[[1]] %>%
  str_trim() %>%
  .[. != ""] %>%
  unique()

length(genes_pathway)
head(genes_pathway, 30)


## Extracting DESeq2 Log2FC and Padj information for the selected genes 

### Getting clean Log2FC/Padj table

deseq_clean <- deseq_tbl %>%
  transmute(
    symbol = SYMBOL_clean,
    log2FoldChange,
    padj
  ) %>%
  separate_rows(symbol, sep = ";") %>%
  mutate(symbol = str_trim(symbol)) %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct(symbol, .keep_all = TRUE)

### Filtering Log2FC/Padj table by the genes of interest

pathway_gene_stats <- deseq_clean %>%
  filter(symbol %in% genes_pathway) %>%
  arrange(padj)

pathway_gene_stats

## Saving the final result

readr::write_csv(pathway_gene_stats, "REAC_R-HSA-909733_genes_DESeq2.csv") # Can rename CV file if needed
