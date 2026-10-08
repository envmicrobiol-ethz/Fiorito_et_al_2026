#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Extended Data Figure 4C: a heatmap of average
# MT_coverage_per_cell values for selected methanogenesis and associated
# metabolic markers across target methanogenic genera.
#
# Expression values are averaged for each genus and KEGG orthologue.
# Positive values are divided into Low, Medium and High classes using the
# 25th and 75th percentiles calculated across all positive values.
#
# INPUT
# 1. MAG-level KEGG expression-profile TSV containing:
#      Genus, KEGG_ko, definition and MT_coverage_per_cell
# 2. Output directory.
#
# OUTPUT
# Extended Data Figure 4C as SVG and PDF.
#
# USAGE
# Rscript ED_Figure_4C.R \
#   selected_methanogen_MAGs_KEGG_expression.tsv \
#   output_directory

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(colorspace)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    paste(
      "Usage: Rscript ED_Figure_4C.R",
      "<KEGG_expression.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

infile <- args[[1]]
output_dir <- args[[2]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "ED_Figure_4C_methanogen_metabolic_expression_heatmap.svg"
)

out_pdf <- file.path(
  output_dir,
  "ED_Figure_4C_methanogen_metabolic_expression_heatmap.pdf"
)

df <- readr::read_tsv(
  infile,
  show_col_types = FALSE
)

###############################
### 3. MARKER SETS
###############################
marker_sets <- list(
  "Lactate" = c("K00016","K21836","K00782","K00104"),
  "Nitrogen" = c("K02588","K02586","K02591"),
  "mcrABCDG" = c("K00399","K00401","K03421","K00402","K03422"),
  "cdhDE" = c("K00194","K00197"),
  "ack, pta, acs" = c("K00925","K00625","K01895"),
  "mtaABC" = c("K14081","K14080","K04480"),
  "mtbABC" = c("K16178","K16179"),
  "mttBC" = c("K14083","K14084"),
  "mtmBC" = c("K16176"),
  "mvh" = c("K14126","K14127","K14128"),
  "frh" = c("K00440","K00443","K00441"),
  "fmdABCDE" = c("K00200","K00201","K00202","K00203","K11261"),
  "hdrABC" = c("K03388","K03389","K03390","K08264","K08265")
)

all_markers <- unique(unlist(marker_sets))

###############################
### 4. CATEGORY COLORS
###############################
cat_base_cols <- c(
  "mcrABCDG"      = "#F26D5B",
  "hdrABC"        = "#D94B45",
  "fmdABCDE"      = "#B83A3A",
  "frh"           = "#D3D3D3",
  "mvh"           = "darkgrey",
  "ack, pta, acs" = "darkorange",
  "cdhDE"         = "darkred",
  "mtmBC"         = "#E4B3C8",
  "mtbABC"        = "#D88FB0",
  "mtaABC"        = "#E56AA6",
  "mttBC"         = "#E8A8D8",
  "Nitrogen"      = "#8BCF7A",
  "Lactate"       = "#A85A1A"
)

###############################
### 5. DESIRED TOP-DOWN ORDER
###############################
desired_order <- c(
  "cdhDE",
  "fmdABCDE",
  "hdrABC",
  "mcrABCDG",
  "frh",
  "mvh",
  "ack, pta, acs",
  "Lactate",
  "mtmBC",
  "mtbABC",
  "mtaABC",
  "mttBC",
  "Nitrogen"
)

###############################
### 6. TARGET GENERA
###############################
target_genera <- c(
  "Bog-38","Methanocella","Methanoregula","Methanosarcina",
  "MVRE01","Methanomassiliicoccaceae","Methanobacterium_A",
  "JACTUC01","JAJRAG01","JAJRAL01"
)

###############################
### 7. CLEAN INPUT DF
###############################
expr_col <- "MT_coverage_per_cell"

df_clean <- df %>%
  mutate(
    Genus   = str_remove(Genus, "^g__") |> str_trim(),
    KEGG_ko = toupper(str_trim(KEGG_ko)),
    expr    = .data[[expr_col]]
  )

###############################
### 8. KO DEFINITIONS
###############################
kegg_def <- df_clean %>%
  filter(KEGG_ko %in% all_markers) %>%
  group_by(KEGG_ko) %>%
  summarise(definition = first(na.omit(definition)), .groups = "drop")

###############################
### 9. MEAN EXPRESSION PER GENUS × KO
###############################
df_sum <- df_clean %>%
  filter(
    KEGG_ko %in% all_markers,
    Genus %in% target_genera,
    !is.na(expr)
  ) %>%
  group_by(Genus, KEGG_ko) %>%
  summarise(mean_expr = mean(expr), .groups = "drop")

###############################
### 10. FULL GRID
###############################
df_sum_full <- expand_grid(
  Genus   = target_genera,
  KEGG_ko = all_markers
) %>%
  left_join(df_sum, by = c("Genus", "KEGG_ko")) %>%
  left_join(kegg_def, by = "KEGG_ko") %>%
  mutate(
    mean_expr  = replace_na(mean_expr, 0),
    definition = if_else(is.na(definition), KEGG_ko, definition)
  )

###############################
### 11. ADD CATEGORY
###############################
ko_to_cat <- enframe(marker_sets) %>%
  unnest_longer(value) %>%
  rename(Category = name, KEGG_ko = value) %>%
  mutate(KEGG_ko = toupper(KEGG_ko))

df_sum2 <- df_sum_full %>%
  left_join(ko_to_cat, by = "KEGG_ko") %>%
  filter(!is.na(Category))

###############################
### 12. BUILD HEAT_DF
###############################
heat_df <- df_sum2 %>%
  mutate(
    Gene     = paste0(KEGG_ko, " — ", definition),
    Category = factor(Category, levels = desired_order)
  ) %>%
  arrange(Category, KEGG_ko) %>%
  select(Gene, Category, Genus, mean_expr) %>%
  pivot_wider(names_from = Genus, values_from = mean_expr, values_fill = 0)

# ordine manuale dei KO dentro ogni categoria
ko_order_within_cat <- c(
  "K00399", "K00401", "K03421", "K03422", "K00402",
  "K03388", "K03389", "K03390", "K08264", "K08265",
  "K00200", "K00201", "K00202", "K00203", "K11261",
  "K00440", "K00441", "K00443",
  "K14126", "K14127", "K14128",
  "K00925", "K00625", "K01895",
  "K00194", "K00197",
  "K16176",
  "K16178", "K16179",
  "K14081", "K14080", "K04480",
  "K14083", "K14084",
  "K02588", "K02586", "K02591",
  "K00016", "K21836", "K00782", "K00104"
)

row_order <- heat_df %>%
  mutate(
    KEGG_ko = str_extract(Gene, "^K\\d+"),
    ko_rank = match(KEGG_ko, ko_order_within_cat)
  ) %>%
  arrange(Category, ko_rank, KEGG_ko) %>%
  pull(Gene)

heat_df <- heat_df %>%
  slice(match(row_order, Gene)) %>%
  filter(!is.na(Category)) %>%
  mutate(Category = droplevels(Category))

###############################
### 13. FIX CATEGORIES
###############################
present_cats <- unique(heat_df$Category) |> as.character()
cat_levels   <- desired_order[desired_order %in% present_cats]


cat_cols <- cat_base_cols[cat_levels]

###############################
### 14. CATEGORY SEPARATORS
###############################
lens_by_cat <- heat_df %>%
  count(Category) %>%
  arrange(factor(Category, levels = cat_levels)) %>%
  pull(n)

starts <- c(1, head(cumsum(lens_by_cat) + 1, -1))
ends   <- cumsum(lens_by_cat)
mids   <- (starts + ends) / 2

n_rows    <- nrow(heat_df)
gap_y     <- cumsum(lens_by_cat)[-length(lens_by_cat)] + 0.5
gap_y_rev <- (n_rows + 1) - gap_y

###############################
### 15. DISCRETE LEVELS
### 0 / Low / Medium / High
###############################
positive_vals <- df_sum2$mean_expr[df_sum2$mean_expr > 0]


q25 <- quantile(positive_vals, 0.25, na.rm = TRUE)
q75 <- quantile(positive_vals, 0.75, na.rm = TRUE)


###############################
### 16. COLORS FOR LOW/MEDIUM/HIGH
###############################
hex_to_hcl <- function(hex) {
  rgb <- col2rgb(hex)
  luv <- convertColor(t(rgb), from = "sRGB", to = "Luv", scale.in = 255)
  H <- (atan2(luv[,3], luv[,2]) * 180 / pi) %% 360
  C <- sqrt(luv[,2]^2 + luv[,3]^2)
  L <- luv[,1]
  out <- c(H = H, C = C, L = L)
  names(out) <- c("H","C","L")
  out
}

hcl_mat <- do.call(
  rbind,
  lapply(cat_levels, function(cat) hex_to_hcl(cat_cols[[cat]]))
)
rownames(hcl_mat) <- cat_levels

cat_map <- tibble(
  Category = factor(cat_levels, levels = cat_levels),
  hex_base = unname(cat_cols),
  H0 = hcl_mat[, "H"],
  C0 = hcl_mat[, "C"],
  L0 = hcl_mat[, "L"]
) %>%
  mutate(
    neutral  = C0 < 5,
    col_zero = "#F6F1E8",
    col_low = if_else(
      neutral,
      "#E0E0E0",
      lighten(hex_base, 0.45)
    ),
    col_medium = if_else(
      neutral,
      "#B5B5B5",
      lighten(hex_base, 0.20)
    ),
    col_high = if_else(
      neutral,
      "#6F6F6F",
      hex_base
    )
  )

###############################
### 17. LONG DF WITH BINNING
###############################
genus_cols <- intersect(target_genera, colnames(heat_df))
plot_gene_levels <- rev(heat_df$Gene)

heat_mat <- heat_df %>%
  select(all_of(genus_cols)) %>%
  as.matrix()

rownames(heat_mat) <- heat_df$Gene

df_long <- as.data.frame(heat_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Genus", values_to = "mean_expr") %>%
  left_join(heat_df %>% select(Gene, Category), by = "Gene") %>%
  mutate(
    Gene     = factor(Gene, levels = plot_gene_levels),
    Genus    = factor(Genus, levels = genus_cols),
    Expr_bin = case_when(
      mean_expr == 0 ~ "0",
      mean_expr <= q25 ~ "Low",
      mean_expr <= q75 ~ "Medium",
      TRUE ~ "High"
    ),
    Expr_bin = factor(Expr_bin, levels = c("0","Low","Medium","High"))
  ) %>%
  left_join(cat_map, by = "Category") %>%
  mutate(
    fill_col = case_when(
      Expr_bin == "0"      ~ col_zero,
      Expr_bin == "Low"    ~ col_low,
      Expr_bin == "Medium" ~ col_medium,
      Expr_bin == "High"   ~ col_high
    )
  )

###############################
### 18. LEFT CATEGORY STRIP
###############################
strip_df <- heat_df %>%
  select(Gene, Category) %>%
  mutate(
    Gene = factor(Gene, levels = plot_gene_levels),
    x = "Category"
  )

p_strip <- ggplot(strip_df, aes(x = x, y = Gene, fill = Category)) +
  geom_tile(color = "black", linewidth = 0.2, width = 0.9) +
  scale_fill_manual(values = cat_cols, guide = "none") +
  scale_y_discrete(limits = plot_gene_levels) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_text(size = 11),
    axis.title.x = element_blank()
  ) +
  geom_hline(yintercept = gap_y_rev, color = "grey30", linewidth = 0.4)

###############################
### 19. RIGHT CATEGORY BLOCKS
###############################
right_blocks <- tibble(
  Category = factor(cat_levels, levels = cat_levels),
  ymin = starts - 0.5,
  ymax = ends + 0.5,
  ymid = mids
) %>%
  left_join(cat_map, by = "Category") %>%
  mutate(
    label    = str_wrap(as.character(Category), width = 18),
    text_col = ifelse(L0 < 60, "white", "black")
  )

p_catcol <- ggplot() +
  geom_rect(
    data = right_blocks,
    aes(xmin = 0.55, xmax = 1.05, ymin = ymin, ymax = ymax, fill = Category),
    color = "black", linewidth = 0.2
  ) +
  geom_text(
    data = right_blocks,
    aes(x = 1.15, y = ymid, label = label, color = text_col),
    hjust = 0, vjust = 0.5, size = 4.5
  ) +
  scale_fill_manual(values = cat_cols, guide = "none") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.5, 3), expand = c(0, 0)) +
  scale_y_reverse(limits = c(n_rows + 0.5, 0.5), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text  = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

###############################
### 20. MAIN HEATMAP
###############################
p_main <- ggplot(df_long, aes(x = Genus, y = Gene, fill = fill_col)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_identity() +
  scale_x_discrete(limits = genus_cols) +
  scale_y_discrete(limits = plot_gene_levels, position = "right") +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    axis.title  = element_blank(),
    panel.grid  = element_blank()
  ) +
  geom_hline(yintercept = gap_y_rev, color = "grey30", linewidth = 0.4)

###############################
### 21. INTENSITY LEGEND
###############################
legend_df <- tibble(
  level = factor(c("0", "Low", "Medium", "High"), levels = c("0","Low","Medium","High")),
  x = c(1, 2, 3, 4),
  y = 1,
  fill = c("#F6F1E8", "#D8D8D8", "#A6A6A6", "#5C5C5C")
)

p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = fill)) +
  geom_tile(color = "black", linewidth = 0.2, width = 0.9, height = 0.9) +
  geom_text(aes(label = level), y = 1.7, size = 5) +
  scale_fill_identity() +
  coord_fixed(clip = "off") +
  theme_void() +
  xlim(0.4, 4.6) +
  ylim(0.4, 2.2)

###############################
### 22. FINAL ASSEMBLY
###############################
main_panel <- p_strip + p_main + p_catcol + plot_layout(widths = c(1.8, 10.5, 3.4))
final_plot <- main_panel / p_legend + plot_layout(heights = c(20, 1.8))

print(final_plot)

###############################
### 23. SAVE
###############################

ggplot2::ggsave(
  filename = out_svg,
  plot = final_plot,
  width = 11.5,
  height = 18,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

if (capabilities("cairo")) {
  ggplot2::ggsave(
    filename = out_pdf,
    plot = final_plot,
    device = grDevices::cairo_pdf,
    width = 11.5,
    height = 18,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}
