#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Extended Data Figure 4B: normalized metatranscriptomic expression
# of methane-cycling MAGs.
#
# Taxa are displayed at genus level when a genus assignment is available.
# For MAGs without an assigned genus, the family-level assignment is used.
# Bubble size represents MT_coverage_per_cell and fill represents phylum.
#
# INPUT
# 1. TSV table containing normalized expression and taxonomy for
#    methane-cycling MAGs.
# 2. Output directory.
#
# OUTPUT
# Extended Data Figure 4B as SVG.
#
# USAGE
# Rscript ED_Figure_4B.R \
#   methane_cycling_MAGs_normalized_expression.tsv \
#   output_directory

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    paste(
      "Usage: Rscript ED_Figure_4B.R",
      "<methane_cycling_MAGs_normalized_expression.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

file_path <- args[[1]]
output_dir <- args[[2]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "ED_Figure_4B_methane_cycling_MAG_expression.svg"
)

data <- readr::read_tsv(
  file_path,
  show_col_types = FALSE
)

mags_to_keep <- c(
  "Havsjomossen_1_0-25cm_2024-06_iMG_concoct_bin.141_sub",
  "Bjorsmossen_2_site3_30cm_2024-06_iMG_metabat2_bin.93",
  "Bjorsmossen_3_site3_60cm_2024-06_iMG_metabat2_bin.171",
  "Bjorsmossen_4_0-25cm_2024-06_iMG_metabat2_bin.78",
  "Bjorsmossen_9_525-540cm_2024-06_iMG_concoct_bin.28",
  "Lungsmossen_5_100-125cm_2024-06_iMG_metabat2_bin.89",
  "Lungsmossen_9_525-550cm_2024-06_iMG_concoct_bin.115",
  "Bjorsmossen_9_525-540cm_2024-06_iMG_metabat2_bin.35",
  "Havsjomossen_2_100-125cm_2024-06_iMG_concoct_bin.148",
  "Lungsmossen_11_675-700cm_2024-06_iMG_concoct_bin.2_sub",
  "Lungsmossen_9_525-550cm_2024-06_iMG_concoct_bin.47",
  "NorraRomyren_2_site3_30cm_2024-06_iMG_metabat2_bin.72",
  "NorraRomyren_3_site3_60cm_2024-06_iMG_concoct_bin.135",
  "NorraRomyren_4_0-25cm_2024-06_iMG_maxbin2_bin.002",
  "NorraRomyren_4_0-25cm_2024-06_iMG_metabat2_bin.63",
  "NorraRomyren_7_275-300cm_2024-06_iMG_metabat2_bin.119",
  "NorraRomyren_5_100-125cm_2024-06_iMG_metabat2_bin.196",
  "Lungsmossen_2_site3_30cm_2024-06_iMG_metabat2_bin.113",
  "Lungsmossen_4_0-25cm_2024-06_iMG_metabat2_bin.21",
  "Lungsmossen_4_0-25cm_2024-06_iMG_metabat2_bin.85",
  "Lungsmossen_5_100-125cm_2024-06_iMG_concoct_bin.79",
  "Havsjomossen_1_0-25cm_2024-06_iMG_concoct_bin.41",
  "Havsjomossen_1_0-25cm_2024-06_iMG_metabat2_bin.217",
  "NorraRomyren_2_site3_30cm_2024-06_iMG_metabat2_bin.229",
  "NorraRomyren_5_100-125cm_2024-06_iMG_concoct_bin.190_sub",
  "NorraRomyren_5_100-125cm_2024-06_iMG_metabat2_bin.216",
  "Lungsmossen_10_650-675cm_2024-06_iMG_concoct_bin.47_sub",
  "Lungsmossen_11_675-700cm_2024-06_iMG_concoct_bin.22",
  "Lungsmossen_11_675-700cm_2024-06_iMG_concoct_bin.87",
  "Lungsmossen_9_525-550cm_2024-06_iMG_metabat2_bin.75",
  "NorraRomyren_2_site3_30cm_2024-06_iMG_maxbin2_bin.169_sub",
  "NorraRomyren_2_site3_30cm_2024-06_iMG_metabat2_bin.135",
  "Havsjomossen_4_265-280cm_2024-06_iMG_concoct_bin.84",
  "NorraRomyren_1_site3_10cm_2024-06_iMG_concoct_bin.259_sub"
)

ordered_samples <- c(
  "Bjorsmossen_1_site3_10cm_2024-06_iMT",
  "Bjorsmossen_2_site3_30cm_2024-06_iMT",
  "Bjorsmossen_3_site3_60cm_2024-06_iMT",
  "Bjorsmossen_4_200-225cm_2024-06_iMT",
  "Bjorsmossen_5_325-350cm_2024-06_iMT",
  "Bjorsmossen_6_525-540cm_2024-06_iMT",
  "Havsjomossen_1_200-225cm_2024-06_iMT",
  "Havsjomossen_2_265-280cm_2024-06_iMT",
  "Lungsmossen_1_site3_10cm_2024-06_iMT",
  "Lungsmossen_2_site3_30cm_2024-06_iMT",
  "Lungsmossen_3_site3_60cm_2024-06_iMT",
  "Lungsmossen_4_200-225cm_2024-06_iMT",
  "Lungsmossen_5_325-350cm_2024-06_iMT",
  "Lungsmossen_6_525-550cm_2024-06_iMT",
  "Lungsmossen_7_650-675cm_2024-06_iMT",
  "NorraRomyren_1_site3_10cm_2024-06_iMT",
  "NorraRomyren_2_site3_30cm_2024-06_iMT",
  "NorraRomyren_3_site3_60cm_2024-06_iMT",
  "NorraRomyren_4_200-225cm_2024-06_iMT",
  "NorraRomyren_5_325-350cm_2024-06_iMT"
)

custom_order <- c(
  "Methanocella",
  "Methanosarcina",
  "Methanoregula",
  "MVRE01",
  "Bog-38",
  "Methanobacterium_A",
  "JACTUC01",
  "JAJRAG01",
  "JAJRAL01",
  "Methanomassiliicoccaceae",
  "Binataceae",
  "Rhodomicrobium",
  "Rhodoblastus",
  "Methylocystis",
  "Methylovirgula",
  "Methylocella"
)

clean_taxon <- function(x, prefix) {
  x <- as.character(x)
  x <- stringr::str_remove(x, paste0("^", prefix))
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "N/A", "Unclassified", "unclassified")] <- NA_character_
  x
}

long_data <- data %>%
  dplyr::mutate(
    Sample = .data$MetaT_sample,
    Abundance = suppressWarnings(
      as.numeric(.data$MT_coverage_per_cell)
    ),
    Genus_clean = clean_taxon(.data$Genus, "g__"),
    Family_clean = clean_taxon(.data$Family, "f__"),
    Phylum = clean_taxon(.data$Phylum, "p__"),
    Taxa = dplyr::coalesce(
      .data$Genus_clean,
      .data$Family_clean
    )
  ) %>%
  dplyr::filter(
    .data$genome %in% mags_to_keep,
    !is.na(.data$Taxa),
    !is.na(.data$Abundance),
    .data$Abundance > 0
  ) %>%
  dplyr::mutate(
    Sample = factor(
      .data$Sample,
      levels = ordered_samples
    ),
    Peatland = sub(
      "_.*",
      "",
      as.character(.data$Sample)
    ),
    Depth = sub(
      ".*_(\\d+(?:-\\d+)?cm)_.*",
      "\\1",
      as.character(.data$Sample)
    ),
    Depth_label = .data$Depth
  )

genus_counts_raw <- long_data %>%
  dplyr::distinct(
    .data$genome,
    .data$Taxa
  ) %>%
  dplyr::count(
    .data$Taxa,
    name = "n_genomes"
  )

plot_df <- long_data %>%
  dplyr::left_join(
    genus_counts_raw,
    by = "Taxa"
  ) %>%
  dplyr::mutate(
    Taxa_label = paste0(
      .data$Taxa,
      " (",
      .data$n_genomes,
      ")"
    )
  ) %>%
  dplyr::group_by(
    .data$Peatland,
    .data$Sample,
    .data$Depth_label,
    .data$Taxa,
    .data$Taxa_label,
    .data$Phylum
  ) %>%
  dplyr::summarise(
    Abundance = sum(
      .data$Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    .data$Peatland,
    .data$Sample,
    .data$Depth_label,
    .data$Taxa,
    .data$Taxa_label
  ) %>%
  dplyr::slice_max(
    order_by = .data$Abundance,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Taxa = factor(
      .data$Taxa,
      levels = custom_order
    )
  )

label_levels <- plot_df %>%
  dplyr::distinct(
    .data$Taxa,
    .data$Taxa_label
  ) %>%
  dplyr::arrange(
    .data$Taxa
  ) %>%
  dplyr::pull(
    .data$Taxa_label
  )

plot_df <- plot_df %>%
  dplyr::mutate(
    Taxa_label = factor(
      .data$Taxa_label,
      levels = rev(label_levels)
    )
  )

custom_colors <- c(
  "Desulfobacterota_B" = "#FF66CC",
  "Halobacteriota" = "#2171B5",
  "Methanobacteriota" = "#C6DBEF",
  "Pseudomonadota" = "#8C510A",
  "Thermoplasmatota" = "#3690C0",
  "Thermoproteota" = "#08306B"
)

expr_breaks <- c(
  0.0001,
  0.001,
  0.01,
  0.1,
  1,
  3.5
)

expr_labels <- c(
  "0.0001–0.001",
  "0.001–0.01",
  "0.01–0.1",
  "0.1–1",
  "1–3.5",
  ">3.5"
)

max_cap <- max(expr_breaks)
min_cap <- min(expr_breaks)

plot_df <- plot_df %>%
  dplyr::mutate(
    Abundance_plot = pmin(
      pmax(
        .data$Abundance,
        min_cap
      ),
      max_cap
    )
  )

bubble_plot <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = .data$Sample,
    y = .data$Taxa_label,
    size = .data$Abundance_plot,
    fill = .data$Phylum
  )
) +
  ggplot2::geom_point(
    shape = 21,
    color = "grey20",
    stroke = 0.25,
    alpha = 0.95
  ) +
  ggplot2::scale_fill_manual(
    values = custom_colors,
    drop = FALSE
  ) +
  ggplot2::scale_x_discrete(
    labels = function(x) {
      depth_map <- unique(
        plot_df[, c("Sample", "Depth_label")]
      )

      stats::setNames(
        depth_map$Depth_label,
        depth_map$Sample
      )[x]
    }
  ) +
  ggplot2::scale_size_continuous(
    breaks = expr_breaks,
    labels = expr_labels,
    range = c(1.5, 14),
    name = "Expression level"
  ) +
  ggplot2::labs(
    x = "Depth",
    y = "Genus or family (n MAGs)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      size = 14
    ),
    axis.text.x = ggplot2::element_text(
      size = 9,
      angle = 45,
      hjust = 1
    ),
    legend.text = ggplot2::element_text(
      size = 12
    ),
    legend.title = ggplot2::element_text(
      size = 14
    ),
    legend.position = "right"
  ) +
  ggplot2::facet_wrap(
    ~ Peatland,
    scales = "free_x",
    ncol = 2
  )

print(bubble_plot)

ggplot2::ggsave(
  filename = out_svg,
  plot = bubble_plot,
  width = 13,
  height = 9,
  units = "in",
  device = svglite::svglite,
  bg = "white",
  limitsize = FALSE
)
