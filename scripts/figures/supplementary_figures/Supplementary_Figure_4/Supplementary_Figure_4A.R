#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Extended Data Figure 4A: abundance of methane-cycling MAGs
# across metagenomic samples.
#
# Taxa are displayed at genus level when a genus assignment is available.
# For MAGs without an assigned genus, the family-level assignment is used.
# Bubble size represents abundance and fill represents phylum.
#
# INPUT
# 1. TSV table containing methane-cycling MAG taxonomy and abundance.
# 2. Output directory.
#
# OUTPUT
# Extended Data Figure 4A as SVG.
#
# USAGE
# Rscript ED_Figure_4A.R \
#   methane_cycling_MAGs_taxonomy_abundance.tsv \
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
      "Usage: Rscript ED_Figure_4A.R",
      "<methane_cycling_MAGs_taxonomy_abundance.tsv>",
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
  "ED_Figure_4A_methane_cycling_MAG_abundance.svg"
)

data <- readr::read_tsv(
  file_path,
  show_col_types = FALSE
)

ordered_samples <- c(
  "Bjorsmossen_1_site3_10cm_2024-06_iMG",
  "Bjorsmossen_2_site3_30cm_2024-06_iMG",
  "Bjorsmossen_3_site3_60cm_2024-06_iMG",
  "Bjorsmossen_4_0-25cm_2024-06_iMG",
  "Bjorsmossen_5_100-125cm_2024-06_iMG",
  "Bjorsmossen_6_200-225cm_2024-06_iMG",
  "Bjorsmossen_7_275-300cm_2024-06_iMG",
  "Bjorsmossen_8_325-350cm_2024-06_iMG",
  "Bjorsmossen_9_525-540cm_2024-06_iMG",
  "Havsjomossen_1_0-25cm_2024-06_iMG",
  "Havsjomossen_2_100-125cm_2024-06_iMG",
  "Havsjomossen_3_200-225cm_2024-06_iMG",
  "Havsjomossen_4_265-280cm_2024-06_iMG",
  "Lungsmossen_1_site3_10cm_2024-06_iMG",
  "Lungsmossen_2_site3_30cm_2024-06_iMG",
  "Lungsmossen_3_site3_60cm_2024-06_iMG",
  "Lungsmossen_4_0-25cm_2024-06_iMG",
  "Lungsmossen_5_100-125cm_2024-06_iMG",
  "Lungsmossen_6_200-225cm_2024-06_iMG",
  "Lungsmossen_7_275-300cm_2024-06_iMG",
  "Lungsmossen_8_325-350cm_2024-06_iMG",
  "Lungsmossen_9_525-550cm_2024-06_iMG",
  "Lungsmossen_10_650-675cm_2024-06_iMG",
  "Lungsmossen_11_675-700cm_2024-06_iMG",
  "NorraRomyren_1_site3_10cm_2024-06_iMG",
  "NorraRomyren_2_site3_30cm_2024-06_iMG",
  "NorraRomyren_3_site3_60cm_2024-06_iMG",
  "NorraRomyren_4_0-25cm_2024-06_iMG",
  "NorraRomyren_5_100-125cm_2024-06_iMG",
  "NorraRomyren_6_200-225cm_2024-06_iMG",
  "NorraRomyren_7_275-300cm_2024-06_iMG",
  "NorraRomyren_8_325-350cm_2024-06_iMG"
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

data_clean <- data %>%
  dplyr::mutate(
    Genus_clean = clean_taxon(.data$Genus, "g__"),
    Family_clean = clean_taxon(.data$Family, "f__"),
    Phylum = clean_taxon(.data$Phylum, "p__"),
    Taxa = dplyr::coalesce(
      .data$Genus_clean,
      .data$Family_clean
    )
  )

genus_counts_raw <- data_clean %>%
  dplyr::filter(
    !is.na(.data$Taxa)
  ) %>%
  dplyr::distinct(
    .data$genome,
    .data$Taxa
  ) %>%
  dplyr::count(
    .data$Taxa,
    name = "n_genomes"
  )

long_data <- data_clean %>%
  tidyr::pivot_longer(
    cols = 9:ncol(.),
    names_to = "Sample",
    values_to = "Abundance"
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
      ".*_(\\d+(?:-\\d+)?cm).*",
      "\\1",
      as.character(.data$Sample)
    ),
    Depth_label = .data$Depth,
    Abundance = suppressWarnings(
      as.numeric(.data$Abundance)
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$Taxa),
    !is.na(.data$Abundance),
    .data$Abundance > 0
  ) %>%
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
  )

plot_df <- long_data %>%
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

abundance_breaks <- c(
  0.001,
  0.01,
  0.1,
  1.5
)

abundance_labels <- c(
  "0.001–0.01",
  "0.01–0.1",
  "0.1–1.5",
  ">1.5"
)

bubble_plot <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = .data$Sample,
    y = .data$Taxa_label,
    size = .data$Abundance,
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
    values = custom_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = function(x) {
      depth_map <- unique(
        long_data[, c("Sample", "Depth_label")]
      )

      stats::setNames(
        depth_map$Depth_label,
        depth_map$Sample
      )[x]
    }
  ) +
  ggplot2::scale_size_continuous(
    breaks = abundance_breaks,
    labels = abundance_labels,
    range = c(1.5, 14),
    name = "Abundance (%)"
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
