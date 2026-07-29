#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Extended Data Figure 5:
#   A. Donut chart showing the 15 most frequent PFAM annotations among
#      Tier1 MHC proteins.
#   B. Mini pie charts showing PSORT localization within each of the
#      15 most frequent PFAM annotations.
#
# The first PFAM annotation associated with each protein is retained.
#
# INPUT
# 1. Excel table containing Tier1 MHC proteins, PSORT localization and
#    PFAM annotation columns named pfam_1, pfam_2, ...
# 2. Output directory.
#
# OUTPUT
# Extended Data Figure 5 as SVG, PDF and PNG.
#
# USAGE
# Rscript ED_Figure_5.R \
#   MHC_overview_Tier1_with_PFAM_and_PSORT.xlsx \
#   output_directory

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
  library(ggh4x)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    paste(
      "Usage: Rscript ED_Figure_5.R",
      "<Tier1_MHC_overview.xlsx>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

file_overview <- args[[1]]
output_dir <- args[[2]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "ED_Figure_5_Tier1_PFAM_annotations_and_PSORT.svg"
)

out_pdf <- file.path(
  output_dir,
  "ED_Figure_5_Tier1_PFAM_annotations_and_PSORT.pdf"
)

out_png <- file.path(
  output_dir,
  "ED_Figure_5_Tier1_PFAM_annotations_and_PSORT.png"
)

mhc <- readxl::read_xlsx(
  file_overview
)

protein_id_col <- "Member"
psort_col <- "PSORT_Localization_member"

top_n <- 15
keep_other <- TRUE

pfam_cols <- names(mhc) %>%
  stringr::str_subset("^pfam_\\d+$")

pfam_cols <- pfam_cols[
  order(
    as.integer(
      stringr::str_extract(
        pfam_cols,
        "\\d+"
      )
    )
  )
]

psort_levels <- c(
  "Extracellular",
  "OuterMembrane",
  "Periplasmic",
  "CytoplasmicMembrane",
  "Cytoplasmic",
  "Unknown"
)

psort_palette <- c(
  "Cytoplasmic" = "#45B6B8",
  "CytoplasmicMembrane" = "#6771B8",
  "Extracellular" = "#0084A7",
  "OuterMembrane" = "#F1DE81",
  "Periplasmic" = "#80146E",
  "Unknown" = "#F2F1E4"
)

clean_psort <- function(x) {
  x <- stringr::str_trim(
    tidyr::replace_na(
      as.character(x),
      "Unknown"
    )
  )

  x <- dplyr::if_else(
    x == "",
    "Unknown",
    x
  )

  dplyr::if_else(
    x %in% psort_levels,
    x,
    "Unknown"
  )
}

annotation_palette_15 <- c(
  "#4E79A7",
  "#F28E2B",
  "#59A14F",
  "#E15759",
  "#76B7B2",
  "#EDC948",
  "#B07AA1",
  "#FF9DA7",
  "#9C755F",
  "#BAB0AC",
  "#1F77B4",
  "#FF7F0E",
  "#2CA02C",
  "#D62728",
  "#9467BD"
)

other_annotation_color <- "grey90"

theme_clean <- function(base_size = 10) {
  ggplot2::theme_minimal(
    base_size = base_size
  ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.key = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      plot.margin = ggplot2::margin(
        6,
        10,
        6,
        6
      )
    )
}

df0 <- mhc %>%
  dplyr::transmute(
    protein_id = as.character(
      .data[[protein_id_col]]
    ),
    PSORT = clean_psort(
      .data[[psort_col]]
    ),
    dplyr::across(
      dplyr::all_of(pfam_cols),
      as.character
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$protein_id),
    .data$protein_id != ""
  ) %>%
  dplyr::mutate(
    PSORT = factor(
      .data$PSORT,
      levels = psort_levels
    )
  )

protein_pfam_long <- df0 %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(pfam_cols),
    names_to = "pfam_col",
    values_to = "pfam_entry"
  ) %>%
  dplyr::mutate(
    pfam_entry = stringr::str_trim(
      tidyr::replace_na(
        .data$pfam_entry,
        ""
      )
    ),
    pfam_rank = as.integer(
      stringr::str_extract(
        .data$pfam_col,
        "\\d+"
      )
    )
  ) %>%
  dplyr::filter(
    .data$pfam_entry != "",
    .data$pfam_entry != "0"
  ) %>%
  dplyr::mutate(
    PFAM_ID = stringr::str_extract(
      .data$pfam_entry,
      "PF\\d{5}"
    ),
    PFAM_description = stringr::str_trim(
      stringr::str_replace(
        .data$pfam_entry,
        "^PF\\d{5}\\s*\\|\\s*",
        ""
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$PFAM_ID)
  ) %>%
  dplyr::mutate(
    PFAM_description = dplyr::if_else(
      .data$PFAM_description == "" |
        is.na(.data$PFAM_description),
      .data$PFAM_ID,
      .data$PFAM_description
    )
  ) %>%
  dplyr::distinct(
    .data$protein_id,
    .data$PSORT,
    .data$pfam_rank,
    .data$PFAM_ID,
    .data$PFAM_description
  )

first_annotation <- protein_pfam_long %>%
  dplyr::arrange(
    .data$protein_id,
    .data$pfam_rank
  ) %>%
  dplyr::group_by(
    .data$protein_id
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup()

total_selected <- dplyr::n_distinct(
  first_annotation$protein_id
)

annotation_stats <- first_annotation %>%
  dplyr::count(
    .data$PFAM_description,
    name = "n_proteins"
  ) %>%
  dplyr::mutate(
    pct = 100 *
      .data$n_proteins /
      total_selected
  ) %>%
  dplyr::arrange(
    dplyr::desc(.data$pct)
  )

top_levels <- annotation_stats %>%
  dplyr::slice_head(
    n = top_n
  ) %>%
  dplyr::pull(
    .data$PFAM_description
  )

hits_labeled <- first_annotation %>%
  dplyr::mutate(
    annotation = dplyr::if_else(
      .data$PFAM_description %in% top_levels,
      .data$PFAM_description,
      "Other annotations"
    )
  )

if (!keep_other) {
  hits_labeled <- hits_labeled %>%
    dplyr::filter(
      .data$annotation != "Other annotations"
    )
}

donut_df <- hits_labeled %>%
  dplyr::count(
    .data$annotation,
    name = "n_proteins"
  ) %>%
  dplyr::mutate(
    pct = 100 *
      .data$n_proteins /
      sum(.data$n_proteins)
  ) %>%
  dplyr::arrange(
    dplyr::desc(.data$pct)
  ) %>%
  dplyr::mutate(
    annotation = factor(
      .data$annotation,
      levels = .data$annotation
    )
  )

annotation_levels <- levels(
  donut_df$annotation
)

annotation_fill_map <- stats::setNames(
  rep(
    other_annotation_color,
    length(annotation_levels)
  ),
  annotation_levels
)

top_annotations_present <- intersect(
  annotation_levels,
  top_levels
)

annotation_fill_map[
  top_annotations_present
] <- annotation_palette_15[
  seq_along(top_annotations_present)
]

annotation_fill_map[
  "Other annotations"
] <- other_annotation_color

p_donut <- ggplot2::ggplot(
  donut_df,
  ggplot2::aes(
    x = 2,
    y = .data$pct,
    fill = .data$annotation
  )
) +
  ggplot2::geom_col(
    width = 0.95,
    color = "grey15",
    linewidth = 0.25
  ) +
  ggplot2::coord_polar(
    theta = "y"
  ) +
  ggplot2::xlim(
    0.8,
    2.55
  ) +
  ggplot2::scale_fill_manual(
    values = annotation_fill_map,
    name = "Annotation"
  ) +
  ggplot2::labs(
    title = "Annotation composition (Top 15)",
    y = NULL,
    x = NULL
  ) +
  theme_clean(
    base_size = 11
  ) +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "grey20",
      fill = NA,
      linewidth = 0.5
    ),
    legend.position = "right",
    legend.text = ggplot2::element_text(
      size = 8
    )
  )

psort_df <- hits_labeled %>%
  dplyr::filter(
    .data$annotation != "Other annotations"
  ) %>%
  dplyr::mutate(
    annotation = factor(
      .data$annotation,
      levels = top_levels
    ),
    PSORT = factor(
      as.character(.data$PSORT),
      levels = psort_levels
    )
  ) %>%
  dplyr::count(
    .data$annotation,
    .data$PSORT,
    name = "n"
  ) %>%
  dplyr::group_by(
    .data$annotation
  ) %>%
  dplyr::mutate(
    pct = 100 *
      .data$n /
      sum(.data$n)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    !is.na(.data$annotation)
  )

annotations_in_psort <- psort_df %>%
  dplyr::distinct(
    .data$annotation
  ) %>%
  dplyr::pull(
    .data$annotation
  ) %>%
  as.character()

strip_colors <- unname(
  annotation_fill_map[
    annotations_in_psort
  ]
)

psort_df <- psort_df %>%
  dplyr::mutate(
    annotation = factor(
      as.character(.data$annotation),
      levels = annotations_in_psort
    )
  )

p_psort <- ggplot2::ggplot(
  psort_df,
  ggplot2::aes(
    x = 1,
    y = .data$pct,
    fill = .data$PSORT
  )
) +
  ggplot2::geom_col(
    width = 1,
    color = "grey15",
    linewidth = 0.22
  ) +
  ggplot2::coord_polar(
    theta = "y"
  ) +
  ggh4x::facet_wrap2(
    ~ annotation,
    ncol = 5,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = strip_colors,
        colour = "grey25",
        linewidth = 0.45
      ),
      text_x = ggh4x::elem_list_text(
        colour = rep(
          "black",
          length(strip_colors)
        ),
        face = "bold",
        size = 7
      )
    )
  ) +
  ggplot2::scale_fill_manual(
    values = psort_palette,
    drop = FALSE,
    name = "PSORT"
  ) +
  ggplot2::labs(
    title = "PSORT localization within each annotation (Top 15)",
    y = NULL,
    x = NULL
  ) +
  theme_clean(
    base_size = 9
  ) +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.35
    ),
    panel.spacing = grid::unit(
      0.8,
      "lines"
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 9
    ),
    legend.text = ggplot2::element_text(
      size = 8
    )
  )

final_plot <- (
  p_donut +
    p_psort
) +
  patchwork::plot_layout(
    widths = c(
      1,
      2
    )
  ) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      )
    )
  )

print(final_plot)

ggplot2::ggsave(
  filename = out_svg,
  plot = final_plot,
  width = 14,
  height = 7,
  units = "in",
  device = svglite::svglite,
  bg = "white",
  limitsize = FALSE
)

if (capabilities("cairo")) {
  ggplot2::ggsave(
    filename = out_pdf,
    plot = final_plot,
    width = 14,
    height = 7,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
}

ggplot2::ggsave(
  filename = out_png,
  plot = final_plot,
  width = 14,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
