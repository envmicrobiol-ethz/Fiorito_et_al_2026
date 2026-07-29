#!/usr/bin/env Rscript

###############################################################################
# FIGURE 1A — GASES AND TERMINAL ELECTRON ACCEPTORS
#
# DESCRIPTION
#
# Generates the depth profiles of nitrate, sulfate, total iron, total
# manganese, methane and carbon dioxide for Björsmossen and Norra Romyren.
#
# INPUT
#
# No external input file is required.
#
# The biogeochemical measurements used for this panel are defined in the
# `raw_data` table below.
#
# Measurements reported as non-detected in the original source table are
# encoded as NA_real_. These individual measurements are removed before
# plotting and are never converted to zero.
#
# Genuine measured concentrations close to zero are retained.
#
# UNITS
#
# Nitrate, sulfate, total Fe and total Mn:
#   µM
#
# Methane and carbon dioxide:
#   source values are µmol L^-1 and are divided by 1000 for plotting in mM
#
# OUTPUT
#
# Figure1A_BM_gases_and_terminal_electron_acceptors.svg
# Figure1A_NR_gases_and_terminal_electron_acceptors.svg
#
# PDF versions are also generated when Cairo support is available.
#
# USAGE
#
# Rscript Figure1A_gases_and_terminal_electron_acceptors.R
#
# Optional output directory:
#
# Rscript Figure1A_gases_and_terminal_electron_acceptors.R path/to/output
###############################################################################

###############################################################################
# 1. REQUIRED PACKAGES
###############################################################################

required_packages <- c(
  "dplyr",
  "ggplot2",
  "patchwork",
  "scales",
  "svglite",
  "tibble",
  "tidyr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "Missing required R packages:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

###############################################################################
# 2. OUTPUT DIRECTORY
###############################################################################

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) > 1) {
  stop(
    "Provide zero or one command-line argument: [output_directory]."
  )
}

output_dir <- if (length(args) == 1) {
  args[[1]]
} else {
  "."
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
  stop(
    paste(
      "Could not create output directory:",
      output_dir
    )
  )
}

bm_svg <- file.path(
  output_dir,
  "Figure1A_BM_gases_and_terminal_electron_acceptors.svg"
)

nr_svg <- file.path(
  output_dir,
  "Figure1A_NR_gases_and_terminal_electron_acceptors.svg"
)

bm_pdf <- file.path(
  output_dir,
  "Figure1A_BM_gases_and_terminal_electron_acceptors.pdf"
)

nr_pdf <- file.path(
  output_dir,
  "Figure1A_NR_gases_and_terminal_electron_acceptors.pdf"
)

session_info_file <- file.path(
  output_dir,
  "Figure1A_sessionInfo.txt"
)

###############################################################################
# 3. SOURCE DATA
#
# Dashes in the original source table are represented as NA_real_.
###############################################################################

raw_data <- tibble::tribble(
  ~Peatland, ~Site, ~Depth,
  ~CH4_umol_L, ~CO2_umol_L,
  ~NO3_uM, ~SO4_uM,
  ~Mn_total_uM, ~Fe_total_uM,

  # BJÖRSMOSSEN — SITE 1
  "Björsmossen", 1, "0–10",
  199.0, 776.9,
  6.0, 1.2,
  0.09, 2.7,

  "Björsmossen", 1, "20–30",
  542.2, 1993.4,
  7.9, NA_real_,
  0.10, 2.2,

  "Björsmossen", 1, "50–60",
  634.5, 2256.6,
  6.3, 0.6,
  0.09, 2.3,

  # BJÖRSMOSSEN — SITE 2
  "Björsmossen", 2, "0–10",
  367.9, 624.4,
  1.7, 0.6,
  0.21, 7.0,

  "Björsmossen", 2, "20–30",
  422.6, 867.5,
  8.7, NA_real_,
  0.08, 1.8,

  "Björsmossen", 2, "50–60",
  295.1, 796.2,
  2.7, 1.6,
  0.06, 1.7,

  # BJÖRSMOSSEN — SITE 3
  "Björsmossen", 3, "0–10",
  263.4, 539.1,
  7.2, NA_real_,
  0.16, 3.6,

  "Björsmossen", 3, "20–30",
  517.8, 1657.7,
  6.1, 1.2,
  0.03, 1.6,

  "Björsmossen", 3, "50–60",
  381.6, 1376.1,
  1.8, 1.0,
  0.03, 1.8,

  # NORRA ROMYREN — SITE 1
  "Norra Romyren", 1, "0–10",
  260.6, 1615.6,
  4.5, NA_real_,
  0.16, 4.6,

  "Norra Romyren", 1, "20–30",
  373.5, 1985.4,
  6.8, 0.6,
  0.16, 3.7,

  "Norra Romyren", 1, "50–60",
  606.5, 3027.4,
  6.5, 0.2,
  0.09, 2.2,

  # NORRA ROMYREN — SITE 2
  "Norra Romyren", 2, "0–10",
  341.5, 1909.8,
  5.7, 1.1,
  0.03, 1.3,

  "Norra Romyren", 2, "20–30",
  403.5, 2388.6,
  5.6, 1.1,
  0.04, 1.9,

  "Norra Romyren", 2, "50–60",
  619.6, 3600.7,
  7.6, 2.7,
  0.07, 2.0,

  # NORRA ROMYREN — SITE 3
  "Norra Romyren", 3, "0–10",
  180.3, 1663.5,
  15.8, 2.7,
  0.42, 11.9,

  "Norra Romyren", 3, "20–30",
  191.6, 1990.7,
  10.3, 3.0,
  0.16, 5.2,

  "Norra Romyren", 3, "50–60",
  177.9, 2928.5,
  8.9, 1.0,
  0.21, 19.2
)

###############################################################################
# 4. VALIDATE SOURCE DATA
###############################################################################

if (nrow(raw_data) != 18) {
  stop(
    "Expected exactly 18 peatland × site × depth rows."
  )
}

if (sum(is.na(raw_data$SO4_uM)) != 4) {
  stop(
    "Expected exactly four non-detected sulfate measurements."
  )
}

columns_without_expected_missing_values <- c(
  "CH4_umol_L",
  "CO2_umol_L",
  "NO3_uM",
  "Mn_total_uM",
  "Fe_total_uM"
)

if (any(is.na(raw_data[columns_without_expected_missing_values]))) {
  stop(
    "Unexpected missing values were found outside the sulfate column."
  )
}

expected_missing_so4 <- tibble::tribble(
  ~Peatland, ~Site, ~Depth,

  "Björsmossen", 1, "20–30",
  "Björsmossen", 2, "20–30",
  "Björsmossen", 3, "0–10",
  "Norra Romyren", 1, "0–10"
) |>
  dplyr::arrange(
    .data$Peatland,
    .data$Site,
    .data$Depth
  )

observed_missing_so4 <- raw_data |>
  dplyr::filter(
    is.na(.data$SO4_uM)
  ) |>
  dplyr::select(
    .data$Peatland,
    .data$Site,
    .data$Depth
  ) |>
  dplyr::arrange(
    .data$Peatland,
    .data$Site,
    .data$Depth
  )

if (!identical(
  observed_missing_so4,
  expected_missing_so4
)) {
  stop(
    paste(
      "The positions of the non-detected sulfate measurements",
      "do not match the original source table."
    )
  )
}

###############################################################################
# 5. PREPARE SAMPLE FACTORS
###############################################################################

# Reverse the levels so that 0–10 cm appears at the top.
depth_levels <- c(
  "50–60",
  "20–30",
  "0–10"
)

site_order <- c(
  "site 1",
  "site 2",
  "site 3"
)

prepared <- raw_data |>
  dplyr::mutate(
    Depth = factor(
      .data$Depth,
      levels = depth_levels
    ),
    Site_label = factor(
      paste(
        "site",
        .data$Site
      ),
      levels = site_order
    )
  )

###############################################################################
# 6. PREPARE SOLUTE DATA
#
# Individual non-detected measurements are removed here.
# The other measurements from the same sample are retained.
###############################################################################

solutes_long <- prepared |>
  dplyr::select(
    .data$Peatland,
    .data$Site_label,
    .data$Depth,
    .data$NO3_uM,
    .data$SO4_uM,
    .data$Mn_total_uM,
    .data$Fe_total_uM
  ) |>
  tidyr::pivot_longer(
    cols = c(
      "NO3_uM",
      "SO4_uM",
      "Mn_total_uM",
      "Fe_total_uM"
    ),
    names_to = "Analyte",
    values_to = "Value_uM"
  ) |>
  dplyr::filter(
    !is.na(.data$Value_uM)
  ) |>
  dplyr::mutate(
    Analyte = dplyr::recode(
      .data$Analyte,
      "NO3_uM" = "NO3",
      "SO4_uM" = "SO4",
      "Mn_total_uM" = "Mn_total",
      "Fe_total_uM" = "Fe_total"
    )
  )

###############################################################################
# 7. PREPARE GAS DATA
#
# Gas concentrations are converted from µmol L^-1 to mM.
###############################################################################

gases_long <- prepared |>
  dplyr::select(
    .data$Peatland,
    .data$Site_label,
    .data$Depth,
    .data$CH4_umol_L,
    .data$CO2_umol_L
  ) |>
  tidyr::pivot_longer(
    cols = c(
      "CH4_umol_L",
      "CO2_umol_L"
    ),
    names_to = "Analyte",
    values_to = "Value_umol_L"
  ) |>
  dplyr::filter(
    !is.na(.data$Value_umol_L)
  ) |>
  dplyr::mutate(
    Analyte = dplyr::recode(
      .data$Analyte,
      "CH4_umol_L" = "CH4",
      "CO2_umol_L" = "CO2"
    ),
    Value_mM = .data$Value_umol_L / 1000
  )

###############################################################################
# 8. APPEARANCE
###############################################################################

plot_font <- Sys.getenv(
  "FIGURE_FONT",
  unset = "Helvetica"
)

analyte_colors <- c(
  "NO3" = "#1B9E77",
  "SO4" = "#D95F02",
  "Fe_total" = "#8C510A",
  "Mn_total" = "#756BB1",
  "CH4" = "#66A61E",
  "CO2" = "#636363"
)

analyte_order <- c(
  "NO3",
  "SO4",
  "Fe_total",
  "Mn_total",
  "CH4",
  "CO2"
)

analyte_labels <- c(
  expression(NO[3]^"-"),
  expression(SO[4]^"2-"),
  "total Fe",
  "total Mn",
  expression(CH[4]),
  expression(CO[2])
)

site_shapes <- c(
  "site 1" = 15,
  "site 2" = 17,
  "site 3" = 16
)

point_position <- ggplot2::position_jitter(
  width = 0,
  height = 0.055,
  seed = 2026
)

###############################################################################
# 9. THEME
###############################################################################

theme_publication <- function(base_size = 15) {

  ggplot2::theme_classic(
    base_size = base_size,
    base_family = plot_font
  ) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(
        linewidth = 0.8,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.8,
        colour = "black"
      ),
      axis.ticks.length = grid::unit(
        3,
        "pt"
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(
        13,
        "pt"
      ),
      legend.key.width = grid::unit(
        20,
        "pt"
      ),
      legend.text = ggplot2::element_text(
        size = 11
      ),
      text = ggplot2::element_text(
        colour = "black",
        family = plot_font
      )
    )
}

###############################################################################
# 10. SCALES AND GUIDES
###############################################################################

colour_scale <- ggplot2::scale_colour_manual(
  values = analyte_colors,
  breaks = analyte_order,
  labels = analyte_labels,
  drop = FALSE,
  name = NULL
)

shape_scale <- ggplot2::scale_shape_manual(
  values = site_shapes,
  breaks = site_order,
  labels = site_order,
  drop = FALSE,
  name = NULL
)

common_guides <- ggplot2::guides(
  colour = ggplot2::guide_legend(
    order = 1,
    override.aes = list(
      shape = 16,
      size = 4.3,
      alpha = 1
    )
  ),
  shape = ggplot2::guide_legend(
    order = 2,
    override.aes = list(
      colour = "black",
      size = 5.2,
      alpha = 1
    )
  )
)

###############################################################################
# 11. SOLUTE PANEL
###############################################################################

make_solute_plot <- function(data) {

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$Value_uM,
      y = .data$Depth,
      colour = .data$Analyte,
      shape = .data$Site_label
    )
  ) +
    ggplot2::geom_point(
      size = 4.8,
      alpha = 0.95,
      stroke = 0.7,
      position = point_position
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        0,
        20
      ),
      breaks = c(
        0,
        5,
        10,
        15,
        20
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0
        )
      ),
      oob = scales::oob_censor
    ) +
    colour_scale +
    shape_scale +
    common_guides +
    ggplot2::labs(
      x = expression(
        Concentration~"(" * mu * "M)"
      ),
      y = "Depth (cm)"
    ) +
    theme_publication()
}

###############################################################################
# 12. GAS PANEL
###############################################################################

make_gas_plot <- function(data) {

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$Value_mM,
      y = .data$Depth,
      colour = .data$Analyte,
      shape = .data$Site_label
    )
  ) +
    ggplot2::geom_point(
      size = 4.8,
      alpha = 0.95,
      stroke = 0.7,
      position = point_position
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        0,
        4
      ),
      breaks = 0:4,
      expand = ggplot2::expansion(
        mult = c(
          0,
          0
        )
      ),
      oob = scales::oob_censor
    ) +
    colour_scale +
    shape_scale +
    common_guides +
    ggplot2::labs(
      x = "Concentration (mM)",
      y = NULL
    ) +
    theme_publication() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
}

###############################################################################
# 13. COMBINE PANELS
###############################################################################

combine_panels <- function(
    solute_plot,
    gas_plot
) {

  solute_plot +
    gas_plot +
    patchwork::plot_layout(
      widths = c(
        1,
        1
      ),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "top",
      text = ggplot2::element_text(
        family = plot_font
      )
    )
}

###############################################################################
# 14. CREATE FIGURES
###############################################################################

bm_solutes <- solutes_long |>
  dplyr::filter(
    .data$Peatland == "Björsmossen"
  )

bm_gases <- gases_long |>
  dplyr::filter(
    .data$Peatland == "Björsmossen"
  )

nr_solutes <- solutes_long |>
  dplyr::filter(
    .data$Peatland == "Norra Romyren"
  )

nr_gases <- gases_long |>
  dplyr::filter(
    .data$Peatland == "Norra Romyren"
  )

p_bm <- combine_panels(
  make_solute_plot(
    bm_solutes
  ),
  make_gas_plot(
    bm_gases
  )
)

p_nr <- combine_panels(
  make_solute_plot(
    nr_solutes
  ),
  make_gas_plot(
    nr_gases
  )
)

print(
  p_bm
)

print(
  p_nr
)

###############################################################################
# 15. SAVE SVG FILES
###############################################################################

ggplot2::ggsave(
  filename = bm_svg,
  plot = p_bm,
  device = svglite::svglite,
  width = 8.8,
  height = 5.0,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  filename = nr_svg,
  plot = p_nr,
  device = svglite::svglite,
  width = 8.8,
  height = 5.0,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# 16. SAVE PDF FILES
###############################################################################

if (capabilities("cairo")) {

  ggplot2::ggsave(
    filename = bm_pdf,
    plot = p_bm,
    device = grDevices::cairo_pdf,
    width = 8.8,
    height = 5.0,
    units = "in",
    bg = "white"
  )

  ggplot2::ggsave(
    filename = nr_pdf,
    plot = p_nr,
    device = grDevices::cairo_pdf,
    width = 8.8,
    height = 5.0,
    units = "in",
    bg = "white"
  )

} else {

  warning(
    "PDF files were not saved because Cairo support is unavailable."
  )
}

###############################################################################
# 17. SESSION INFORMATION
###############################################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_info_file
)

cat(
  "\nFigure 1A generation completed successfully.\n"
)

cat(
  "Non-detected measurements were excluded from plotting.\n"
)
