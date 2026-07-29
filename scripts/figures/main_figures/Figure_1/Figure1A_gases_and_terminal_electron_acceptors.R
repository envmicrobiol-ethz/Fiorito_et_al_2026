#!/usr/bin/env Rscript

###############################################################################
# FIGURE 1A — GASES AND TERMINAL ELECTRON ACCEPTORS
#
# Dashes or blank cells in the source table are encoded as NA and are removed
# before plotting. They are never converted to zero. Genuine measured values
# close to zero are retained.
#
# Usage:
#   Rscript 01_Figure1A_gases_and_terminal_electron_acceptors.R [output_dir]
###############################################################################

required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "scales", "svglite", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1) stop("Provide zero or one argument: [output_dir].")
output_dir <- if (length(args) == 1) args[[1]] else "."
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bm_svg <- file.path(output_dir, "Figure1A_BM_gases_and_TEA.svg")
nr_svg <- file.path(output_dir, "Figure1A_NR_gases_and_TEA.svg")
combined_svg <- file.path(output_dir, "Figure1A_gases_and_TEA.svg")
bm_pdf <- file.path(output_dir, "Figure1A_BM_gases_and_TEA.pdf")
nr_pdf <- file.path(output_dir, "Figure1A_NR_gases_and_TEA.pdf")
combined_pdf <- file.path(output_dir, "Figure1A_gases_and_TEA.pdf")
source_data_file <- file.path(output_dir, "Figure1A_source_data.tsv")
non_detect_file <- file.path(output_dir, "Figure1A_non_detected_measurements.tsv")

# Source data. A dash in the source table is represented as NA_real_.
raw_data <- tibble::tribble(
  ~Peatland, ~Site, ~Depth,
  ~CH4_umol_L, ~CO2_umol_L,
  ~NO3_uM, ~SO4_uM,
  ~Mn_total_uM, ~Fe_total_uM,

  "Björsmossen", 1, "0–10",  199.0,  776.9,  6.0, 1.2,      0.09,  2.7,
  "Björsmossen", 1, "20–30", 542.2, 1993.4,  7.9, NA_real_, 0.10,  2.2,
  "Björsmossen", 1, "50–60", 634.5, 2256.6,  6.3, 0.6,      0.09,  2.3,

  "Björsmossen", 2, "0–10",  367.9,  624.4,  1.7, 0.6,      0.21,  7.0,
  "Björsmossen", 2, "20–30", 422.6,  867.5,  8.7, NA_real_, 0.08,  1.8,
  "Björsmossen", 2, "50–60", 295.1,  796.2,  2.7, 1.6,      0.06,  1.7,

  "Björsmossen", 3, "0–10",  263.4,  539.1,  7.2, NA_real_, 0.16,  3.6,
  "Björsmossen", 3, "20–30", 517.8, 1657.7,  6.1, 1.2,      0.03,  1.6,
  "Björsmossen", 3, "50–60", 381.6, 1376.1,  1.8, 1.0,      0.03,  1.8,

  "Norra Romyren", 1, "0–10",  260.6, 1615.6,  4.5, NA_real_, 0.16,  4.6,
  "Norra Romyren", 1, "20–30", 373.5, 1985.4,  6.8, 0.6,      0.16,  3.7,
  "Norra Romyren", 1, "50–60", 606.5, 3027.4,  6.5, 0.2,      0.09,  2.2,

  "Norra Romyren", 2, "0–10",  341.5, 1909.8,  5.7, 1.1,      0.03,  1.3,
  "Norra Romyren", 2, "20–30", 403.5, 2388.6,  5.6, 1.1,      0.04,  1.9,
  "Norra Romyren", 2, "50–60", 619.6, 3600.7,  7.6, 2.7,      0.07,  2.0,

  "Norra Romyren", 3, "0–10",  180.3, 1663.5, 15.8, 2.7,      0.42, 11.9,
  "Norra Romyren", 3, "20–30", 191.6, 1990.7, 10.3, 3.0,      0.16,  5.2,
  "Norra Romyren", 3, "50–60", 177.9, 2928.5,  8.9, 1.0,      0.21, 19.2
)

if (nrow(raw_data) != 18) stop("Expected 18 peatland × site × depth rows.")
if (sum(is.na(raw_data$SO4_uM)) != 4) stop("Expected four non-detected sulfate measurements.")
if (any(is.na(raw_data[c("CH4_umol_L", "CO2_umol_L", "NO3_uM", "Mn_total_uM", "Fe_total_uM")]))) {
  stop("Unexpected missing values outside sulfate.")
}

expected_missing_so4 <- tibble::tribble(
  ~Peatland, ~Site, ~Depth,
  "Björsmossen", 1, "20–30",
  "Björsmossen", 2, "20–30",
  "Björsmossen", 3, "0–10",
  "Norra Romyren", 1, "0–10"
) |>
  dplyr::arrange(.data$Peatland, .data$Site, .data$Depth)

observed_missing_so4 <- raw_data |>
  dplyr::filter(is.na(.data$SO4_uM)) |>
  dplyr::select(.data$Peatland, .data$Site, .data$Depth) |>
  dplyr::arrange(.data$Peatland, .data$Site, .data$Depth)

if (!identical(observed_missing_so4, expected_missing_so4)) {
  stop("The positions of non-detected sulfate measurements do not match the source table.")
}

utils::write.table(raw_data, source_data_file, sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "NA")
utils::write.table(observed_missing_so4, non_detect_file, sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "NA")

depth_levels <- c("50–60", "20–30", "0–10")
site_order <- c("site 1", "site 2", "site 3")
prepared <- raw_data |>
  dplyr::mutate(
    Depth = factor(.data$Depth, levels = depth_levels),
    Site_label = factor(paste("site", .data$Site), levels = site_order)
  )

# Remove only the individual non-detected measurements; never remove an entire
# sample when other analytes were measured for that sample.
solutes_long <- prepared |>
  dplyr::select(.data$Peatland, .data$Site_label, .data$Depth,
                .data$NO3_uM, .data$SO4_uM, .data$Mn_total_uM, .data$Fe_total_uM) |>
  tidyr::pivot_longer(
    cols = c("NO3_uM", "SO4_uM", "Mn_total_uM", "Fe_total_uM"),
    names_to = "Analyte", values_to = "Value_uM"
  ) |>
  dplyr::filter(!is.na(.data$Value_uM)) |>
  dplyr::mutate(
    Analyte = dplyr::recode(
      .data$Analyte,
      "NO3_uM" = "NO3", "SO4_uM" = "SO4",
      "Mn_total_uM" = "Mn_total", "Fe_total_uM" = "Fe_total"
    )
  )

gases_long <- prepared |>
  dplyr::select(.data$Peatland, .data$Site_label, .data$Depth,
                .data$CH4_umol_L, .data$CO2_umol_L) |>
  tidyr::pivot_longer(
    cols = c("CH4_umol_L", "CO2_umol_L"),
    names_to = "Analyte", values_to = "Value_umol_L"
  ) |>
  dplyr::filter(!is.na(.data$Value_umol_L)) |>
  dplyr::mutate(
    Analyte = dplyr::recode(.data$Analyte,
                            "CH4_umol_L" = "CH4", "CO2_umol_L" = "CO2"),
    Value_mM = .data$Value_umol_L / 1000
  )

plot_font <- Sys.getenv("FIGURE_FONT", unset = "Helvetica")
analyte_colors <- c(
  "NO3" = "#1B9E77", "SO4" = "#D95F02", "Fe_total" = "#8C510A",
  "Mn_total" = "#756BB1", "CH4" = "#66A61E", "CO2" = "#636363"
)
analyte_order <- c("NO3", "SO4", "Fe_total", "Mn_total", "CH4", "CO2")
analyte_labels <- c(
  expression(NO[3]^"-"), expression(SO[4]^"2-"), "total Fe", "total Mn",
  expression(CH[4]), expression(CO[2])
)
site_shapes <- c("site 1" = 15, "site 2" = 17, "site 3" = 16)
raw_point_position <- ggplot2::position_jitter(width = 0, height = 0.055, seed = 2026)

theme_pub <- function(base_size = 15) {
  ggplot2::theme_classic(base_size = base_size, base_family = plot_font) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.8, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.8, colour = "black"),
      axis.ticks.length = grid::unit(3, "pt"),
      panel.grid.major.y = ggplot2::element_line(colour = "grey85", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(13, "pt"),
      legend.key.width = grid::unit(20, "pt"),
      legend.text = ggplot2::element_text(size = 11),
      text = ggplot2::element_text(colour = "black", family = plot_font)
    )
}

colour_scale <- ggplot2::scale_colour_manual(
  values = analyte_colors, breaks = analyte_order, labels = analyte_labels,
  drop = FALSE, name = NULL
)
shape_scale <- ggplot2::scale_shape_manual(
  values = site_shapes, breaks = site_order, labels = site_order,
  drop = FALSE, name = NULL
)
common_guides <- ggplot2::guides(
  colour = ggplot2::guide_legend(order = 1,
    override.aes = list(shape = 16, size = 4.3, alpha = 1)),
  shape = ggplot2::guide_legend(order = 2,
    override.aes = list(colour = "black", size = 5.2, alpha = 1))
)

make_solute_plot <- function(data, peatland_label = NULL) {
  ggplot2::ggplot(data, ggplot2::aes(
    x = .data$Value_uM, y = .data$Depth,
    colour = .data$Analyte, shape = .data$Site_label
  )) +
    ggplot2::geom_point(size = 4.8, alpha = 0.95, stroke = 0.7,
                        position = raw_point_position) +
    ggplot2::scale_x_continuous(
      limits = c(0, 20), breaks = c(0, 5, 10, 15, 20),
      expand = ggplot2::expansion(mult = c(0, 0)),
      oob = scales::oob_censor
    ) +
    colour_scale + shape_scale + common_guides +
    ggplot2::labs(title = peatland_label,
                  x = expression(Concentration~"(" * mu * "M)"), y = "Depth (cm)") +
    theme_pub() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, face = "plain"))
}

make_gas_plot <- function(data) {
  ggplot2::ggplot(data, ggplot2::aes(
    x = .data$Value_mM, y = .data$Depth,
    colour = .data$Analyte, shape = .data$Site_label
  )) +
    ggplot2::geom_point(size = 4.8, alpha = 0.95, stroke = 0.7,
                        position = raw_point_position) +
    ggplot2::scale_x_continuous(
      limits = c(0, 4), breaks = 0:4,
      expand = ggplot2::expansion(mult = c(0, 0)),
      oob = scales::oob_censor
    ) +
    colour_scale + shape_scale + common_guides +
    ggplot2::labs(x = "Concentration (mM)", y = NULL) +
    theme_pub() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
}

make_peatland_row <- function(peatland, label) {
  p_solute <- make_solute_plot(
    dplyr::filter(solutes_long, .data$Peatland == peatland), label
  )
  p_gas <- make_gas_plot(
    dplyr::filter(gases_long, .data$Peatland == peatland)
  )
  p_solute + p_gas + patchwork::plot_layout(widths = c(1, 1), guides = "collect")
}

p_bm <- make_peatland_row("Björsmossen", "BM")
p_nr <- make_peatland_row("Norra Romyren", "NR")
p_combined <- (p_bm / p_nr) + patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "top")

print(p_bm)
print(p_nr)
print(p_combined)

ggplot2::ggsave(bm_svg, p_bm, device = svglite::svglite,
                 width = 8.8, height = 5.0, units = "in", bg = "white")
ggplot2::ggsave(nr_svg, p_nr, device = svglite::svglite,
                 width = 8.8, height = 5.0, units = "in", bg = "white")
ggplot2::ggsave(combined_svg, p_combined, device = svglite::svglite,
                 width = 8.8, height = 9.2, units = "in", bg = "white")

if (capabilities("cairo")) {
  ggplot2::ggsave(bm_pdf, p_bm, device = grDevices::cairo_pdf,
                  width = 8.8, height = 5.0, units = "in", bg = "white")
  ggplot2::ggsave(nr_pdf, p_nr, device = grDevices::cairo_pdf,
                  width = 8.8, height = 5.0, units = "in", bg = "white")
  ggplot2::ggsave(combined_pdf, p_combined, device = grDevices::cairo_pdf,
                  width = 8.8, height = 9.2, units = "in", bg = "white")
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure1A_sessionInfo.txt"))
cat("Figure 1A completed. Non-detected measurements were excluded from plotting.\n")
