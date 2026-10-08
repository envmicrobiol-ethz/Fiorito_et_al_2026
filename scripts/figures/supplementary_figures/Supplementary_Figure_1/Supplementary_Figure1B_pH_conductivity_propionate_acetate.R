#!/usr/bin/env Rscript

###############################################################################
# EXTENDED DATA FIGURE 1B — BIOGEOCHEMICAL DEPTH PROFILES
#
# DESCRIPTION
# Generates pH, conductivity, propionate and acetate depth profiles for
# Björsmossen and Norra Romyren.
#
# INPUT
# No external input file. Source measurements are defined below.
#
# OUTPUT
# EDFigure1B_BM_pH_conductivity_propionate_acetate.svg
# EDFigure1B_NR_pH_conductivity_propionate_acetate.svg
#
# USAGE
# Rscript EDFigure1B_pH_conductivity_propionate_acetate.R [output_directory]
###############################################################################

required_packages <- c(
  "dplyr",
  "ggplot2",
  "patchwork",
  "scales",
  "svglite",
  "tibble"
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
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

###############################################################################
# INPUT AND OUTPUT
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 1) {
  stop(
    "Usage: Rscript EDFigure1B_pH_conductivity_propionate_acetate.R [output_directory]"
  )
}

output_dir <- if (length(args) == 1) args[[1]] else "."

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

bm_svg <- file.path(
  output_dir,
  "EDFigure1B_BM_pH_conductivity_propionate_acetate.svg"
)

lm_svg <- file.path(
  output_dir,
  "EDFigure1B_LM_pH_conductivity_propionate_acetate.svg"
)

nr_svg <- file.path(
  output_dir,
  "EDFigure1B_NR_pH_conductivity_propionate_acetate.svg"
)

###############################################################################
# SOURCE DATA
#
# Propionate and acetate are expressed in mM.
# Non-detected or unavailable measurements are represented as NA.
###############################################################################

raw_data <- tibble::tribble(
  ~Peatland, ~Site, ~Depth, ~pH,
  ~Conductivity_uS_cm, ~Propionate_mM, ~Acetate_mM,

  # Björsmossen — site 1
  "Björsmossen", 1, "0–10",  4.2, 44.7, 0.7, NA_real_,
  "Björsmossen", 1, "20–30", 4.1, 40.7, 1.0, NA_real_,
  "Björsmossen", 1, "50–60", 4.0, 43.7, 0.7, NA_real_,

  # Björsmossen — site 2
  "Björsmossen", 2, "0–10",  4.1, 46.7, 0.7, NA_real_,
  "Björsmossen", 2, "20–30", 4.2, 30.0, 0.7, NA_real_,
  "Björsmossen", 2, "50–60", 4.4, 21.5, 0.7, NA_real_,

  # Björsmossen — site 3
  "Björsmossen", 3, "0–10",  4.0, 47.7, 0.7, 0.02,
  "Björsmossen", 3, "20–30", 4.1, 41.4, 0.9, NA_real_,
  "Björsmossen", 3, "50–60", 4.1, 39.7, 0.8, NA_real_,

  # Norra Romyren — site 1
  "Norra Romyren", 1, "0–10",  4.1, 44.4, 0.7, NA_real_,
  "Norra Romyren", 1, "20–30", 4.0, 55.7, 1.2, NA_real_,
  "Norra Romyren", 1, "50–60", 4.0, 49.0, 0.7, NA_real_,

  # Norra Romyren — site 2
  "Norra Romyren", 2, "0–10",  4.3, 34.3, 1.3, NA_real_,
  "Norra Romyren", 2, "20–30", 3.9, 53.3, 0.7, NA_real_,
  "Norra Romyren", 2, "50–60", 4.0, 53.6, 0.7, NA_real_,

  # Norra Romyren — site 3
  "Norra Romyren", 3, "0–10",  4.4, 35.1, NA_real_, 0.20,
  "Norra Romyren", 3, "20–30", 4.1, 41.3, 2.9, 0.10,
  "Norra Romyren", 3, "50–60", 4.0, 46.3, NA_real_, 0.01,

  # Lungsmossen — site 1
  "Lungsmossen", 1, "0–10",  4.1, 45.0, NA_real_, NA_real_,
  "Lungsmossen", 1, "20–30", 4.0, 60.8, NA_real_, NA_real_,
  "Lungsmossen", 1, "50–60", 4.0, 64.0, NA_real_, NA_real_,

  # Lungsmossen — site 2
  "Lungsmossen", 2, "0–10",  4.0, 51.8, NA_real_, NA_real_,
  "Lungsmossen", 2, "20–30", 4.1, 44.2, NA_real_, NA_real_,
  "Lungsmossen", 2, "50–60", 4.1, 43.9, NA_real_, NA_real_,

  # Lungsmossen — site 3
  "Lungsmossen", 3, "0–10",  4.0, 52.0, NA_real_, NA_real_,
  "Lungsmossen", 3, "20–30", 4.1, 58.1, NA_real_, NA_real_,
  "Lungsmossen", 3, "50–60", 4.0, 59.4, NA_real_, NA_real_
)

if (nrow(raw_data) != 27) {
  stop("Expected 27 peatland × site × depth rows.")
}

if (anyNA(raw_data$pH) || anyNA(raw_data$Conductivity_uS_cm)) {
  stop("Unexpected missing pH or conductivity values.")
}

###############################################################################
# FACTORS AND APPEARANCE
###############################################################################

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

site_shapes <- c(
  "site 1" = 15,
  "site 2" = 17,
  "site 3" = 16
)

prepared <- raw_data |>
  dplyr::mutate(
    Depth = factor(
      .data$Depth,
      levels = depth_levels
    ),
    Site_label = factor(
      paste("site", .data$Site),
      levels = site_order
    )
  )

analyte_colors <- c(
  "pH" = "#B2182B",
  "Conductivity" = "#2171B5",
  "Propionate" = "#CC2F8A",
  "Acetate" = "#D9A400"
)

point_position <- ggplot2::position_jitter(
  width = 0,
  height = 0.055,
  seed = 2026
)

theme_publication <- function() {
  ggplot2::theme_classic(
    base_size = 15,
    base_family = "Helvetica"
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
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = 13,
        face = "plain",
        hjust = 0.5
      ),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = 11
      )
    )
}

###############################################################################
# PANEL FUNCTION
###############################################################################

make_depth_panel <- function(
    data,
    value_column,
    analyte_name,
    x_label,
    x_limits,
    x_breaks,
    show_y_axis = TRUE
) {

  axis_template <- tibble::tibble(
    x = rep(
      x_limits[1],
      length(depth_levels)
    ),
    Depth = factor(
      depth_levels,
      levels = depth_levels
    )
  )

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[value_column]],
      y = .data$Depth,
      shape = .data$Site_label
    )
  ) +
    ggplot2::geom_blank(
      data = axis_template,
      ggplot2::aes(
        x = .data$x,
        y = .data$Depth
      ),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      colour = analyte_colors[[analyte_name]],
      size = 4.8,
      alpha = 0.95,
      stroke = 0.7,
      position = point_position,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks,
      expand = ggplot2::expansion(
        mult = c(0.02, 0.04)
      ),
      oob = scales::oob_censor
    ) +
    ggplot2::scale_y_discrete(
      limits = depth_levels,
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = site_shapes,
      breaks = site_order,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = analyte_name,
      x = x_label,
      y = if (show_y_axis) "Depth (cm)" else NULL,
      shape = NULL
    ) +
    theme_publication()

  if (!show_y_axis) {
    plot <- plot +
      ggplot2::theme(
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      )
  }

  plot
}

###############################################################################
# COMBINE PANELS
###############################################################################

make_four_panel_figure <- function(data) {

  p_ph <- make_depth_panel(
    data,
    "pH",
    "pH",
    "pH",
    c(3.8, 4.6),
    c(3.8, 4.0, 4.2, 4.4, 4.6),
    TRUE
  )

  p_conductivity <- make_depth_panel(
    data,
    "Conductivity_uS_cm",
    "Conductivity",
    expression(Conductivity~"(" * mu * "S " * cm^{-1} * ")"),
    c(20, 70),
    c(20, 30, 40, 50, 60, 70),
    FALSE
  )

  p_propionate <- make_depth_panel(
    data,
    "Propionate_mM",
    "Propionate",
    "Concentration (mM)",
    c(0, 3),
    c(0, 1, 2, 3),
    FALSE
  )

  p_acetate <- make_depth_panel(
    data,
    "Acetate_mM",
    "Acetate",
    "Concentration (mM)",
    c(0, 0.3),
    c(0, 0.1, 0.2, 0.3),
    FALSE
  )

  p_ph +
    p_conductivity +
    p_propionate +
    p_acetate +
    patchwork::plot_layout(
      nrow = 1,
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "top"
    )
}

###############################################################################
# CREATE AND SAVE FIGURES
###############################################################################

bm_plot <- make_four_panel_figure(
  dplyr::filter(
    prepared,
    .data$Peatland == "Björsmossen"
  )
)

lm_plot <- make_four_panel_figure(
  dplyr::filter(
    prepared,
    .data$Peatland == "Lungsmossen"
  )
)

nr_plot <- make_four_panel_figure(
  dplyr::filter(
    prepared,
    .data$Peatland == "Norra Romyren"
  )
)

ggplot2::ggsave(
  bm_svg,
  bm_plot,
  device = svglite::svglite,
  width = 12,
  height = 4.8,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  lm_svg,
  lm_plot,
  device = svglite::svglite,
  width = 12,
  height = 4.8,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  nr_svg,
  nr_plot,
  device = svglite::svglite,
  width = 12,
  height = 4.8,
  units = "in",
  bg = "white"
)

cat(
  "\nExtended Data Figure 1B completed.\n",
  sep = ""
)
