#!/usr/bin/env Rscript

required_packages <- c("dplyr", "ggplot2", "ggh4x", "readxl", "svglite", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript SCRIPT.R <combined_relative_abundance.xlsx> [output_dir]")
}
input_xlsx <- args[[1]]
output_dir <- if (length(args) == 2) args[[2]] else "."
if (!file.exists(input_xlsx)) stop("Input Excel file not found: ", input_xlsx)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_font <- Sys.getenv("FIGURE_FONT", unset = "Arial")

samples_spatial <- c(
  "BM.1.10", "BM.1.30", "BM.1.60",
  "BM.2.10", "BM.2.30", "BM.2.60",
  "BM.3.30", "BM.3.60",  # BM.3.10 intentionally excluded
  "NR.1.10", "NR.1.30", "NR.1.60",
  "NR.2.10", "NR.2.30", "NR.2.60",
  "NR.3.10", "NR.3.30", "NR.3.60",
  "LM.1.10", "LM.1.30", "LM.1.60",
  "LM.2.10", "LM.2.30", "LM.2.60",
  "LM.3.10", "LM.3.30", "LM.3.60"
)

samples_deep <- c(
  "BM.0.50", "BM.50", "BM.100", "BM.150", "BM.200", "BM.250",
  "BM.300", "BM.350", "BM.400", "BM.450", "BM.490.540",
  "HV.0.25", "HV.50", "HV.100", "HV.150", "HV.200", "HV.250", "HV.265.280",
  "LM.0.25", "LM.50", "LM.100", "LM.150", "LM.200", "LM.250",
  "LM.300", "LM.350", "LM.400", "LM.450", "LM.550", "LM.600",
  "LM.650", "LM.675.700",
  "NR.0.25", "NR.50", "NR.100", "NR.150", "NR.200", "NR.250",
  "NR.300", "NR.325.350"
)

otu_plot_raw <- readxl::read_excel(input_xlsx, sheet = "relative_abundance")
names(otu_plot_raw) <- trimws(sub("^\\ufeff", "", names(otu_plot_raw)))
required_columns <- c("Sample", "Phylum", "Abundance")
missing_columns <- setdiff(required_columns, names(otu_plot_raw))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

df_raw <- otu_plot_raw |>
  dplyr::transmute(
    Sample = trimws(as.character(.data$Sample)),
    Phylum = trimws(tidyr::replace_na(as.character(.data$Phylum), "NA")),
    Abundance = suppressWarnings(as.numeric(.data$Abundance))
  ) |>
  dplyr::mutate(Abundance = dplyr::if_else(is.na(.data$Abundance), 0, .data$Abundance))

archaeal_phyla <- c("p__Halobacterota", "p__Crenarchaeota")
bacterial_phyla <- c(
  "p__Acidobacteriota", "p__Proteobacteria",
  "p__Verrucomicrobiota", "p__Desulfobacterota"
)
main_phyla <- c(archaeal_phyla, bacterial_phyla)
phylum_order <- c(main_phyla, "Others")
phylum_colors <- c(
  "p__Halobacterota" = "#2171B5",
  "p__Crenarchaeota" = "#08306B",
  "p__Acidobacteriota" = "#B2182B",
  "p__Proteobacteria" = "#8C510A",
  "p__Verrucomicrobiota" = "#D8B365",
  "p__Desulfobacterota" = "#C77526",
  "Others" = "grey60"
)

prepare_grouped_data <- function(selected_samples) {
  selected <- df_raw |>
    dplyr::filter(.data$Sample %in% selected_samples)
  missing_samples <- setdiff(selected_samples, unique(selected$Sample))
  if (length(missing_samples) > 0) {
    stop("Samples missing from the plotting table: ", paste(missing_samples, collapse = ", "))
  }
  grouped <- selected |>
    dplyr::mutate(
      Phylum_plot = dplyr::if_else(.data$Phylum %in% main_phyla, .data$Phylum, "Others")
    ) |>
    dplyr::group_by(.data$Sample, .data$Phylum_plot) |>
    dplyr::summarise(Abundance = sum(.data$Abundance, na.rm = TRUE), .groups = "drop") |>
    tidyr::complete(
      Sample = selected_samples,
      Phylum_plot = phylum_order,
      fill = list(Abundance = 0)
    ) |>
    dplyr::mutate(
      Phylum_plot = factor(.data$Phylum_plot, levels = phylum_order),
      Site = sub("\\..*$", "", .data$Sample)
    )
  totals <- grouped |>
    dplyr::group_by(.data$Sample) |>
    dplyr::summarise(total = sum(.data$Abundance), .groups = "drop")
  if (any(totals$total <= 0)) stop("At least one sample has zero total abundance.")
  grouped
}

fill_scale <- ggplot2::scale_fill_manual(
  values = phylum_colors, breaks = phylum_order,
  labels = function(x) gsub("^p__", "", x), drop = FALSE
)
abundance_scale <- ggplot2::scale_x_continuous(
  breaks = seq(0, 1, by = 0.25), labels = function(x) x * 100,
  limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0))
)
common_theme <- ggplot2::theme_bw(base_size = 10, base_family = plot_font) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(colour = "grey85", linewidth = 0.35),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.7),
    strip.background = ggplot2::element_rect(fill = "grey90", colour = "black", linewidth = 0.6),
    axis.text = ggplot2::element_text(size = 8),
    axis.title = ggplot2::element_text(size = 9),
    legend.title = ggplot2::element_blank(),
    legend.position = "right",
    plot.margin = ggplot2::margin(6, 6, 6, 6, unit = "pt")
  )


###############################################################################
# FIGURE 1D — 16S rRNA TAXONOMIC COMPOSITION, DEEP PROFILING
###############################################################################

output_svg <- file.path(output_dir, "Figure1D_16S_deep_composition.svg")
output_pdf <- file.path(output_dir, "Figure1D_16S_deep_composition.pdf")
source_tsv <- file.path(output_dir, "Figure1D_source_data.tsv")

deep_depth_label <- function(sample_name) {
  without_site <- sub("^[A-Z]{2}\\.", "", sample_name)
  parts <- strsplit(without_site, "\\.")[[1]]
  if (length(parts) == 1) return(parts[1])
  if (length(parts) == 2 && parts[1] == "0") return(paste0("0–", parts[2]))
  paste0(parts[length(parts) - 1], "–", parts[length(parts)])
}
deep_depth_labels <- stats::setNames(
  vapply(samples_deep, deep_depth_label, character(1)), samples_deep
)

df_deep <- prepare_grouped_data(samples_deep) |>
  dplyr::mutate(
    Site = factor(.data$Site, levels = c("BM", "HV", "LM", "NR")),
    Sample = factor(.data$Sample, levels = rev(samples_deep))
  )

p_deep <- ggplot2::ggplot(
  df_deep,
  ggplot2::aes(x = .data$Abundance, y = .data$Sample, fill = .data$Phylum_plot)
) +
  ggplot2::geom_col(
    width = 0.92, colour = "black", linewidth = 0.15,
    orientation = "y", position = ggplot2::position_fill(reverse = TRUE)
  ) +
  ggplot2::facet_wrap(
    facets = ggplot2::vars(Site), nrow = 1, scales = "free_y"
  ) +
  abundance_scale +
  ggplot2::scale_y_discrete(labels = deep_depth_labels, drop = TRUE) +
  fill_scale +
  ggplot2::labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  common_theme +
  ggplot2::theme(
    panel.spacing.x = grid::unit(10, "pt"),
    strip.text.x = ggplot2::element_text(face = "plain", size = 10)
  )

utils::write.table(df_deep, source_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
print(p_deep)
ggplot2::ggsave(output_svg, p_deep, device = svglite::svglite,
                 width = 12.0, height = 2.6, units = "in", bg = "white",
                 system_fonts = list(sans = plot_font))
if (capabilities("cairo")) {
  ggplot2::ggsave(output_pdf, p_deep, device = grDevices::cairo_pdf,
                  width = 12.0, height = 2.6, units = "in", bg = "white")
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure1D_sessionInfo.txt"))
cat("Figure 1D completed.\n")
