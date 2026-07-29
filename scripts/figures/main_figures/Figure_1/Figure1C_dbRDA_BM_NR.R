#!/usr/bin/env Rscript

required_packages <- c(
  "dplyr", "ggplot2", "ggrepel", "readr", "readxl",
  "scales", "svglite", "tibble", "tidyr", "vegan"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript SCRIPT.R <microbial_data.xlsx> <Supplementary_Table_2.xlsx> [output_dir]")
}
microbial_file <- args[[1]]
metadata_file <- args[[2]]
output_dir <- if (length(args) == 3) args[[3]] else "."
if (!file.exists(microbial_file)) stop("Microbial data file not found: ", microbial_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

microbial_data <- readxl::read_excel(microbial_file)
required_microbial_columns <- c("Sample", "OTU", "Abundance")
missing_microbial_columns <- setdiff(required_microbial_columns, names(microbial_data))
if (length(missing_microbial_columns) > 0) {
  stop("Missing microbial columns: ", paste(missing_microbial_columns, collapse = ", "))
}
microbial_data <- microbial_data |>
  dplyr::transmute(
    Sample = trimws(as.character(.data$Sample)),
    OTU = trimws(as.character(.data$OTU)),
    Abundance = suppressWarnings(as.numeric(.data$Abundance))
  )
if (anyNA(microbial_data$Sample) || anyNA(microbial_data$OTU) || anyNA(microbial_data$Abundance)) {
  stop("Missing or non-numeric values were found in the microbial table.")
}
duplicate_microbe_rows <- microbial_data |>
  dplyr::count(.data$Sample, .data$OTU, name = "n") |>
  dplyr::filter(.data$n > 1)
if (nrow(duplicate_microbe_rows) > 0) {
  stop("Duplicated Sample–OTU rows were found in the microbial table.")
}

env_data_raw <- readxl::read_excel(
  metadata_file, sheet = "Supplementary_Table_2", skip = 4, na = ""
)
names(env_data_raw) <- trimws(names(env_data_raw))
required_metadata_columns <- c(
  "Sampling", "Peatland", "Site", "Depth (cm)", "Sample name", "pH",
  "Temperature (°C)", "Electrical conductivity (μS cm⁻¹)",
  "CO₂:CH₄ (mol mol⁻¹)", "NO₃⁻ (µM)", "SO₄²⁻(µM)",
  "total Mn (µM)", "total Fe (µM)", "Propionate (mM)", "Acetate (mM)"
)
missing_metadata_columns <- setdiff(required_metadata_columns, names(env_data_raw))
if (length(missing_metadata_columns) > 0) {
  stop("Missing metadata columns: ", paste(missing_metadata_columns, collapse = ", "))
}

# This reproduces the original dbRDA preparation: a dash in a measured
# chemical variable denotes non-detection/below detection and is entered as 0.
# Blank cells remain NA and are excluded by complete-case filtering.
to_numeric_metadata <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("-", "–", "—")] <- "0"
  x[x == ""] <- NA_character_
  suppressWarnings(as.numeric(x))
}

numeric_source_columns <- c(
  "Temperature (°C)", "pH", "Electrical conductivity (μS cm⁻¹)",
  "CO₂:CH₄ (mol mol⁻¹)", "NO₃⁻ (µM)", "SO₄²⁻(µM)",
  "Propionate (mM)", "Acetate (mM)", "total Mn (µM)", "total Fe (µM)"
)
non_detect_report <- dplyr::bind_rows(lapply(numeric_source_columns, function(column_name) {
  values <- trimws(as.character(env_data_raw[[column_name]]))
  tibble::tibble(
    source_column = column_name,
    number_of_dashes_converted_to_zero = sum(values %in% c("-", "–", "—"), na.rm = TRUE),
    number_of_blank_or_missing_cells = sum(is.na(values) | values == "", na.rm = TRUE)
  )
}))

env_data <- env_data_raw |>
  tidyr::fill(Sampling, Peatland, Site, .direction = "down") |>
  dplyr::filter(.data$Sampling == "Spatial profiling") |>
  dplyr::mutate(
    Group = dplyr::case_when(
      .data$Peatland == "Björsmossen" ~ "BM",
      .data$Peatland == "Norra Romyren" ~ "NR",
      .data$Peatland == "Lungsmossen" ~ "LM",
      .data$Peatland == "Havsjömossen" ~ "HV",
      TRUE ~ NA_character_
    ),
    Site = as.integer(.data$Site),
    Depth_end = as.integer(sub(".*-", "", as.character(.data$`Depth (cm)`))),
    Sample = paste(.data$Group, .data$Site, .data$Depth_end, sep = ".")
  ) |>
  dplyr::transmute(
    Sample, Group,
    Temperature = to_numeric_metadata(.data$`Temperature (°C)`),
    pH = to_numeric_metadata(.data$pH),
    Conductivity = to_numeric_metadata(.data$`Electrical conductivity (μS cm⁻¹)`),
    CO2_CH4 = to_numeric_metadata(.data$`CO₂:CH₄ (mol mol⁻¹)`),
    NO3 = to_numeric_metadata(.data$`NO₃⁻ (µM)`),
    SO4 = to_numeric_metadata(.data$`SO₄²⁻(µM)`),
    Propionate = to_numeric_metadata(.data$`Propionate (mM)`),
    Acetate = to_numeric_metadata(.data$`Acetate (mM)`),
    total_Mn = to_numeric_metadata(.data$`total Mn (µM)`),
    total_Fe = to_numeric_metadata(.data$`total Fe (µM)`)
  )

if (anyNA(env_data$Sample) || anyNA(env_data$Group)) {
  stop("At least one spatial-profile metadata row could not be assigned a sample ID or peatland code.")
}
if (anyDuplicated(env_data$Sample)) {
  stop("Duplicated sample identifiers were generated from Supplementary Table 2.")
}

matching_samples <- intersect(unique(microbial_data$Sample), unique(env_data$Sample))
if (length(matching_samples) == 0) stop("No samples match between microbial data and metadata.")
filtered_microbial <- microbial_data |>
  dplyr::filter(.data$Sample %in% matching_samples, .data$Sample != "BM.3.10")
filtered_metadata <- env_data |>
  dplyr::filter(.data$Sample %in% matching_samples, .data$Sample != "BM.3.10")

variable_labels <- c(
  Temperature = "Temperature", pH = "pH", Conductivity = "Conductivity",
  CO2_CH4 = "CO2:CH4", NO3 = "NO3-", SO4 = "SO4²-",
  Propionate = "Propionate", Acetate = "Acetate",
  total_Mn = "total Mn", total_Fe = "total Fe"
)

run_dbrda <- function(groups, selected_variables, analysis_name, output_prefix,
                      plot_title, group_colors) {
  metadata_selected <- filtered_metadata |>
    dplyr::filter(.data$Group %in% groups)
  microbial_selected <- filtered_microbial |>
    dplyr::filter(.data$Sample %in% metadata_selected$Sample)

  community_wide <- microbial_selected |>
    dplyr::select(.data$Sample, .data$OTU, .data$Abundance) |>
    tidyr::pivot_wider(names_from = OTU, values_from = Abundance,
                       values_fill = list(Abundance = 0)) |>
    dplyr::arrange(.data$Sample)
  community_matrix <- community_wide |>
    dplyr::select(-Sample) |>
    as.matrix()
  rownames(community_matrix) <- community_wide$Sample
  storage.mode(community_matrix) <- "numeric"

  env_complete <- metadata_selected |>
    dplyr::select(.data$Sample, .data$Group, dplyr::all_of(selected_variables)) |>
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(selected_variables)))) |>
    dplyr::arrange(.data$Sample)

  if (nrow(env_complete) <= length(selected_variables) + 1) {
    stop("Too few complete samples for a model containing ",
         length(selected_variables), " environmental variables.")
  }
  missing_community_samples <- setdiff(env_complete$Sample, rownames(community_matrix))
  if (length(missing_community_samples) > 0) {
    stop("Complete metadata samples missing from community matrix: ",
         paste(missing_community_samples, collapse = ", "))
  }

  zero_variance <- selected_variables[
    vapply(env_complete[selected_variables], function(x) stats::sd(x, na.rm = TRUE) == 0,
           FUN.VALUE = logical(1))
  ]
  if (length(zero_variance) > 0) {
    stop("Zero-variance environmental variables: ", paste(zero_variance, collapse = ", "))
  }

  community <- community_matrix[env_complete$Sample, , drop = FALSE]
  env_scaled <- env_complete
  env_scaled[selected_variables] <- scale(env_complete[selected_variables])
  rownames(env_scaled) <- env_scaled$Sample
  if (anyNA(env_scaled[selected_variables])) stop("NA values were generated during scaling.")

  community_hell <- vegan::decostand(community, method = "hellinger")
  env_model <- env_scaled[, selected_variables, drop = FALSE]
  formula_adonis <- stats::reformulate(selected_variables, response = "community_hell")

  permanova_overall <- vegan::adonis2(
    formula_adonis, data = env_model, method = "bray", permutations = 999
  )
  permanova_margin <- vegan::adonis2(
    formula_adonis, data = env_model, method = "bray",
    permutations = 999, by = "margin"
  )

  model <- vegan::capscale(
    community_hell ~ ., data = env_model, distance = "bray"
  )
  model_overall <- vegan::anova.cca(model, permutations = 999)
  model_margin <- vegan::anova.cca(model, permutations = 999, by = "margin")

  constrained_eigenvalues <- vegan::eigenvals(model, model = "constrained")
  if (length(constrained_eigenvalues) < 2) stop("The dbRDA produced fewer than two constrained axes.")
  constrained_percent <- constrained_eigenvalues / sum(constrained_eigenvalues) * 100
  axis1_percent <- round(constrained_percent[1], 2)
  axis2_percent <- round(constrained_percent[2], 2)

  site_scores <- as.data.frame(vegan::scores(model, display = "sites", scaling = 2))
  site_scores$Sample <- rownames(site_scores)
  site_scores <- dplyr::left_join(
    site_scores, env_scaled |> dplyr::select(.data$Sample, .data$Group), by = "Sample"
  )

  fitted_vectors <- vegan::envfit(
    model, env_model, permutations = 999
  )
  vector_scores <- as.data.frame(vegan::scores(fitted_vectors, display = "vectors"))
  vector_scores$Variable <- rownames(vector_scores)
  p_values <- tibble::tibble(
    Variable = names(fitted_vectors$vectors$pvals),
    p_value = as.numeric(fitted_vectors$vectors$pvals)
  ) |>
    dplyr::mutate(
      significance = dplyr::case_when(
        .data$p_value < 0.001 ~ "***",
        .data$p_value < 0.01 ~ "**",
        .data$p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  vector_scores <- vector_scores |>
    dplyr::left_join(p_values, by = "Variable") |>
    dplyr::mutate(
      Variable_label = unname(variable_labels[.data$Variable]),
      Variable_label = dplyr::if_else(is.na(.data$Variable_label),
                                      .data$Variable, .data$Variable_label),
      label = paste0(.data$Variable_label, " (", .data$significance, ")"),
      CAP1_plot = .data$CAP1 * 1.3,
      CAP2_plot = .data$CAP2 * 1.3
    )

  plot <- ggplot2::ggplot(
    site_scores, ggplot2::aes(x = .data$CAP1, y = .data$CAP2, colour = .data$Group)
  ) +
    ggplot2::stat_ellipse(level = 0.95, linewidth = 0.9, alpha = 0.35) +
    ggplot2::geom_point(size = 4.5, alpha = 0.9) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$Sample), size = 4.2, colour = "black",
      box.padding = 0.3, point.padding = 0.2, max.overlaps = 100
    ) +
    ggplot2::geom_segment(
      data = vector_scores,
      ggplot2::aes(x = 0, y = 0, xend = .data$CAP1_plot, yend = .data$CAP2_plot),
      arrow = grid::arrow(length = grid::unit(0.22, "cm"), type = "open"),
      colour = scales::alpha("grey20", 0.65), linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = vector_scores,
      ggplot2::aes(x = .data$CAP1_plot * 1.02, y = .data$CAP2_plot * 1.02,
                   label = .data$label),
      size = 5, colour = "black", box.padding = 0.45,
      point.padding = 0.3, force = 2, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = group_colors) +
    ggplot2::labs(
      title = plot_title,
      x = paste0("dbRDA1 (", axis1_percent, "% of constrained variation)"),
      y = paste0("dbRDA2 (", axis2_percent, "% of constrained variation)"),
      colour = "Peatland"
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.title = ggplot2::element_text(size = 18),
      axis.text = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 17),
      legend.text = ggplot2::element_text(size = 15)
    )

  svg_file <- file.path(output_dir, paste0(output_prefix, ".svg"))
  pdf_file <- file.path(output_dir, paste0(output_prefix, ".pdf"))
  print(plot)
  ggplot2::ggsave(svg_file, plot, device = svglite::svglite,
                   width = 8, height = 8, units = "in", bg = "white")
  if (capabilities("cairo")) {
    ggplot2::ggsave(pdf_file, plot, device = grDevices::cairo_pdf,
                    width = 8, height = 8, units = "in", bg = "white")
  }

  readr::write_tsv(env_complete,
    file.path(output_dir, paste0(output_prefix, "_samples_and_environment.tsv")))
  readr::write_tsv(site_scores,
    file.path(output_dir, paste0(output_prefix, "_site_scores.tsv")))
  readr::write_tsv(vector_scores,
    file.path(output_dir, paste0(output_prefix, "_envfit_vectors.tsv")))
  writeLines(capture.output(permanova_overall),
    file.path(output_dir, paste0(output_prefix, "_PERMANOVA_overall.txt")))
  writeLines(capture.output(permanova_margin),
    file.path(output_dir, paste0(output_prefix, "_PERMANOVA_marginal.txt")))
  writeLines(capture.output(model_overall),
    file.path(output_dir, paste0(output_prefix, "_dbRDA_overall_test.txt")))
  writeLines(capture.output(model_margin),
    file.path(output_dir, paste0(output_prefix, "_dbRDA_marginal_tests.txt")))
  writeLines(capture.output(summary(model)),
    file.path(output_dir, paste0(output_prefix, "_model_summary.txt")))

  cat("\n", analysis_name, "\n", sep = "")
  cat("Samples used: ", nrow(env_complete), "\n", sep = "")
  cat("Environmental variables: ", paste(selected_variables, collapse = ", "), "\n", sep = "")
  cat("Saved: ", svg_file, "\n", sep = "")
  invisible(list(plot = plot, model = model, samples = env_complete))
}


###############################################################################
# FIGURE 1C — dbRDA, BJÖRSMOSSEN + NORRA ROMYREN, FULL ENVIRONMENTAL DATASET
###############################################################################

result <- run_dbrda(
  groups = c("BM", "NR"),
  selected_variables = c(
    "Temperature", "pH", "Conductivity", "CO2_CH4", "NO3", "SO4",
    "Propionate", "Acetate", "total_Mn", "total_Fe"
  ),
  analysis_name = "Figure 1C dbRDA: BM + NR, full environmental dataset",
  output_prefix = "Figure1C_dbRDA_BM_NR",
  plot_title = "dbRDA (BM and NR — full environmental variables)",
  group_colors = c("BM" = "#3366CC", "NR" = "#E6863B")
)

readr::write_tsv(non_detect_report,
  file.path(output_dir, "Figure1C_non_detect_conversion_report.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure1C_sessionInfo.txt"))
