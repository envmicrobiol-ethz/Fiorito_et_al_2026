#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Extended Data Figure 2A: a taxonomic-composition comparison
# across 16S rRNA amplicons, SingleM profiles of metagenomic reads, and
# reconstructed MAGs.
#
# Eight phyla are displayed separately. All remaining phyla are grouped as
# "Other phyla", while the fraction of metagenomic reads not represented by
# reconstructed MAGs is shown as "Unmapped" in the MAG panel.
#
# INPUT
# 1. Excel table containing 16S relative-abundance data.
# 2. Directory containing SingleM *_profile.tsv files.
# 3. MAG taxonomy and abundance TSV.
# 4. Output directory.
#
# OUTPUT
# Extended Data Figure 2A as SVG and PDF, plus the plotted values as TSV.
#
# USAGE
# Rscript ED_Figure_2A.R \
#   combined_16S_taxonomy_relative_abundance.xlsx \
#   singlem_profiles_directory \
#   MAG_taxonomy_abundance.tsv \
#   output_directory

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(readr)
  library(stringr)
  library(svglite)
  library(tibble)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript ED_Figure_2A.R",
      "<16S_taxonomy.xlsx>",
      "<SingleM_directory>",
      "<MAG_taxonomy_abundance.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

excel_file_path <- args[[1]]
singlem_dir <- args[[2]]
mags_file <- args[[3]]
output_dir <- args[[4]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "ED_Figure_2A_taxonomy_16S_SingleM_MAGs.svg"
)

out_pdf <- file.path(
  output_dir,
  "ED_Figure_2A_taxonomy_16S_SingleM_MAGs.pdf"
)

out_tsv <- file.path(
  output_dir,
  "ED_Figure_2A_taxonomy_16S_SingleM_MAGs_values.tsv"
)

# ============================================================
# 2) SAMPLE ORDER AND MATCHING
#
# BM.3.10 is included in the common x-axis order because SingleM and MAG data
# are available for it. It is removed only from the 16S data, so the 16S panel
# retains the corresponding empty bar position.
# The two vectors must correspond one-to-one by position.
# ============================================================

custom_order_excel <- c(
  "BM.3.10", "BM.3.30", "BM.3.60", "BM.0.50", "BM.100", "BM.200",
  "BM.250", "BM.300", "BM.490.540",
  
  "HV.0.25", "HV.100", "HV.200", "HV.265.280",
  
  "LM.3.10", "LM.3.30", "LM.3.60", "LM.0.25", "LM.100",
  "LM.200", "LM.250", "LM.300", "LM.550", "LM.650",
  "LM.675.700",
  
  "NR.3.10", "NR.3.30", "NR.3.60", "NR.0.25", "NR.100",
  "NR.200", "NR.250", "NR.325.350"
)

custom_order_singleM <- c(
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


map_samples <- tibble::tibble(
  Sample16S = custom_order_excel,
  SampleMG = custom_order_singleM
) |>
  dplyr::mutate(
    Sample = factor(.data$Sample16S, levels = custom_order_excel),
    
    Peatland = dplyr::case_when(
      stringr::str_starts(.data$SampleMG, "Bjorsmossen") ~ "BM",
      stringr::str_starts(.data$SampleMG, "Havsjomossen") ~ "HV",
      stringr::str_starts(.data$SampleMG, "Lungsmossen") ~ "LM",
      stringr::str_starts(.data$SampleMG, "NorraRomyren") ~ "NR",
      TRUE ~ NA_character_
    ),
    
    Depth = dplyr::case_when(
      stringr::str_detect(.data$SampleMG, "site3_10cm") ~ "10",
      stringr::str_detect(.data$SampleMG, "site3_30cm") ~ "30",
      stringr::str_detect(.data$SampleMG, "site3_60cm") ~ "60",
      TRUE ~ stringr::str_extract(.data$SampleMG, "\\d{1,3}-\\d{1,3}(?=cm)")
    ),
    Depth = stringr::str_replace_all(.data$Depth, "-", "–"),
    
    Peatland = factor(.data$Peatland, levels = c("BM", "HV", "LM", "NR"))
  )


# Named lookups avoid fragile joins/select calls involving SampleMG.
samplemg_to_16s <- stats::setNames(
  as.character(map_samples$Sample16S),
  as.character(map_samples$SampleMG)
)

sample_to_peatland <- stats::setNames(
  as.character(map_samples$Peatland),
  as.character(map_samples$Sample16S)
)

sample_to_depth <- stats::setNames(
  as.character(map_samples$Depth),
  as.character(map_samples$Sample16S)
)

# ============================================================
# 3) TAXONOMY HARMONIZATION
# ============================================================

ensure_p_prefix <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- "Unclassified"
  x <- stringr::str_trim(x)
  
  out <- ifelse(
    stringr::str_detect(x, "p__"),
    stringr::str_extract(x, "p__[^;\\s]+"),
    NA_character_
  )
  
  out <- ifelse(
    is.na(out),
    ifelse(
      x %in% c("Unclassified", "Unassigned", "p", "p_", "p__"),
      "Unclassified",
      paste0("p__", x)
    ),
    out
  )
  
  out[out %in% c("p", "p_", "p__")] <- "Unclassified"
  out
}

phylum_map_16S_to_MAG <- c(
  "p__Actinobacteriota" = "p__Actinomycetota",
  "p__Chloroflexi" = "p__Chloroflexota",
  "p__Firmicutes" = "p__Bacillota",
  "p__Proteobacteria" = "p__Pseudomonadota",
  "p__Crenarchaeota" = "p__Thermoproteota",
  "p__Euryarchaeota" = "p__Methanobacteriota",
  "p__Halobacterota" = "p__Halobacteriota"
)

to_MAG_phylum <- function(ph) {
  ph <- ensure_p_prefix(ph)
  ph <- stringr::str_replace(ph, "(^p__[^_]+)_[^_]+$", "\\1")
  mapped <- unname(phylum_map_16S_to_MAG[ph])
  ifelse(!is.na(mapped), mapped, ph)
}

legend_labels <- function(x) {
  x |>
    stringr::str_remove("^p__")
}

# ============================================================
# 4) EIGHT DISPLAYED PHYLA
#
# Seven phyla match the revised phylogenetic tree; Pseudomonadota is added here
# Archaeal phyla are adjacent, followed by bacterial phyla.
# Pseudomonadota is shown explicitly because it is recurrent across the
# three profiling approaches and remains within the reviewer's suggested range of 7–10 coloured groups.
# ============================================================

archaea_top8 <- c(
  "p__Halobacteriota",
  "p__Thermoproteota"
)

bacteria_top8 <- c(
  "p__Acidobacteriota",
  "p__Actinomycetota",
  "p__Verrucomicrobiota",
  "p__Chloroflexota",
  "p__Desulfobacterota",
  "p__Pseudomonadota"
)

top8_phyla <- c(
  archaea_top8,
  bacteria_top8
)

# Stack order, bottom to top.
stack_order <- c(
  "Other phyla",
  archaea_top8,
  bacteria_top8,
  "Unmapped"
)

# Legend order, grouped by domain.
legend_order <- c(
  archaea_top8,
  bacteria_top8,
  "Other phyla",
  "Unmapped"
)

phylum_colors <- c(
  "p__Halobacteriota" = "#2171B5",
  "p__Thermoproteota" = "#08306B",
  "p__Acidobacteriota" = "#B2182B",
  "p__Actinomycetota" = "#FFF04F",
  "p__Verrucomicrobiota" = "#D8B365",
  "p__Chloroflexota" = "#FF8C00",
  "p__Desulfobacterota" = "#C77526",
  "p__Pseudomonadota" = "#8C510A",
  "Other phyla" = "grey60",
  "Unmapped" = "#E6E6E6"
)

# ============================================================
# 5) INPUT VALIDATION
# ============================================================

if (!file.exists(excel_file_path)) {
  stop("16S Excel file not found: ", excel_file_path)
}

if (!dir.exists(singlem_dir)) {
  stop("SingleM directory not found: ", singlem_dir)
}

if (!file.exists(mags_file)) {
  stop("MAG abundance table not found: ", mags_file)
}

# ============================================================
# 6) 16S DATA
# ============================================================

# BM.3.10 is excluded only from the 16S dataset.
# Its x-axis position is retained through the SingleM and MAG panels.
exclude_16s_samples <- c("BM.3.10")

df_16s <- readxl::read_excel(excel_file_path) |>
  dplyr::mutate(
    Sample = as.character(.data$Sample),
    Phylum = to_MAG_phylum(.data$Phylum),
    Abundance = as.numeric(.data$Abundance)
  ) |>
  dplyr::filter(
    .data$Sample %in% custom_order_excel,
    !.data$Sample %in% exclude_16s_samples
  ) |>
  dplyr::group_by(.data$Sample, .data$Phylum) |>
  dplyr::summarise(
    Abundance = sum(.data$Abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(.data$Sample) |>
  dplyr::mutate(
    sample_total = sum(.data$Abundance, na.rm = TRUE),
    rel_abund = dplyr::if_else(
      .data$sample_total > 0,
      .data$Abundance / .data$sample_total,
      0
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    Sample = .data$Sample,
    Phylum = .data$Phylum,
    rel_abund = .data$rel_abund,
    source = "16S"
  )


# ============================================================
# 7) SingleM RAW-READ DATA
#
# Robust handling of BM.3.10:
#   1. Match expected sample names after removing common read/profile suffixes.
#   2. If the sample column does not match exactly, search the file name too.
#   3. Prefer rows whose taxonomy terminates at phylum level, as in the
#      original script. If a file has no such rows, fall back to every row
#      containing a phylum assignment and collapse coverage by phylum.
# ============================================================

normalize_singlem_sample_name <- function(x) {
  x <- basename(as.character(x))
  x <- stringr::str_trim(x)
  
  # File-name suffixes.
  x <- stringr::str_remove(x, "\\.tsv$")
  x <- stringr::str_remove(x, "_profile$")
  
  # Common paired-end/read suffixes.
  x <- stringr::str_remove(x, "_PE[0-9]+$")
  x <- stringr::str_remove(x, "\\.PE[0-9]+$")
  x <- stringr::str_remove(x, "_R[12](?:_001)?$")
  
  x
}

resolve_singlem_sample <- function(raw_samples, file) {
  
  candidates <- unique(c(
    normalize_singlem_sample_name(raw_samples),
    normalize_singlem_sample_name(basename(file))
  ))
  
  candidates <- candidates[
    !is.na(candidates) & candidates != ""
  ]
  
  # First choice: exact match after normalization.
  exact_hits <- intersect(candidates, custom_order_singleM)
  
  if (length(exact_hits) == 1) {
    return(list(
      SampleMG = exact_hits,
      match_mode = "exact normalized sample/file name",
      candidates = paste(candidates, collapse = " | ")
    ))
  }
  
  # Second choice: one expected sample name is contained in a longer
  # sample/file string. This handles suffixes not anticipated above.
  substring_hits <- custom_order_singleM[
    vapply(
      custom_order_singleM,
      function(expected) {
        any(vapply(
          candidates,
          function(candidate) {
            stringr::str_detect(candidate, stringr::fixed(expected)) ||
              stringr::str_detect(expected, stringr::fixed(candidate))
          },
          logical(1)
        ))
      },
      logical(1)
    )
  ]
  
  substring_hits <- unique(substring_hits)
  
  if (length(substring_hits) == 1) {
    return(list(
      SampleMG = substring_hits,
      match_mode = "substring sample/file name",
      candidates = paste(candidates, collapse = " | ")
    ))
  }
  
  stop(
    paste0(
      "Could not uniquely match SingleM file to an expected sample.\n",
      "File: ", file, "\n",
      "Observed sample/file names: ", paste(candidates, collapse = " | "), "\n",
      "Expected BM.3.10 metagenomic name: ", custom_order_singleM[1]
    )
  )
}

process_singlem_file <- function(file) {
  
  d <- readr::read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  required <- c("taxonomy", "coverage")
  missing <- setdiff(required, names(d))
  
  if (length(missing) > 0) {
    stop(
      "SingleM file ", basename(file),
      " lacks columns: ", paste(missing, collapse = ", ")
    )
  }
  
  raw_samples <- if ("sample" %in% names(d)) {
    unique(as.character(d$sample))
  } else {
    character(0)
  }
  
  sample_resolution <- resolve_singlem_sample(
    raw_samples = raw_samples,
    file = file
  )
  
  taxonomy_text <- as.character(d$taxonomy)
  
  # Original behaviour: use phylum-level rows when they exist.
  terminal_phylum <- stringr::str_detect(
    taxonomy_text,
    "(^|;)\\s*p__[^;\\s]+\\s*$"
  )
  
  # Robust fallback: accept deeper lineages and extract their phylum.
  contains_phylum <- stringr::str_detect(
    taxonomy_text,
    "(^|;)\\s*p__[^;\\s]+"
  )
  
  if (any(terminal_phylum, na.rm = TRUE)) {
    keep_rows <- terminal_phylum
    taxonomy_filter_mode <- "terminal phylum rows"
  } else {
    keep_rows <- contains_phylum
    taxonomy_filter_mode <- "fallback: all rows containing phylum"
  }
  
  if (!any(keep_rows, na.rm = TRUE)) {
    warning(
      "No rows containing a phylum assignment in SingleM file: ",
      basename(file)
    )
    
    return(tibble::tibble(
      SampleMG = character(0),
      Phylum = character(0),
      Abundance = numeric(0),
      Source_file = character(0),
      Raw_sample_names = character(0),
      Sample_match_mode = character(0),
      Taxonomy_filter_mode = character(0)
    ))
  }
  
  d |>
    dplyr::filter(.env$keep_rows) |>
    dplyr::mutate(
      SampleMG = sample_resolution$SampleMG,
      Phylum = stringr::str_extract(
        as.character(.data$taxonomy),
        "p__[^;\\s]+"
      ),
      Phylum = to_MAG_phylum(.data$Phylum),
      Abundance = suppressWarnings(as.numeric(.data$coverage)),
      Source_file = basename(file),
      Raw_sample_names = if (length(raw_samples) > 0) {
        paste(unique(raw_samples), collapse = " | ")
      } else {
        "sample column absent; matched from file name"
      },
      Sample_match_mode = sample_resolution$match_mode,
      Taxonomy_filter_mode = taxonomy_filter_mode
    ) |>
    dplyr::filter(
      !is.na(.data$Phylum),
      !is.na(.data$Abundance)
    ) |>
    dplyr::group_by(
      .data$SampleMG,
      .data$Phylum,
      .data$Source_file,
      .data$Raw_sample_names,
      .data$Sample_match_mode,
      .data$Taxonomy_filter_mode
    ) |>
    dplyr::summarise(
      Abundance = sum(.data$Abundance, na.rm = TRUE),
      .groups = "drop"
    )
}

singlem_files <- list.files(
  path = singlem_dir,
  pattern = "_profile\\.tsv$",
  full.names = TRUE
)

if (length(singlem_files) == 0) {
  stop("No SingleM *_profile.tsv files found in: ", singlem_dir)
}

singlem_processed <- lapply(
  singlem_files,
  process_singlem_file
) |>
  dplyr::bind_rows()

df_singlem <- singlem_processed |>
  dplyr::select(
    dplyr::all_of(c("SampleMG", "Phylum", "Abundance"))
  ) |>
  dplyr::group_by(
    .data$SampleMG,
    .data$Phylum
  ) |>
  dplyr::summarise(
    Abundance = sum(.data$Abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Sample16S = unname(
      samplemg_to_16s[as.character(.data$SampleMG)]
    )
  ) |>
  dplyr::filter(!is.na(.data$Sample16S)) |>
  dplyr::group_by(
    .data$Sample16S,
    .data$Phylum
  ) |>
  dplyr::summarise(
    Abundance = sum(.data$Abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(.data$Sample16S) |>
  dplyr::mutate(
    sample_total = sum(.data$Abundance, na.rm = TRUE),
    rel_abund = dplyr::if_else(
      .data$sample_total > 0,
      .data$Abundance / .data$sample_total,
      0
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    Sample = .data$Sample16S,
    Phylum = .data$Phylum,
    rel_abund = .data$rel_abund,
    source = "SingleM reads"
  )

# ============================================================
# 8) MAG DATA AND UNMAPPED FRACTION
# ============================================================

mag_wide <- readr::read_tsv(
  mags_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  dplyr::mutate(
    Phylum = to_MAG_phylum(.data$Phylum)
  )

required_mag_columns <- c("genome", "Phylum")
missing_mag_columns <- setdiff(required_mag_columns, names(mag_wide))

if (length(missing_mag_columns) > 0) {
  stop(
    "MAG table lacks columns: ",
    paste(missing_mag_columns, collapse = ", ")
  )
}

sample_cols_mag <- intersect(
  names(mag_wide),
  custom_order_singleM
)

if (length(sample_cols_mag) != length(custom_order_singleM)) {
  missing_samples <- setdiff(custom_order_singleM, sample_cols_mag)
  stop(
    "MAG table is missing sample columns: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Fraction of metagenomic reads represented by reconstructed MAGs.
mapped_from_mag <- mag_wide |>
  dplyr::select(dplyr::all_of(sample_cols_mag)) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::everything(),
      ~ sum(as.numeric(.x), na.rm = TRUE)
    )
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "SampleMG",
    values_to = "mapped_frac"
  ) |>
  dplyr::mutate(
    mapped_frac = as.numeric(.data$mapped_frac),
    mapped_frac = dplyr::if_else(
      is.na(.data$mapped_frac),
      0,
      .data$mapped_frac
    ),
    mapped_frac = dplyr::if_else(
      .data$mapped_frac > 1.5,
      .data$mapped_frac / 100,
      .data$mapped_frac
    ),
    mapped_frac = pmin(pmax(.data$mapped_frac, 0), 1)
  ) |>
  dplyr::mutate(
    Sample = unname(
      samplemg_to_16s[as.character(.data$SampleMG)]
    )
  ) |>
  dplyr::filter(!is.na(.data$Sample)) |>
  dplyr::transmute(
    Sample = .data$Sample,
    mapped_frac = .data$mapped_frac
  )


# Composition within the reconstructed MAG fraction.
df_mags <- mag_wide |>
  dplyr::select(
    .data$genome,
    .data$Phylum,
    dplyr::all_of(sample_cols_mag)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols_mag),
    names_to = "SampleMG",
    values_to = "Abundance"
  ) |>
  dplyr::mutate(
    Abundance = as.numeric(.data$Abundance)
  ) |>
  dplyr::filter(
    !is.na(.data$Abundance),
    .data$Abundance > 0
  ) |>
  dplyr::mutate(
    Sample16S = unname(
      samplemg_to_16s[as.character(.data$SampleMG)]
    )
  ) |>
  dplyr::filter(!is.na(.data$Sample16S)) |>
  dplyr::group_by(.data$Sample16S, .data$Phylum) |>
  dplyr::summarise(
    Abundance = sum(.data$Abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(.data$Sample16S) |>
  dplyr::mutate(
    sample_total = sum(.data$Abundance, na.rm = TRUE),
    rel_abund = dplyr::if_else(
      .data$sample_total > 0,
      .data$Abundance / .data$sample_total,
      0
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    Sample = .data$Sample16S,
    Phylum = .data$Phylum,
    rel_abund = .data$rel_abund,
    source = "MAGs"
  )

# ============================================================
# 9) COLLAPSE TO 8 PHYLA + OTHER
#
# IMPORTANT:
#   - all non-selected phyla are retained in "Other phyla"
#   - negative/non-finite values are set to zero
#   - 16S and SingleM are normalized AGAIN after collapsing
#   - MAG taxonomic composition is normalized within the mapped fraction
# This guarantees that no unexplained white residual remains.
# ============================================================

collapse_phyla <- function(data) {
  data |>
    dplyr::mutate(
      rel_abund = suppressWarnings(as.numeric(.data$rel_abund)),
      rel_abund = dplyr::if_else(
        is.finite(.data$rel_abund) & .data$rel_abund > 0,
        .data$rel_abund,
        0
      ),
      Phylum_plot = dplyr::if_else(
        .data$Phylum %in% top8_phyla,
        .data$Phylum,
        "Other phyla"
      )
    ) |>
    dplyr::group_by(
      .data$source,
      .data$Sample,
      .data$Phylum_plot
    ) |>
    dplyr::summarise(
      rel_abund = sum(.data$rel_abund, na.rm = TRUE),
      .groups = "drop"
    )
}

normalize_full_bar <- function(data, dataset_label) {
  
  totals_before <- data |>
    dplyr::group_by(.data$source, .data$Sample) |>
    dplyr::summarise(
      total_before = sum(.data$rel_abund, na.rm = TRUE),
      .groups = "drop"
    )
  
  zero_samples <- totals_before |>
    dplyr::filter(
      !is.finite(.data$total_before) |
        .data$total_before <= 0
    )
  
  if (nrow(zero_samples) > 0) {
    stop(
      dataset_label,
      " contains sample(s) with no positive abundance after processing: ",
      paste(zero_samples$Sample, collapse = ", ")
    )
  }
  
  data |>
    dplyr::left_join(
      totals_before,
      by = c("source", "Sample")
    ) |>
    dplyr::mutate(
      rel_abund = .data$rel_abund / .data$total_before
    ) |>
    dplyr::select(-dplyr::all_of("total_before"))
}

# 16S and SingleM are compositional profiles and must each fill 100%.
df_16s_collapsed <- df_16s |>
  collapse_phyla() |>
  normalize_full_bar("16S")

df_singlem_collapsed <- df_singlem |>
  collapse_phyla() |>
  normalize_full_bar("SingleM reads")

# MAG composition is first normalized to 100% within the reconstructed MAGs.
df_mags_composition <- df_mags |>
  collapse_phyla() |>
  normalize_full_bar("MAG composition")

# ============================================================
# 10) SCALE MAG COMPOSITION AND ADD UNMAPPED
# ============================================================

# Scale the normalized MAG composition to the fraction represented by MAGs.
df_mags_scaled <- df_mags_composition |>
  dplyr::left_join(
    mapped_from_mag,
    by = "Sample"
  ) |>
  dplyr::mutate(
    mapped_frac = dplyr::if_else(
      is.finite(.data$mapped_frac),
      pmin(pmax(.data$mapped_frac, 0), 1),
      0
    ),
    rel_abund = .data$rel_abund * .data$mapped_frac
  ) |>
  dplyr::select(
    .data$source,
    .data$Sample,
    .data$Phylum_plot,
    .data$rel_abund
  )

# Unmapped is created ONLY for the MAG panel.
df_unmapped <- mapped_from_mag |>
  dplyr::mutate(
    mapped_frac = dplyr::if_else(
      is.finite(.data$mapped_frac),
      pmin(pmax(.data$mapped_frac, 0), 1),
      0
    )
  ) |>
  dplyr::transmute(
    source = "MAGs",
    Sample = .data$Sample,
    Phylum_plot = "Unmapped",
    rel_abund = 1 - .data$mapped_frac
  )

# ============================================================
# 11) FINAL COMPLETE PLOTTING TABLE
#
# An explicit source × sample × category grid is used so that:
#   - every sample position is aligned across panels;
#   - BM.3.10 remains blank ONLY in 16S;
#   - missing categories are zero, not missing values;
#   - every non-empty 16S/SingleM bar is exactly 1;
#   - every MAG bar, including Unmapped, is exactly 1.
# ============================================================

observed_values <- dplyr::bind_rows(
  df_16s_collapsed,
  df_singlem_collapsed,
  df_mags_scaled,
  df_unmapped
) |>
  dplyr::group_by(
    .data$source,
    .data$Sample,
    .data$Phylum_plot
  ) |>
  dplyr::summarise(
    rel_abund = sum(.data$rel_abund, na.rm = TRUE),
    .groups = "drop"
  )

source_sample_grid <- dplyr::bind_rows(
  tibble::tibble(
    source = "16S",
    Sample = custom_order_excel
  ),
  tibble::tibble(
    source = "SingleM reads",
    Sample = custom_order_excel
  ),
  tibble::tibble(
    source = "MAGs",
    Sample = custom_order_excel
  )
)

complete_grid <- tidyr::crossing(
  source_sample_grid,
  Phylum_plot = stack_order
)

df_plot <- complete_grid |>
  dplyr::left_join(
    observed_values,
    by = c("source", "Sample", "Phylum_plot")
  ) |>
  dplyr::mutate(
    rel_abund = dplyr::coalesce(.data$rel_abund, 0),
    
    # Unmapped cannot exist outside MAGs.
    rel_abund = dplyr::if_else(
      .data$source != "MAGs" & .data$Phylum_plot == "Unmapped",
      0,
      .data$rel_abund
    ),
    
    # BM.3.10 must be blank only in 16S.
    rel_abund = dplyr::if_else(
      .data$source == "16S" & .data$Sample == "BM.3.10",
      0,
      .data$rel_abund
    )
  )

# Final hard normalization for non-empty 16S and SingleM bars.
nonmag_totals <- df_plot |>
  dplyr::filter(.data$source %in% c("16S", "SingleM reads")) |>
  dplyr::group_by(.data$source, .data$Sample) |>
  dplyr::summarise(
    total = sum(.data$rel_abund, na.rm = TRUE),
    .groups = "drop"
  )

df_plot <- df_plot |>
  dplyr::left_join(
    nonmag_totals |>
      dplyr::rename(nonmag_total = total),
    by = c("source", "Sample")
  ) |>
  dplyr::mutate(
    rel_abund = dplyr::if_else(
      .data$source %in% c("16S", "SingleM reads") &
        .data$nonmag_total > 0,
      .data$rel_abund / .data$nonmag_total,
      .data$rel_abund
    )
  ) |>
  dplyr::select(-dplyr::all_of("nonmag_total"))

# MAG bars should already equal one, but explicitly correct only tiny numerical
# residuals by assigning them to Unmapped. This does not alter the mapped
# taxonomic composition.
mag_totals_before_residual <- df_plot |>
  dplyr::filter(.data$source == "MAGs") |>
  dplyr::group_by(.data$Sample) |>
  dplyr::summarise(
    total = sum(.data$rel_abund, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    residual = 1 - .data$total
  )

df_plot <- df_plot |>
  dplyr::left_join(
    mag_totals_before_residual |>
      dplyr::select(.data$Sample, .data$residual),
    by = "Sample"
  ) |>
  dplyr::mutate(
    rel_abund = dplyr::if_else(
      .data$source == "MAGs" & .data$Phylum_plot == "Unmapped",
      pmax(.data$rel_abund + dplyr::coalesce(.data$residual, 0), 0),
      .data$rel_abund
    )
  ) |>
  dplyr::select(-dplyr::all_of("residual"))

# Add facet metadata and final factors.
df_plot <- df_plot |>
  dplyr::mutate(
    Peatland = unname(
      sample_to_peatland[as.character(.data$Sample)]
    ),
    Depth = unname(
      sample_to_depth[as.character(.data$Sample)]
    ),
    source = factor(
      .data$source,
      levels = c("16S", "SingleM reads", "MAGs")
    ),
    Sample = factor(
      .data$Sample,
      levels = custom_order_excel
    ),
    Peatland = factor(
      .data$Peatland,
      levels = c("BM", "HV", "LM", "NR")
    ),
    Phylum_plot = factor(
      .data$Phylum_plot,
      levels = stack_order
    )
  )


readr::write_tsv(
  df_plot |>
    dplyr::mutate(
      source = as.character(.data$source),
      Sample = as.character(.data$Sample),
      Peatland = as.character(.data$Peatland),
      Phylum_plot = as.character(.data$Phylum_plot)
    ),
  out_tsv
)

# ============================================================
# 12) PLOT
# ============================================================

p_tax <- ggplot2::ggplot(
  df_plot,
  ggplot2::aes(
    x = .data$Sample,
    y = .data$rel_abund,
    fill = .data$Phylum_plot
  )
) +
  ggplot2::geom_col(
    width = 0.92,
    colour = "grey25",
    linewidth = 0.18,
    position = "stack"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$source),
    cols = ggplot2::vars(.data$Peatland),
    scales = "free_x",
    space = "free_x",
    switch = "y",
    drop = TRUE
  ) +
  ggplot2::scale_x_discrete(
    labels = stats::setNames(
      map_samples$Depth,
      map_samples$Sample16S
    ),
    drop = TRUE
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) x * 100,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0, 1),
    clip = "on"
  ) +
  ggplot2::scale_fill_manual(
    values = phylum_colors,
    breaks = legend_order,
    labels = c(
      "p__Halobacteriota" = "Halobacteriota",
      "p__Thermoproteota" = "Thermoproteota",
      "p__Acidobacteriota" = "Acidobacteriota",
      "p__Actinomycetota" = "Actinomycetota",
      "p__Verrucomicrobiota" = "Verrucomicrobiota",
      "p__Chloroflexota" = "Chloroflexota",
      "p__Desulfobacterota" = "Desulfobacterota",
      "p__Pseudomonadota" = "Pseudomonadota",
      "Other phyla" = "Other phyla",
      "Unmapped" = "Unmapped (MAGs only)"
    ),
    drop = FALSE,
    na.value = "#FF00FF",
    name = "Phylum"
  ) +
  ggplot2::labs(
    x = "Depth (cm)",
    y = "Relative abundance (%)"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(
        colour = "grey25",
        linewidth = 0.2
      )
    )
  ) +
  ggplot2::theme_classic(
    base_size = 12,
    base_family = "Helvetica"
  ) +
  ggplot2::theme(
    panel.border = ggplot2::element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.75
    ),
    
    panel.spacing.x = grid::unit(1.5, "mm"),
    panel.spacing.y = grid::unit(1.2, "mm"),
    
    strip.background.x = ggplot2::element_rect(
      fill = "grey94",
      colour = "black",
      linewidth = 0.75
    ),
    
    strip.background.y = ggplot2::element_blank(),
    
    strip.text.x = ggplot2::element_text(
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    strip.text.y.left = ggplot2::element_text(
      face = "bold",
      size = 11,
      angle = 90,
      colour = "black"
    ),
    
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8.5,
      colour = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      colour = "black"
    ),
    
    axis.title = ggplot2::element_text(
      colour = "black"
    ),
    
    legend.position = "right",
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    legend.text = ggplot2::element_text(
      size = 10
    ),
    legend.key.height = grid::unit(0.42, "cm"),
    legend.key.width = grid::unit(0.42, "cm"),
    
    plot.margin = ggplot2::margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8,
      unit = "pt"
    )
  )

print(p_tax)

# ============================================================
# 13) SAVE
# ============================================================

ggplot2::ggsave(
  filename = out_svg,
  plot = p_tax,
  device = svglite::svglite,
  width = 12.0,
  height = 7.4,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

if (capabilities("cairo")) {
  ggplot2::ggsave(
    filename = out_pdf,
    plot = p_tax,
    device = grDevices::cairo_pdf,
    width = 12.0,
    height = 7.4,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
} else {
  warning("PDF not saved because Cairo support is unavailable.")
}
