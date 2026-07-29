#!/usr/bin/env Rscript

###############################################################################
# 16S rRNA ALPHA DIVERSITY ANALYSIS
#
# This script:
#   1. Loads the experimental design, OTU table and taxonomy
#   2. Removes chloroplast and mitochondrial OTUs and retains
#      OTUs assigned to Bacteria or Archaea
#   3. Rarefies the OTU table to the minimum sequencing depth
#      across 100 independent iterations
#   4. Calculates observed richness and Shannon diversity
#      for each rarefied table
#   5. Averages alpha-diversity metrics across rarefactions
#   6. Calculates Pielou's evenness
#   7. Applies bestNormalize transformations
#   8. Performs two-way ANOVA using treatment, system and their interaction
#
# Expected input files:
#   design.txt
#   9.all.OTU_map.txt
#   9.all.OTU_tax.sintax.rf.txt
###############################################################################

###############################################################################
# 1. PACKAGES
###############################################################################

library(bestNormalize)
library(dplyr)
library(gtools)
library(parallel)
library(readr)
library(tibble)
library(vegan)

###############################################################################
# 2. INPUT AND OUTPUT PATHS
###############################################################################

# Run the script from the directory containing the input files,
# or change analysis_dir accordingly.
analysis_dir <- "."

design_file <- file.path(
  analysis_dir,
  "design.txt"
)

otu_table_file <- file.path(
  analysis_dir,
  "9.all.OTU_map.txt"
)

taxonomy_file <- file.path(
  analysis_dir,
  "9.all.OTU_tax.sintax.rf.txt"
)

output_dir <- file.path(
  analysis_dir,
  "alpha_diversity_results"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# 3. ANALYSIS PARAMETERS
###############################################################################

iterations <- 100

threshold_rare <- 0
threshold_sparse <- 1
threshold_outlier <- 0

###############################################################################
# 4. CHECK INPUT FILES
###############################################################################

input_files <- c(
  design_file,
  otu_table_file,
  taxonomy_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "The following input files were not found:",
      paste(missing_files, collapse = "\n"),
      sep = "\n"
    )
  )
}

###############################################################################
# 5. LOAD EXPERIMENTAL DESIGN
###############################################################################

design <- readr::read_tsv(
  design_file,
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    dplyr::across(
      where(is.character),
      trimws
    )
  )

# Accept either "sample" or "Sample" as the sample-column name.
if ("sample" %in% names(design)) {
  names(design)[names(design) == "sample"] <- "Sample"
}

required_design_columns <- c(
  "Sample",
  "date",
  "block",
  "treatment",
  "system"
)

missing_design_columns <- setdiff(
  required_design_columns,
  names(design)
)

if (length(missing_design_columns) > 0) {
  stop(
    paste(
      "Missing columns in design.txt:",
      paste(missing_design_columns, collapse = ", ")
    )
  )
}

design <- design %>%
  dplyr::mutate(
    date = as.factor(.data$date),
    block = as.factor(.data$block),
    treatment = as.factor(.data$treatment),
    system = as.factor(.data$system),
    TxS = interaction(
      .data$treatment,
      .data$system,
      sep = " x ",
      drop = TRUE
    )
  ) %>%
  droplevels()

design <- design[
  gtools::mixedorder(design$Sample),
  ,
  drop = FALSE
]

###############################################################################
# 6. LOAD AND FILTER TAXONOMY
###############################################################################

tax <- read.table(
  file = taxonomy_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)

tax <- tax[
  gtools::mixedorder(rownames(tax)),
  ,
  drop = FALSE
]

# Replace empty taxonomy fields with "unclassified".
tax <- data.frame(
  apply(
    tax,
    2,
    function(x) gsub("^$", "unclassified", x)
  ),
  row.names = rownames(tax),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Add genus information to assigned species names, as in the original analysis.
tax$species <- ifelse(
  grepl("^s:", tax$species),
  paste(
    gsub("g:", "s:", tax$genus),
    "_",
    gsub("s:", "", tax$species),
    sep = ""
  ),
  tax$species
)

# Remove chloroplast and mitochondrial sequences.
tax.select <- subset(
  tax,
  tax$order != "o:Chloroplast" &
    tax$family != "f:Mitochondria"
)

# Retain only bacterial and archaeal OTUs.
tax.select <- subset(
  tax.select,
  tax.select$domain == "d:Bacteria" |
    tax.select$domain == "d:Archaea"
)

###############################################################################
# 7. LOAD OTU COUNT TABLE
###############################################################################

som <- read.table(
  file = otu_table_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  comment.char = "",
  check.names = FALSE
)

som <- som[
  gtools::mixedorder(rownames(som)),
  ,
  drop = FALSE
]

# Retain OTUs present in both the count and taxonomy tables.
extract <- intersect(
  rownames(som),
  rownames(tax.select)
)

som.select <- som[
  extract,
  ,
  drop = FALSE
]

tax.select <- tax.select[
  extract,
  ,
  drop = FALSE
]

# Convert to samples in rows and OTUs in columns.
som.select <- t(som.select)

som.select <- som.select[
  gtools::mixedorder(rownames(som.select)),
  ,
  drop = FALSE
]

storage.mode(som.select) <- "numeric"

###############################################################################
# 8. MATCH SAMPLE DESIGN AND OTU TABLE
###############################################################################

samples_missing_from_design <- setdiff(
  rownames(som.select),
  design$Sample
)

samples_missing_from_otu_table <- setdiff(
  design$Sample,
  rownames(som.select)
)

if (length(samples_missing_from_design) > 0) {
  stop(
    paste(
      "Samples present in the OTU table but missing from design.txt:",
      paste(samples_missing_from_design, collapse = ", ")
    )
  )
}

if (length(samples_missing_from_otu_table) > 0) {
  warning(
    paste(
      "Samples present in design.txt but missing from the OTU table:",
      paste(samples_missing_from_otu_table, collapse = ", ")
    )
  )
}

# Retain and order the design rows corresponding to the OTU table.
design <- design[
  match(rownames(som.select), design$Sample),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(som.select),
    design$Sample
  )
)

stopifnot(
  identical(
    colnames(som.select),
    rownames(tax.select)
  )
)

###############################################################################
# 9. OTU FILTERING
###############################################################################

# Minimum percentage abundance across the complete dataset.
som.select <- som.select[
  ,
  100 / sum(som.select) * colSums(som.select) >= threshold_rare,
  drop = FALSE
]

# Minimum number of samples in which an OTU must occur.
som.select <- som.select[
  ,
  apply(
    som.select,
    2,
    function(x) sum(x > 0)
  ) >= threshold_sparse,
  drop = FALSE
]

# Ratio between the second-most abundant and most abundant sample.
som.select <- som.select[
  ,
  apply(
    som.select,
    2,
    function(x) {
      sorted_x <- sort(x, decreasing = TRUE)
      sorted_x[2] / sorted_x[1]
    }
  ) >= threshold_outlier,
  drop = FALSE
]

tax.select <- tax.select[
  colnames(som.select),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(som.select),
    design$Sample
  )
)

stopifnot(
  identical(
    colnames(som.select),
    rownames(tax.select)
  )
)

###############################################################################
# 10. SEQUENCING DEPTH
###############################################################################

sample_depths <- rowSums(som.select)

min_depth <- min(sample_depths)
min_sample <- names(which.min(sample_depths))

cat(
  "Minimum sequencing depth:",
  min_depth,
  "reads\n"
)

cat(
  "Sample with minimum sequencing depth:",
  min_sample,
  "\n"
)

sequencing_depth_table <- tibble::tibble(
  Sample = names(sample_depths),
  sequencing_depth = as.numeric(sample_depths)
) %>%
  dplyr::arrange(.data$sequencing_depth)

readr::write_tsv(
  sequencing_depth_table,
  file.path(
    output_dir,
    "sequencing_depth_per_sample.tsv"
  )
)

###############################################################################
# 11. RAREFACTION
###############################################################################

cat(
  "Running",
  iterations,
  "rarefaction iterations at",
  min_depth,
  "reads per sample.\n"
)

rarefaction_cores <- if (.Platform$OS.type == "windows") {
  1
} else {
  2
}

# No random seed is set here because the original analysis did not use one.
som.select.sub <- parallel::mclapply(
  as.list(seq_len(iterations)),
  function(x) {
    vegan::rrarefy(
      som.select,
      min_depth
    )
  },
  mc.cores = rarefaction_cores
)

###############################################################################
# 12. CALCULATE ALPHA DIVERSITY FOR EACH RAREFACTION
###############################################################################

use_cores <- if (.Platform$OS.type == "windows") {
  1
} else {
  max(
    1,
    parallel::detectCores() - 1
  )
}

calc_alpha_one <- function(x, sample_ids) {

  x <- as.matrix(x)
  storage.mode(x) <- "numeric"

  # vegan expects samples in rows.
  row_match <- if (!is.null(rownames(x))) {
    sum(rownames(x) %in% sample_ids)
  } else {
    0
  }

  column_match <- if (!is.null(colnames(x))) {
    sum(colnames(x) %in% sample_ids)
  } else {
    0
  }

  if (column_match > row_match) {
    x <- t(x)
  }

  sobs <- vegan::specnumber(x)

  shannon <- vegan::diversity(
    x,
    index = "shannon"
  )

  if (is.null(names(sobs)) || is.null(names(shannon))) {
    stop(
      paste(
        "Alpha-diversity calculations produced unnamed vectors.",
        "Check the row and column names of the OTU table."
      )
    )
  }

  tibble::tibble(
    Sample = names(sobs),
    sobs = as.numeric(sobs),
    shannon = as.numeric(shannon)
  )
}

if (use_cores == 1) {

  alpha_long_df <- dplyr::bind_rows(
    lapply(
      som.select.sub,
      calc_alpha_one,
      sample_ids = design$Sample
    ),
    .id = "raref_id"
  )

} else {

  alpha_parts <- parallel::mclapply(
    som.select.sub,
    calc_alpha_one,
    sample_ids = design$Sample,
    mc.cores = use_cores
  )

  alpha_long_df <- dplyr::bind_rows(
    alpha_parts,
    .id = "raref_id"
  )
}

stopifnot(
  is.data.frame(alpha_long_df)
)

stopifnot(
  all(
    c(
      "Sample",
      "sobs",
      "shannon"
    ) %in% names(alpha_long_df)
  )
)

stopifnot(
  nrow(alpha_long_df) > 0
)

###############################################################################
# 13. MEAN ALPHA DIVERSITY ACROSS THE 100 RAREFACTIONS
###############################################################################

alpha_df <- alpha_long_df %>%
  dplyr::summarise(
    sobs = mean(
      .data$sobs,
      na.rm = TRUE
    ),
    shannon = mean(
      .data$shannon,
      na.rm = TRUE
    ),
    .by = .data$Sample
  ) %>%
  dplyr::mutate(
    evenness = dplyr::if_else(
      .data$sobs > 1,
      .data$shannon / log(.data$sobs),
      NA_real_
    )
  )

stopifnot(
  "Sample" %in% names(alpha_df)
)

stopifnot(
  nrow(alpha_df) > 1
)

###############################################################################
# 14. JOIN ALPHA DIVERSITY WITH EXPERIMENTAL DESIGN
###############################################################################

df_alpha <- design %>%
  dplyr::left_join(
    alpha_df,
    by = "Sample"
  ) %>%
  droplevels()

cat(
  "Matched samples:",
  sum(!is.na(df_alpha$sobs)),
  "out of",
  nrow(df_alpha),
  "\n"
)

if (sum(!is.na(df_alpha$sobs)) < nrow(df_alpha)) {

  warning(
    paste(
      "Some samples in design.txt did not match the alpha-diversity",
      "results. Check sample naming."
    )
  )

  print(
    df_alpha %>%
      dplyr::filter(
        is.na(.data$sobs)
      ) %>%
      dplyr::select(
        .data$Sample,
        .data$date,
        .data$block,
        .data$treatment,
        .data$system
      )
  )
}

df_alpha <- df_alpha %>%
  dplyr::filter(
    !is.na(.data$sobs),
    !is.na(.data$shannon)
  ) %>%
  droplevels()

###############################################################################
# 15. BESTNORMALIZE TRANSFORMATIONS
###############################################################################

df2 <- df_alpha %>%
  dplyr::mutate(
    treatment = factor(.data$treatment),
    system = factor(.data$system)
  ) %>%
  dplyr::filter(
    !is.na(.data$sobs),
    !is.na(.data$shannon),
    !is.na(.data$evenness)
  ) %>%
  droplevels()

bn_sobs <- bestNormalize::bestNormalize(
  df2$sobs
)

bn_shannon <- bestNormalize::bestNormalize(
  df2$shannon
)

bn_evenness <- bestNormalize::bestNormalize(
  df2$evenness
)

df2 <- df2 %>%
  dplyr::mutate(
    sobs_t = predict(bn_sobs),
    shannon_t = predict(bn_shannon),
    evenness_t = predict(bn_evenness)
  )

###############################################################################
# 16. TWO-WAY ANOVA
###############################################################################

model_sobs <- aov(
  sobs_t ~ treatment * system,
  data = df2
)

model_shannon <- aov(
  shannon_t ~ treatment * system,
  data = df2
)

model_evenness <- aov(
  evenness_t ~ treatment * system,
  data = df2
)

cat("\nObserved richness ANOVA\n")
print(summary(model_sobs))

cat("\nShannon diversity ANOVA\n")
print(summary(model_shannon))

cat("\nPielou evenness ANOVA\n")
print(summary(model_evenness))

###############################################################################
# 17. SUMMARY TABLES
###############################################################################

alpha_summary_by_profile <- df2 %>%
  dplyr::group_by(
    .data$treatment
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_sobs = mean(
      .data$sobs,
      na.rm = TRUE
    ),
    median_sobs = median(
      .data$sobs,
      na.rm = TRUE
    ),
    mean_shannon = mean(
      .data$shannon,
      na.rm = TRUE
    ),
    median_shannon = median(
      .data$shannon,
      na.rm = TRUE
    ),
    mean_evenness = mean(
      .data$evenness,
      na.rm = TRUE
    ),
    median_evenness = median(
      .data$evenness,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

alpha_summary_by_peatland <- df2 %>%
  dplyr::group_by(
    .data$system,
    .data$treatment
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_sobs = mean(
      .data$sobs,
      na.rm = TRUE
    ),
    median_sobs = median(
      .data$sobs,
      na.rm = TRUE
    ),
    mean_shannon = mean(
      .data$shannon,
      na.rm = TRUE
    ),
    median_shannon = median(
      .data$shannon,
      na.rm = TRUE
    ),
    mean_evenness = mean(
      .data$evenness,
      na.rm = TRUE
    ),
    median_evenness = median(
      .data$evenness,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

###############################################################################
# 18. WRITE OUTPUTS
###############################################################################

readr::write_tsv(
  alpha_long_df,
  file.path(
    output_dir,
    "alpha_diversity_all_rarefactions.tsv"
  )
)

readr::write_tsv(
  alpha_df,
  file.path(
    output_dir,
    "alpha_diversity_mean_across_rarefactions.tsv"
  )
)

readr::write_tsv(
  df2,
  file.path(
    output_dir,
    "alpha_diversity_with_design_and_transformations.tsv"
  )
)

readr::write_tsv(
  alpha_summary_by_profile,
  file.path(
    output_dir,
    "alpha_diversity_summary_by_profile.tsv"
  )
)

readr::write_tsv(
  alpha_summary_by_peatland,
  file.path(
    output_dir,
    "alpha_diversity_summary_by_peatland.tsv"
  )
)

anova_output <- c(
  "OBSERVED RICHNESS",
  capture.output(summary(model_sobs)),
  "",
  "SHANNON DIVERSITY",
  capture.output(summary(model_shannon)),
  "",
  "PIELOU EVENNESS",
  capture.output(summary(model_evenness))
)

writeLines(
  anova_output,
  file.path(
    output_dir,
    "alpha_diversity_two_way_ANOVA.txt"
  )
)

normalization_output <- c(
  "OBSERVED RICHNESS",
  capture.output(print(bn_sobs)),
  "",
  "SHANNON DIVERSITY",
  capture.output(print(bn_shannon)),
  "",
  "PIELOU EVENNESS",
  capture.output(print(bn_evenness))
)

writeLines(
  normalization_output,
  file.path(
    output_dir,
    "alpha_diversity_bestNormalize_results.txt"
  )
)

# Save the exact rarefied tables and processed objects for the subsequent
# beta-diversity and figure-generation scripts.
saveRDS(
  list(
    design = design,
    taxonomy = tax.select,
    otu_table = som.select,
    rarefied_otu_tables = som.select.sub,
    alpha_diversity = alpha_df,
    alpha_diversity_with_design = df2,
    models = list(
      observed_richness = model_sobs,
      shannon = model_shannon,
      evenness = model_evenness
    )
  ),
  file.path(
    output_dir,
    "alpha_diversity_analysis_objects.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

cat(
  "\nAlpha-diversity analysis completed successfully.\n"
)

cat(
  "Results written to:",
  normalizePath(output_dir),
  "\n"
)
