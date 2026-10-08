#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Extended Data Figure 2B: an NMDS ordination of metagenomic
# read datasets (Värmland, SPRUCE, Stordalen Mire) based on Mash distances.
#
# Duplicate pairwise comparisons are collapsed by retaining the minimum
# distance. Missing pairwise distances are assigned a distance of 1
#
# INPUT
# 1. Metadata CSV containing Sample and Peatland_type columns.
# 2. Mash distance table produced by `mash dist`.
# 3. Output directory.
#
# OUTPUT
# Extended Data Figure 2B as SVG and PDF and the fitted NMDS object as RDS.
#
# USAGE
# Rscript ED_Figure_2B_MASH.R \
#   metadata_MASH_reads.csv \
#   mash_distances.tab \
#   output_directory

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript ED_Figure_2B_MASH.R",
      "<metadata.csv>",
      "<mash_distances.tab>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

metadata_file <- args[[1]]
mash_file <- args[[2]]
output_dir <- args[[3]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "ED_Figure_2B_NMDS_MASH_reads.svg"
)

out_pdf <- file.path(
  output_dir,
  "ED_Figure_2B_NMDS_MASH_reads.pdf"
)

out_rds <- file.path(
  output_dir,
  "ED_Figure_2B_NMDS_MASH_object.rds"
)

clean_name <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.fq\\.gz$", ".fastq", x, ignore.case = TRUE)
  x <- sub("\\.fastq\\.gz$", ".fastq", x, ignore.case = TRUE)
  x <- sub("\\.(fa|fna|fasta)(\\.gz)?$", "", x, ignore.case = TRUE)
  x
}

metadata <- read.csv(
  metadata_file,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

metadata$Sample <- clean_name(
  metadata$Sample
)

distances <- data.table::fread(
  mash_file,
  header = FALSE,
  sep = "\t"
)

colnames(distances) <- c(
  "genome1",
  "genome2",
  "distance",
  "pval",
  "shared"
)

distances$genome1 <- clean_name(
  distances$genome1
)

distances$genome2 <- clean_name(
  distances$genome2
)

all_samples <- unique(
  metadata$Sample
)

distances <- distances[
  genome1 %in% all_samples &
    genome2 %in% all_samples,
]

distances <- distances[
  genome1 != genome2,
]

distances_unique <- distances %>%
  dplyr::mutate(
    pair_1 = pmin(.data$genome1, .data$genome2),
    pair_2 = pmax(.data$genome1, .data$genome2)
  ) %>%
  dplyr::group_by(
    .data$pair_1,
    .data$pair_2
  ) %>%
  dplyr::slice_min(
    order_by = .data$distance,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    -.data$pair_1,
    -.data$pair_2
  )

all_genomes <- sort(
  unique(
    c(
      distances_unique$genome1,
      distances_unique$genome2
    )
  )
)

dist_mat_full <- matrix(
  NA_real_,
  nrow = length(all_genomes),
  ncol = length(all_genomes),
  dimnames = list(
    all_genomes,
    all_genomes
  )
)

for (i in seq_len(nrow(distances_unique))) {
  genome_1 <- distances_unique$genome1[i]
  genome_2 <- distances_unique$genome2[i]
  mash_distance <- distances_unique$distance[i]

  dist_mat_full[
    genome_1,
    genome_2
  ] <- mash_distance

  dist_mat_full[
    genome_2,
    genome_1
  ] <- mash_distance
}

diag(dist_mat_full) <- 0

dist_mat_full[
  is.na(dist_mat_full)
] <- 1

dist_obj <- stats::as.dist(
  dist_mat_full
)

set.seed(123)

nmds <- vegan::metaMDS(
  dist_obj,
  k = 2,
  trymax = 3,
  autotransform = FALSE
)

saveRDS(
  nmds,
  out_rds
)

points <- as.data.frame(
  nmds$points
)

points$Sample <- rownames(
  points
)

plot_data <- merge(
  points,
  metadata,
  by = "Sample"
)

plot_data$Peatland_type <- factor(
  plot_data$Peatland_type,
  levels = c(
    "Varmland bog",
    "Stordalen bog",
    "Stordalen fen",
    "Stordalen palsa",
    "Spruce bog"
  )
)

p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$MDS1,
    y = .data$MDS2,
    fill = .data$Peatland_type
  )
) +
  ggplot2::geom_point(
    size = 3.6,
    alpha = 0.9,
    shape = 21,
    color = "black",
    stroke = 0.3
  ) +
  ggplot2::labs(
    title = "NMDS on metagenomic reads (MASH distances)",
    subtitle = paste(
      "Stress =",
      round(
        nmds$stress,
        3
      )
    ),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Varmland bog" = "#D73027",
      "Spruce bog" = "#56B4E9",
      "Stordalen bog" = "#00441B",
      "Stordalen fen" = "#1B7837",
      "Stordalen palsa" = "#A1D99B"
    )
  ) +
  ggplot2::coord_fixed(
    ratio = 1
  ) +
  ggplot2::theme_classic(
    base_size = 15
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 17,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = ggplot2::element_text(
      size = 13,
      hjust = 0.5,
      color = "gray30"
    ),
    axis.title = ggplot2::element_text(
      size = 14
    ),
    axis.text = ggplot2::element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "right",
    legend.title = ggplot2::element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(
      size = 12
    ),
    legend.key.size = grid::unit(
      0.9,
      "lines"
    ),
    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.4
    ),
    axis.ticks = ggplot2::element_line(
      color = "black",
      linewidth = 0.4
    ),
    plot.margin = ggplot2::margin(
      10,
      10,
      10,
      10
    )
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      title = "Peatland type",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(
        size = 4
      )
    )
  )

print(p)

ggplot2::ggsave(
  filename = out_svg,
  plot = p,
  width = 6,
  height = 5,
  units = "in",
  dpi = 300,
  device = svglite::svglite,
  bg = "white"
)

if (capabilities("cairo")) {
  ggplot2::ggsave(
    filename = out_pdf,
    plot = p,
    width = 6,
    height = 5,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}
