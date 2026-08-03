#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Supplementary Figure 3, showing the summed relative abundance
# and transcriptional activity of EET-bearing Acidobacteriota and
# Verrucomicrobiota across peatland depths.
#
# INPUT
# EET_depth_plot_input.tsv
#
# Required columns:
#   metric
#   omics
#   site
#   depth_label
#   depth_mid_cm
#   phylum
#   total_signal
#
# OUTPUT
# Supplementary_Figure_3_plot_data.tsv
# Supplementary_Figure_3.svg
# Supplementary_Figure_3.pdf
#
# USAGE
# Rscript Supplementary_Figure_3.R \
#   EET_depth_plot_input.tsv \
#   output_directory

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    paste(
      "Usage: Rscript Supplementary_Figure_3.R",
      "<EET_depth_plot_input.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

INPUT_TSV <- args[[1]]
OUTPUT_DIR <- args[[2]]

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_PLOT_DATA <- file.path(
  OUTPUT_DIR,
  "Supplementary_Figure_3_plot_data.tsv"
)

OUT_SVG <- file.path(
  OUTPUT_DIR,
  "Supplementary_Figure_3.svg"
)

OUT_PDF <- file.path(
  OUTPUT_DIR,
  "Supplementary_Figure_3.pdf"
)


###############################################################################
# COLORS
###############################################################################

phylum_colors <- c(
  "Acidobacteriota"   = "#B2182B",
  "Verrucomicrobiota" = "#D8B365"
)

phylum_levels <- c(
  "Acidobacteriota",
  "Verrucomicrobiota"
)

###############################################################################
# PEATLAND ORDER AND LABELS
###############################################################################

site_short <- c(
  "Lungsmossen"  = "LM",
  "NorraRomyren" = "NR",
  "Havsjomossen" = "HV",
  "Bjorsmossen"  = "BM"
)

peatland_order <- c(
  "LM",
  "NR",
  "HV",
  "BM"
)

###############################################################################
# DEPTH ORDER
###############################################################################

# Used only to sort depths that truly exist within each peatland.
depth_order <- c(
  "10",
  "30",
  "60",
  "0-25",
  "100-125",
  "200-225",
  "265-280",
  "275-300",
  "325-350",
  "525-540",
  "525-550",
  "650-675",
  "675-700"
)

###############################################################################
# AXIS SETTINGS
###############################################################################

ABUNDANCE_TICK_STEP <- 5
ACTIVITY_TICK_STEP <- 0.05

next_axis_tick <- function(max_value, tick_step) {
  (floor(max_value / tick_step) + 1) * tick_step
}

###############################################################################
# READ DATA
###############################################################################

if (!file.exists(INPUT_TSV)) {
  stop("Input file not found: ", INPUT_TSV)
}

dat <- read_tsv(
  INPUT_TSV,
  show_col_types = FALSE
)

required_columns <- c(
  "metric",
  "omics",
  "site",
  "depth_label",
  "depth_mid_cm",
  "phylum",
  "total_signal"
)

missing_columns <- setdiff(required_columns, colnames(dat))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

dat <- dat %>%
  mutate(
    site = as.character(site),
    depth_label = as.character(depth_label),
    depth_mid_cm = as.numeric(depth_mid_cm),
    phylum = as.character(phylum),
    total_signal = as.numeric(total_signal)
  )


###############################################################################
# KEEP TARGET PHYLA AND ADD PEATLAND ABBREVIATIONS
###############################################################################

dat <- dat %>%
  filter(
    phylum %in% phylum_levels
  ) %>%
  mutate(
    peatland = unname(site_short[site])
  )


dat <- dat %>%
  mutate(
    peatland = factor(
      peatland,
      levels = peatland_order
    ),
    phylum = factor(
      phylum,
      levels = phylum_levels
    )
  )

###############################################################################
# SPLIT ABUNDANCE AND ACTIVITY
###############################################################################

abundance_df <- dat %>%
  filter(metric == "Relative abundance")

activity_df <- dat %>%
  filter(metric == "MAG transcriptional activity")



###############################################################################
# MASTER ROW STRUCTURE FROM ACTUAL metaG SAMPLES
###############################################################################

# Only peatland-depth combinations that truly exist in metaG are kept.
master_rows <- abundance_df %>%
  distinct(
    peatland,
    depth_label
  ) %>%
  mutate(
    depth_rank = match(depth_label, depth_order),
    row_id = paste(
      as.character(peatland),
      depth_label,
      sep = "__"
    )
  ) %>%
  arrange(
    peatland,
    depth_rank
  )

master_row_levels <- rev(master_rows$row_id)

master_rows <- master_rows %>%
  mutate(
    row_id = factor(
      row_id,
      levels = master_row_levels
    )
  )

###############################################################################
# COMPLETE EACH PANEL AGAINST MASTER ROWS
###############################################################################

prepare_panel_data <- function(source_df, target_phylum, panel_name) {
  
  source_sub <- source_df %>%
    filter(phylum == target_phylum) %>%
    select(
      peatland,
      depth_label,
      total_signal
    )
  
  duplicate_rows <- source_sub %>%
    count(
      peatland,
      depth_label
    ) %>%
    filter(n > 1)
  
  if (nrow(duplicate_rows) > 0) {
    stop(
      "Duplicated peatland-depth rows found for panel ",
      panel_name
    )
  }
  
  master_rows %>%
    left_join(
      source_sub,
      by = c("peatland", "depth_label")
    ) %>%
    mutate(
      phylum = factor(target_phylum, levels = phylum_levels),
      panel = panel_name,
      row_id = factor(
        as.character(row_id),
        levels = master_row_levels
      ),
      missing_sample = is.na(total_signal)
    )
}

abundance_acido <- prepare_panel_data(
  source_df = abundance_df,
  target_phylum = "Acidobacteriota",
  panel_name = "Acidobacteriota abundance"
)

abundance_verruco <- prepare_panel_data(
  source_df = abundance_df,
  target_phylum = "Verrucomicrobiota",
  panel_name = "Verrucomicrobiota abundance"
)

activity_acido <- prepare_panel_data(
  source_df = activity_df,
  target_phylum = "Acidobacteriota",
  panel_name = "Acidobacteriota activity"
)

activity_verruco <- prepare_panel_data(
  source_df = activity_df,
  target_phylum = "Verrucomicrobiota",
  panel_name = "Verrucomicrobiota activity"
)

###############################################################################
# COMMON X LIMITS
###############################################################################

abundance_max <- max(
  c(
    abundance_acido$total_signal,
    abundance_verruco$total_signal
  ),
  na.rm = TRUE
)

activity_max <- max(
  c(
    activity_acido$total_signal,
    activity_verruco$total_signal
  ),
  na.rm = TRUE
)


abundance_limit <- next_axis_tick(
  abundance_max,
  ABUNDANCE_TICK_STEP
)

activity_limit <- next_axis_tick(
  activity_max,
  ACTIVITY_TICK_STEP
)

###############################################################################
# SAVE PLOT DATA
###############################################################################

plot_input <- bind_rows(
  abundance_acido,
  abundance_verruco,
  activity_acido,
  activity_verruco
)

write_tsv(
  plot_input,
  OUT_PLOT_DATA
)

###############################################################################
# COMMON THEME
###############################################################################

common_theme <- theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 10.5,
      hjust = 0
    ),
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 8.5,
      color = "black",
      margin = margin(t = 6)
    ),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.10, "cm"),
    strip.placement = "outside",
    strip.background.y = element_rect(
      fill = "grey85",
      color = "black",
      linewidth = 0.45
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 9,
      color = "black",
      margin = margin(4, 5, 4, 5)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 9,
      color = "black"
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

###############################################################################
# PLOT FUNCTION
###############################################################################

make_panel <- function(
    df,
    panel_title,
    x_title,
    x_limit,
    tick_step,
    label_accuracy,
    show_y_labels = TRUE,
    show_strips = TRUE,
    show_missing_asterisks = FALSE
) {
  
  p <- ggplot(
    df,
    aes(
      x = total_signal,
      y = row_id,
      fill = phylum
    )
  ) +
    geom_vline(
      xintercept = 0,
      color = "black",
      linewidth = 0.45
    ) +
    geom_col(
      color = "black",
      linewidth = 0.4,
      width = 0.76,
      na.rm = TRUE,
      show.legend = TRUE
    ) +
    facet_grid(
      rows = vars(peatland),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      drop = TRUE
    ) +
    scale_y_discrete(
      drop = TRUE,
      labels = function(x) sub("^.*__", "", x)
    ) +
    scale_x_continuous(
      limits = c(0, x_limit),
      breaks = seq(0, x_limit, by = tick_step),
      labels = label_number(
        accuracy = label_accuracy,
        trim = TRUE
      ),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(
      values = phylum_colors,
      limits = phylum_levels,
      drop = FALSE
    ) +
    labs(
      title = panel_title,
      x = x_title,
      y = NULL
    ) +
    common_theme
  
  if (show_missing_asterisks) {
    p <- p +
      geom_text(
        data = df %>% filter(missing_sample),
        aes(
          x = x_limit * 0.05,
          y = row_id
        ),
        label = "*",
        inherit.aes = FALSE,
        size = 7,
        fontface = "bold",
        color = "black"
      )
  }
  
  if (!show_y_labels) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  if (!show_strips) {
    p <- p +
      theme(
        strip.background.y = element_blank(),
        strip.text.y.left = element_blank()
      )
  }
  
  p
}

###############################################################################
# BUILD PANELS IN THE REQUESTED ORDER
###############################################################################

p1 <- make_panel(
  df = abundance_acido,
  panel_title = "A  Acidobacteriota abundance (n = 96)",
  x_title = "Summed relative abundance (%)",
  x_limit = abundance_limit,
  tick_step = ABUNDANCE_TICK_STEP,
  label_accuracy = 1,
  show_y_labels = TRUE,
  show_strips = TRUE,
  show_missing_asterisks = FALSE
)

p2 <- make_panel(
  df = abundance_verruco,
  panel_title = "B  Verrucomicrobiota abundance (n = 24)",
  x_title = "Summed relative abundance (%)",
  x_limit = abundance_limit,
  tick_step = ABUNDANCE_TICK_STEP,
  label_accuracy = 1,
  show_y_labels = FALSE,
  show_strips = FALSE,
  show_missing_asterisks = FALSE
)

p3 <- make_panel(
  df = activity_acido,
  panel_title = "C  Acidobacteriota activity (n = 96)",
  x_title = "Summed transcriptional activity",
  x_limit = activity_limit,
  tick_step = ACTIVITY_TICK_STEP,
  label_accuracy = 0.01,
  show_y_labels = FALSE,
  show_strips = FALSE,
  show_missing_asterisks = TRUE
)

p4 <- make_panel(
  df = activity_verruco,
  panel_title = "D  Verrucomicrobiota activity (n = 24)",
  x_title = "Summed transcriptional activity",
  x_limit = activity_limit,
  tick_step = ACTIVITY_TICK_STEP,
  label_accuracy = 0.01,
  show_y_labels = FALSE,
  show_strips = FALSE,
  show_missing_asterisks = TRUE
)

###############################################################################
# COMBINE
###############################################################################

final_plot <- (
  p1 | p2 | p3 | p4
) +
  plot_layout(
    widths = c(1.18, 1, 1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

###############################################################################
# SAVE
###############################################################################

ggsave(
  filename = OUT_SVG,
  plot = final_plot,
  width = 18,
  height = 9.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

ggsave(
  filename = OUT_PDF,
  plot = final_plot,
  width = 9,
  height = 6.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

print(final_plot)
