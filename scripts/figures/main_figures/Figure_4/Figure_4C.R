#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Figure 4C, comparing metagenomic and
# metatranscriptomic coverage per cell for the EET1-EET4 gene subgroups
# across peatlands and depth layers.
#
# A common set of activity thresholds is calculated from the pooled positive
# MetaG and MetaT values. MetaT depths without a corresponding sample are
# displayed as not determined.
#
# INPUT
# 1. MetaT EET-gene table.
# 2. MetaG EET-gene table.
# 3. Output directory.
#
# OUTPUT
# Figure 4C as SVG and PDF and the plotted values as TSV.
#
# USAGE
# Rscript Figure_4C.R \
#   EET_genes_MetaT.tsv \
#   EET_genes_MetaG.tsv \
#   output_directory

suppressPackageStartupMessages({
  library(tidyverse)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript Figure_4C.R",
      "<MetaT_EET_table.tsv>",
      "<MetaG_EET_table.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

infile_mt <- args[[1]]
infile_mg <- args[[2]]
output_dir <- args[[3]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "Figure_4C_EET_subgroups_MetaG_MetaT.svg"
)

out_pdf <- file.path(
  output_dir,
  "Figure_4C_EET_subgroups_MetaG_MetaT.pdf"
)

out_tsv <- file.path(
  output_dir,
  "Figure_4C_EET_subgroups_MetaG_MetaT_values.tsv"
)

peatland_levels <- c("LM", "NR", "HV", "BM")
macro_levels <- c("EET1", "EET2", "EET3", "EET4")
dataset_levels <- c("MG", "MT")

depth_levels <- c(
  "0–10 cm", "10–30 cm", "30–60 cm", "0–25 cm",
  "100–125 cm", "200–225 cm", "265–280 cm",
  "275–300 cm", "325–350 cm",
  "525–540 cm", "525–550 cm", "650–675 cm", "675–700 cm"
)

y_levels <- c(
  "675–700 cm", "650–675 cm", "525–550 cm", "525–540 cm",
  "325–350 cm", "275–300 cm",
  "265–280 cm", "200–225 cm", "100–125 cm",
  "0–25 cm",
  "30–60 cm", "10–30 cm", "0–10 cm"
)

cream <- "#FCF8F2"
nd_fill <- "#FCF8F2"

fill_map <- c(
  "ND"          = nd_fill,
  "Zero"        = cream,
  "EET1_Low"    = "#F2E2DD",
  "EET1_Medium" = "#D9A79B",
  "EET1_High"   = "#9B4C45",
  "EET2_Low"    = "#E4F0E3",
  "EET2_Medium" = "#A8CFA3",
  "EET2_High"   = "#4F8A44",
  "EET3_Low"    = "#E3EDF8",
  "EET3_Medium" = "#A8C4E6",
  "EET3_High"   = "#4B7DB8",
  "EET4_Low"    = "#EEE5F7",
  "EET4_Medium" = "#C5A8DD",
  "EET4_High"   = "#7B57A5"
)

clean_depth <- function(x) {
  x %>%
    stringr::str_replace_all("–", "-") %>%
    stringr::str_replace_all("—", "-") %>%
    stringr::str_replace_all(" ", "") %>%
    stringr::str_to_lower()
}

depth_map <- tibble::tribble(
  ~Depth_raw,   ~Depth_layer,
  "10cm",       "0–10 cm",
  "0-10cm",     "0–10 cm",
  "30cm",       "10–30 cm",
  "60cm",       "30–60 cm",
  "0-25cm",     "0–25 cm",
  "100-125cm",  "100–125 cm",
  "200-225cm",  "200–225 cm",
  "265-280cm",  "265–280 cm",
  "275-300cm",  "275–300 cm",
  "325-350cm",  "325–350 cm",
  "525-540cm",  "525–540 cm",
  "525-550cm",  "525–550 cm",
  "650-675cm",  "650–675 cm",
  "675-700cm",  "675–700 cm"
) %>%
  dplyr::mutate(
    Depth_raw = clean_depth(.data$Depth_raw)
  )

extract_peatland <- function(sample_name) {
  dplyr::case_when(
    stringr::str_detect(
      sample_name,
      stringr::fixed("Bjorsmossen", ignore_case = TRUE)
    ) ~ "BM",
    stringr::str_detect(
      sample_name,
      stringr::fixed("NorraRomyren", ignore_case = TRUE)
    ) ~ "NR",
    stringr::str_detect(
      sample_name,
      stringr::fixed("Havsjomossen", ignore_case = TRUE)
    ) ~ "HV",
    stringr::str_detect(
      sample_name,
      stringr::fixed("Lungsmossen", ignore_case = TRUE)
    ) ~ "LM",
    TRUE ~ NA_character_
  )
}


get_subgroup_order <- function(df_subgroups) {
  df_subgroups %>%
    distinct(Macro_EET, Subgroup) %>%
    mutate(
      Subgroup_short = str_remove(Subgroup, "^EET[1-4]_"),
      subgroup_rank = case_when(
        Macro_EET == "EET1" & Subgroup_short == "porin"  ~ 1,
        Macro_EET == "EET1" & Subgroup_short == "NHL"    ~ 2,
        Macro_EET == "EET1" & Subgroup_short == "cytc"   ~ 3,
        
        Macro_EET == "EET2" & Subgroup_short == "mtrC"   ~ 1,
        Macro_EET == "EET2" & Subgroup_short == "dmsE"   ~ 2,
        Macro_EET == "EET2" & Subgroup_short == "porin"  ~ 3,
        Macro_EET == "EET2" & Subgroup_short == "cytc"   ~ 4,
        
        Macro_EET == "EET3" & Subgroup_short == "9-heme" ~ 1,
        Macro_EET == "EET3" & Subgroup_short == "8-heme" ~ 2,
        Macro_EET == "EET3" & Subgroup_short == "porin"  ~ 3,
        
        Macro_EET == "EET4" & Subgroup_short == "9-heme" ~ 1,
        Macro_EET == "EET4" & Subgroup_short == "8-heme" ~ 2,
        Macro_EET == "EET4" & Subgroup_short == "porin"  ~ 3,
        
        TRUE ~ 999
      )
    ) %>%
    arrange(Macro_EET, subgroup_rank, Subgroup)
}


df_mt <- readr::read_tsv(
  infile_mt,
  show_col_types = FALSE
)

mt_cols <- grep(
  "^TranscriptM_",
  names(df_mt),
  value = TRUE
)

mt_long <- df_mt %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(mt_cols),
    names_to = "Sample_col",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Value = tidyr::replace_na(.data$Value, 0),
    Sample = .data$Sample_col %>%
      stringr::str_remove("^TranscriptM_") %>%
      stringr::str_remove("_iMT_all_genes_with_MAG_info_tpm$"),
    Macro_EET = stringr::str_extract(
      .data$EET_type,
      "^EET[1-4]"
    ),
    Subgroup = .data$EET_type,
    Peatland = extract_peatland(.data$Sample),
    Depth_raw = clean_depth(
      stringr::str_extract(
        .data$Sample,
        "\\d+-\\d+cm|\\d+cm"
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$Macro_EET),
    !is.na(.data$Peatland)
  ) %>%
  dplyr::left_join(
    depth_map,
    by = "Depth_raw"
  ) %>%
  dplyr::mutate(
    Peatland = factor(
      .data$Peatland,
      levels = peatland_levels
    ),
    Depth_layer = factor(
      .data$Depth_layer,
      levels = depth_levels
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$Depth_layer)
  )

mt_agg <- mt_long %>%
  dplyr::group_by(
    .data$Peatland,
    .data$Depth_layer,
    .data$Macro_EET,
    .data$Subgroup
  ) %>%
  dplyr::summarise(
    Value_mean = mean(
      .data$Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Dataset = "MT"
  )

mt_pd_present <- mt_agg %>%
  dplyr::distinct(
    .data$Peatland,
    .data$Depth_layer
  )



df_mg <- readr::read_tsv(
  infile_mg,
  show_col_types = FALSE
)

mg_meta_cols <- c(
  "TXT_protein",
  "EET_type"
)

mg_sample_cols <- setdiff(
  names(df_mg),
  mg_meta_cols
)

mg_long <- df_mg %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(mg_sample_cols),
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Value = tidyr::replace_na(.data$Value, 0),
    Macro_EET = stringr::str_extract(
      .data$EET_type,
      "^EET[1-4]"
    ),
    Subgroup = .data$EET_type,
    Peatland = extract_peatland(.data$Sample),
    Depth_raw = clean_depth(
      stringr::str_extract(
        .data$Sample,
        "\\d+-\\d+cm|\\d+cm"
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$Macro_EET),
    !is.na(.data$Peatland)
  ) %>%
  dplyr::left_join(
    depth_map,
    by = "Depth_raw"
  ) %>%
  dplyr::mutate(
    Peatland = factor(
      .data$Peatland,
      levels = peatland_levels
    ),
    Depth_layer = factor(
      .data$Depth_layer,
      levels = depth_levels
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$Depth_layer)
  )

mg_agg <- mg_long %>%
  dplyr::group_by(
    .data$Peatland,
    .data$Depth_layer,
    .data$Macro_EET,
    .data$Subgroup
  ) %>%
  dplyr::summarise(
    Value_mean = mean(
      .data$Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Dataset = "MG"
  )

mg_pd_present <- mg_agg %>%
  dplyr::distinct(
    .data$Peatland,
    .data$Depth_layer
  )


############################################################
### 5) COMMON GRID (use MetaG depths as reference)
############################################################

subgroup_order_df <- bind_rows(
  mt_agg %>% select(Macro_EET, Subgroup),
  mg_agg %>% select(Macro_EET, Subgroup)
) %>%
  get_subgroup_order()

subgroup_grid <- subgroup_order_df %>%
  select(Macro_EET, Subgroup)

mg_plot <- mg_pd_present %>%
  crossing(subgroup_grid) %>%
  left_join(
    mg_agg,
    by = c("Peatland", "Depth_layer", "Macro_EET", "Subgroup")
  ) %>%
  mutate(
    Value_mean = replace_na(Value_mean, 0),
    Dataset = "MG",
    Status = "Measured"
  ) %>%
  filter(str_detect(Subgroup, paste0("^", Macro_EET, "_")))

mt_plot <- mg_pd_present %>%
  crossing(subgroup_grid) %>%
  left_join(
    mt_agg,
    by = c("Peatland", "Depth_layer", "Macro_EET", "Subgroup")
  ) %>%
  left_join(
    mt_pd_present %>% mutate(MT_depth_present = TRUE),
    by = c("Peatland", "Depth_layer")
  ) %>%
  mutate(
    Dataset = "MT",
    MT_depth_present = replace_na(MT_depth_present, FALSE),
    Status = case_when(
      !MT_depth_present ~ "ND",
      TRUE ~ "Measured"
    ),
    Value_mean = case_when(
      Status == "ND" ~ NA_real_,
      TRUE ~ replace_na(Value_mean, 0)
    )
  ) %>%
  filter(str_detect(Subgroup, paste0("^", Macro_EET, "_")))

plot_df <- bind_rows(mg_plot, mt_plot) %>%
  mutate(
    Peatland = factor(Peatland, levels = peatland_levels),
    Depth_layer = factor(Depth_layer, levels = depth_levels),
    Depth_layer_plot = factor(as.character(Depth_layer), levels = y_levels),
    Depth_label = factor(
      str_remove(as.character(Depth_layer_plot), " cm$"),
      levels = str_remove(y_levels, " cm$")
    ),
    Macro_EET = factor(Macro_EET, levels = macro_levels),
    Dataset = factor(Dataset, levels = dataset_levels),
    Panel = paste0(Macro_EET, "\n", Dataset),
    Panel = factor(
      Panel,
      levels = c(
        "EET1\nMG", "EET1\nMT",
        "EET2\nMG", "EET2\nMT",
        "EET3\nMG", "EET3\nMT",
        "EET4\nMG", "EET4\nMT"
      )
    ),
    Subgroup = factor(Subgroup, levels = subgroup_order_df$Subgroup),
    Subgroup_short = str_remove(as.character(Subgroup), "^EET[1-4]_")
  )

############################################################
### 6) COMMON BINNING ACROSS MG + MT
############################################################

positive_vals <- plot_df %>%
  filter(Status == "Measured", !is.na(Value_mean), Value_mean > 0) %>%
  pull(Value_mean)

thr_low <- unname(quantile(positive_vals, probs = 0.25, na.rm = TRUE))
thr_med <- unname(quantile(positive_vals, probs = 0.75, na.rm = TRUE))

bin_fun <- function(x, low, med, status) {
  case_when(
    status == "ND" ~ "ND",
    is.na(x) ~ NA_character_,
    x == 0 ~ "Zero",
    x > 0 & x <= low ~ "Low",
    x > low & x <= med ~ "Medium",
    x > med ~ "High"
  )
}

plot_df <- plot_df %>%
  mutate(
    Bin = bin_fun(Value_mean, thr_low, thr_med, Status),
    Bin = factor(Bin, levels = c("ND", "Zero", "Low", "Medium", "High")),
    Fill_group = case_when(
      Bin == "ND"     ~ "ND",
      Bin == "Zero"   ~ "Zero",
      Bin == "Low"    ~ paste0(Macro_EET, "_Low"),
      Bin == "Medium" ~ paste0(Macro_EET, "_Medium"),
      Bin == "High"   ~ paste0(Macro_EET, "_High"),
      TRUE ~ NA_character_
    )
  )

############################################################
### 7) PLOT
############################################################

p <- ggplot(
  plot_df,
  aes(x = Subgroup, y = Depth_label, fill = Fill_group)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_hline(yintercept = 4.5, colour = "black", linewidth = 0.35) +
  geom_text(
    data = plot_df %>% filter(Bin == "ND"),
    aes(label = "*"),
    size = 3.8,
    color = "black"
  ) +
  facet_grid(
    rows = vars(Peatland),
    cols = vars(Panel),
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  scale_x_discrete(labels = function(x) str_remove(x, "^EET[1-4]_")) +
  scale_fill_manual(
    values = fill_map,
    drop = FALSE,
    na.value = "white",
    name = "Coverage per cell",
    breaks = c("ND", "Zero", "EET1_Low", "EET1_Medium", "EET1_High"),
    labels = c("n.d.", "Zero", "Low", "Medium", "High"),
    guide = guide_legend(
      override.aes = list(fill = c("white", cream, "grey75", "grey45", "black"))
    )
  ) +
  labs(
    x = "EET subgroup",
    y = "Depth",
    title = "EET subgroup coverage per cell across peatlands and depths",
    subtitle = paste0(
      "* n.d. = no MetaT sample at the corresponding MetaG depth | ",
      "Common thresholds across MG + MT: Low ≤ ", signif(thr_low, 3),
      "; Medium ≤ ", signif(thr_med, 3),
      "; High > ", signif(thr_med, 3)
    )
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    
    strip.background = element_rect(fill = "grey95", colour = "grey40"),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.text.y.left = element_text(size = 13, angle = 0, face = "bold"),
    strip.placement = "outside",
    
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 9),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    
    panel.spacing.x = unit(0.35, "lines"),
    panel.spacing.y = unit(0.45, "lines"),
    legend.position = "right"
  )

print(p)


ggplot2::ggsave(
  filename = out_svg,
  plot = p,
  device = svglite::svglite,
  width = 10,
  height = 10,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

if (capabilities("cairo")) {
  ggplot2::ggsave(
    filename = out_pdf,
    plot = p,
    device = grDevices::cairo_pdf,
    width = 10,
    height = 10,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}

plot_df %>%
  dplyr::select(
    .data$Peatland,
    .data$Depth_layer,
    .data$Macro_EET,
    .data$Dataset,
    .data$Subgroup,
    .data$Value_mean,
    .data$Bin,
    .data$Status
  ) %>%
  dplyr::arrange(
    .data$Peatland,
    .data$Depth_layer,
    .data$Macro_EET,
    .data$Dataset,
    .data$Subgroup
  ) %>%
  readr::write_tsv(
    out_tsv
  )
