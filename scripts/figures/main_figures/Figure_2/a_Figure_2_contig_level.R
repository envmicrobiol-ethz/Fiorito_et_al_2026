#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Figure 2 contig-level metatranscriptomic functional heatmap
# based on MT_coverage_per_cell values from KEGG and METABOLIC annotations.
#
# Values are summed across samples and markers within each peatland, depth and
# metabolic function, and divided by the number of markers assigned to that
# function. Activity is represented using four fixed intensity classes.
#
# INPUT
# 1. Gene-catalogue KEGG expression-profile TSV.
# 2. Gene-catalogue METABOLIC expression-profile TSV.
# 3. Output directory.
#
# OUTPUT
# Figure 2A as SVG and PDF, an SVG without legends, and the plotted values TSV.
#
# USAGE
# Rscript Figure_2A.R \
#   genes_reps.KEGG.expression_profile.tsv \
#   genes_reps.METABOLIC.expression_profile.tsv \
#   figure_output_directory

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript Figure_2A.R",
      "<KEGG_expression.tsv>",
      "<METABOLIC_expression.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

kegg_file <- args[[1]]
metab_file <- args[[2]]
output_dir <- args[[3]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "Figure_2A_metatranscriptomic_functional_heatmap.svg"
)

out_svg_no_legend <- file.path(
  output_dir,
  "Figure_2A_metatranscriptomic_functional_heatmap_no_legend.svg"
)

out_pdf <- file.path(
  output_dir,
  "Figure_2A_metatranscriptomic_functional_heatmap.pdf"
)

out_values <- file.path(
  output_dir,
  "Figure_2A_metatranscriptomic_functional_heatmap_values.tsv"
)


############################################################
### 0b) SETTINGS
############################################################

peatland_levels <- c("LM","NR","HV","BM")

depth_levels <- c(
  "0–10 cm","10–30 cm","30–60 cm",
  "200–225 cm","265–280 cm","325–350 cm",
  "525–540 cm","525–550 cm","650–675 cm"
)

shallow_levels <- c("0–10 cm","10–30 cm","30–60 cm")

# order within shallow (0–10 on top)
shallow_y_levels <- c("30–60 cm","10–30 cm","0–10 cm")

# order within deep (200–225 on top of deep block; deepest at bottom)
deep_y_levels <- c("650–675 cm","525–550 cm","525–540 cm","325–350 cm","265–280 cm","200–225 cm")

# global y levels (deep block below + shallow block above)
y_levels <- c(deep_y_levels, shallow_y_levels)

# gap between shallow and deep + between peatlands (pt)
gap_pt <- 1.2

############################################################
### 1) Load MetaT — explicitly use MT_coverage_per_cell
############################################################

kegg_raw <- read_tsv(
  kegg_file,
  show_col_types = FALSE
)

metab_raw <- read_tsv(
  metab_file,
  show_col_types = FALSE
)


df_kegg <- kegg_raw %>%
  dplyr::transmute(
    Marker_ID = as.character(.data$KEGG_ko),
    Sample = as.character(.data$MetaT_sample),
    Coverage_per_cell = suppressWarnings(
      as.numeric(.data$MT_coverage_per_cell)
    ),
    Source = "KEGG"
  )

df_metab <- metab_raw %>%
  dplyr::transmute(
    Marker_ID = as.character(.data$METABOLIC_hmm),
    Sample = as.character(.data$MetaT_sample),
    Coverage_per_cell = suppressWarnings(
      as.numeric(.data$MT_coverage_per_cell)
    ),
    Source = "METABOLIC"
  )

df_MT <- dplyr::bind_rows(
  df_kegg,
  df_metab
) %>%
  dplyr::filter(
    !is.na(.data$Sample),
    !is.na(.data$Marker_ID),
    !is.na(.data$Coverage_per_cell)
  )


############################################################
### 1b) Extract peatland from Sample
############################################################

df_MT <- df_MT %>%
  dplyr::mutate(
    Peatland = dplyr::case_when(
      str_detect(Sample, regex("bj(ö|o)rsmossen|bjorsmossen", ignore_case = TRUE)) ~ "BM",
      str_detect(Sample, regex("norra[_-]?romyren|norraromyren|romyren", ignore_case = TRUE)) ~ "NR",
      str_detect(Sample, regex("havsj(ö|o)mossen|havsjomossen", ignore_case = TRUE)) ~ "HV",
      str_detect(Sample, regex("lungsmossen|lungs?mossen", ignore_case = TRUE)) ~ "LM",
      TRUE ~ NA_character_
    ),
    Peatland = factor(Peatland, levels = peatland_levels)
  ) %>%
  dplyr::filter(!is.na(Peatland))

############################################################
### 2) Depth cleaning + mapping
############################################################

clean_depth <- function(x){
  x %>%
    str_replace_all("–", "-") %>%
    str_replace_all("—", "-") %>%
    str_replace_all(" ", "") %>%
    str_to_lower()
}

df_MT <- df_MT %>%
  dplyr::mutate(Depth_raw = clean_depth(str_extract(Sample, "\\d+-\\d+cm|\\d+cm")))

depth_map <- tribble(
  ~Depth_raw,      ~Depth_layer,
  "10cm",          "0–10 cm",
  "0-10cm",        "0–10 cm",
  "30cm",          "10–30 cm",
  "60cm",          "30–60 cm",
  "200-225cm",     "200–225 cm",
  "265-280cm",     "265–280 cm",
  "325-350cm",     "325–350 cm",

  "525-540cm",     "525–540 cm",
  "525-550cm",     "525–550 cm",
  "650-675cm",     "650–675 cm"
) %>% dplyr::mutate(Depth_raw = clean_depth(Depth_raw))

df_MT <- df_MT %>%
  dplyr::left_join(depth_map, by = "Depth_raw") %>%
  dplyr::mutate(Depth_layer = factor(Depth_layer, levels = depth_levels)) %>%
  dplyr::filter(!is.na(Depth_layer))

############################################################
### 3) Marker table: Macro + Category + KO/HMM (ORDER MATTERS)
############################################################

# Macro-categories are used only for ordering and visual separation.
# They no longer receive different colours.

macro_levels <- c(
  "Complex Carbon metabolism",
  "Fermentation",
  "C1 metabolism",
  "Carbon fixation",
  "Methane metabolism",
  "Hydrogenases",
  "Oxidative phosphorylation",
  "Oxygen metabolism",
  "Iron cycling",
  "Manganese cycling",
  "Nitrogen cycling",
  "Sulfur cycling"
)

marker_tbl <- tribble(
  ~Macro, ~Function, ~Marker_ID_raw,
  
  "Complex Carbon metabolism","Amino acid utilization","K00823, K07250, K13524, K14268, K03918",
  "Complex Carbon metabolism","Amino acid utilization","K05825",
  "Complex Carbon metabolism","Amino acid utilization","K00831",
  "Complex Carbon metabolism","Amino acid utilization","K00819, K00821, K05830, K00840",
  "Complex Carbon metabolism","Amino acid utilization","K00826, K02619, K03342",
  "Complex Carbon metabolism","Amino acid utilization","K00812, K00813, K11358, K00832",
  "Complex Carbon metabolism","Amino acid utilization","K00817",
  "Complex Carbon metabolism","Amino acid utilization","K00830",
  
  "Complex Carbon metabolism","Aromatics degradation","K03381",
  "Complex Carbon metabolism","Aromatics degradation","K03186",
  "Complex Carbon metabolism","Aromatics degradation","K01612",
  "Complex Carbon metabolism","Aromatics degradation","K04112",
  "Complex Carbon metabolism","Aromatics degradation","K04113",
  "Complex Carbon metabolism","Aromatics degradation","K04114",
  "Complex Carbon metabolism","Aromatics degradation","K04115",
  
  "Complex Carbon metabolism","Complex carbon degradation","K19668",
  "Complex Carbon metabolism","Complex carbon degradation","K01179, K20542",
  "Complex Carbon metabolism","Complex carbon degradation","K01188, K05349, K05350",
  "Complex Carbon metabolism","Complex carbon degradation","K01209, K15921",
  "Complex Carbon metabolism","Complex carbon degradation","K01195",
  "Complex Carbon metabolism","Complex carbon degradation","K05989",
  "Complex Carbon metabolism","Complex carbon degradation","K01218, K19355",
  "Complex Carbon metabolism","Complex carbon degradation","K01811",
  "Complex Carbon metabolism","Complex carbon degradation","K01198",
  "Complex Carbon metabolism","Complex carbon degradation","K01192",
  "Complex Carbon metabolism","Complex carbon degradation","K01190, K12111, K12308",
  "Complex Carbon metabolism","Complex carbon degradation","K01176, K05343, K07405",
  "Complex Carbon metabolism","Complex carbon degradation","K01178",
  "Complex Carbon metabolism","Complex carbon degradation","K01200",
  "Complex Carbon metabolism","Complex carbon degradation","K01214",
  "Complex Carbon metabolism","Complex carbon degradation","K01183, K13381",
  "Complex Carbon metabolism","Complex carbon degradation","K01207, K12373, K14459",
  
  "Fermentation","Alcohol metabolism","K00129, K00138",
  "Fermentation","Alcohol metabolism","K13954, K00001, K04072, K00114, K00002, K04022, K22473",
  "Fermentation","Alcohol metabolism","K00001",
  "Fermentation","Lactate metabolism","K00016",
  "Fermentation","Acetate metabolism","K01905",
  "Fermentation","Acetate metabolism","K00925",
  "Fermentation","Acetate metabolism","K00625",
  "Fermentation","Acetate metabolism","K01895",
  
  "C1 metabolism","Methanol oxidation","K14028",
  "C1 metabolism","Methanol oxidation","K00093",
  "C1 metabolism","Methylamine oxidation","K15228",
  "C1 metabolism","Methylamine oxidation","K15229",
  "C1 metabolism","Formaldehyde oxidation","K00148",
  "C1 metabolism","Formaldehyde oxidation","K01070",
  "C1 metabolism","Formaldehyde oxidation","K00121",
  "C1 metabolism","Formaldehyde oxidation","K00153",
  "C1 metabolism","Formaldehyde oxidation","K10713",
  "C1 metabolism","Formate metabolism","K00123",
  "C1 metabolism","Formate metabolism","K22515",
  "C1 metabolism","Formate metabolism","K00124",
  "C1 metabolism","Formate metabolism","K22516",
  "C1 metabolism","Formate metabolism","K00125",
  "C1 metabolism","CO (aerobic)","K03518",
  "C1 metabolism","CO (aerobic)","K03519",
  "C1 metabolism","CO (aerobic)","K03520",
  
  "Carbon fixation","CBB cycle - Rubisco","K01601",
  "Carbon fixation","CBB cycle - Rubisco","K01601",
  "Carbon fixation","3 Hydroxypropionate cycle","K14469",
  "Carbon fixation","3 Hydroxypropionate cycle","K14468",
  "Carbon fixation","Wood Ljungdahl pathway","K00194",
  "Carbon fixation","Wood Ljungdahl pathway","K00197",
  "Carbon fixation","Wood Ljungdahl pathway","K00198",
  "Carbon fixation","Reverse TCA cycle","K15230",
  "Carbon fixation","Reverse TCA cycle","K15231",
  
  "Methane metabolism","Methanotrophy","K10944",
  "Methane metabolism","Methanotrophy","K10945",
  "Methane metabolism","Methanotrophy","K10946",
  "Methane metabolism","Methanotrophy","K16160",
  "Methane metabolism","Methanotrophy","K16162",
  "Methane metabolism","Methanogenesis","K00399",
  "Methane metabolism","Methanogenesis","K00401",
  "Methane metabolism","Methanogenesis","K03421",
  
  "Hydrogenases","FeFe hydrogenase","fefe-group-a13.hmm",
  "Hydrogenases","FeFe hydrogenase","fefe-group-a2.hmm",
  "Hydrogenases","FeFe hydrogenase","fefe-group-a4.hmm",
  "Hydrogenases","FeFe hydrogenase","fefe-group-b.hmm",
  "Hydrogenases","FeFe hydrogenase","fefe-group-c1.hmm",
  "Hydrogenases","FeFe hydrogenase","fefe-group-c2.hmm",
  "Hydrogenases","FeFe hydrogenase","fefe-group-c3.hmm",
  "Hydrogenases","Ni-Fe Hydrogenase","nife-group-1.hmm",
  "Hydrogenases","Ni-Fe Hydrogenase","nife-group-2ade.hmm",
  "Hydrogenases","Ni-Fe Hydrogenase","nife-group-2bc.hmm",
  "Hydrogenases","Ni-Fe Hydrogenase","nife-group-3abd.hmm",
  "Hydrogenases","Ni-Fe Hydrogenase","nife-group-3c.hmm",
  "Hydrogenases","Ni-Fe Hydrogenase","nife-group-4a-g.hmm",
  "Hydrogenases","Ni-Fe Hydrogenase","nife-group-4hi.hmm",
  
  "Oxidative phosphorylation","Complex I","K00330",
  "Oxidative phosphorylation","Complex I","K00331",
  "Oxidative phosphorylation","Complex I","K00332",
  "Oxidative phosphorylation","Complex I","K05572",
  "Oxidative phosphorylation","Complex I","K05573",
  "Oxidative phosphorylation","Complex I","K05574",
  "Oxidative phosphorylation","Complex II","K00241",
  "Oxidative phosphorylation","Complex II","K00242",
  "Oxidative phosphorylation","Complex III","K00411",
  "Oxidative phosphorylation","Complex III","K00412",
  "Oxidative phosphorylation","Complex III","K00410",
  "Oxidative phosphorylation","Complex V","K02117",
  "Oxidative phosphorylation","Complex V","K02118",
  "Oxidative phosphorylation","Complex V","K02111",
  "Oxidative phosphorylation","Complex V","K02112",
  
  "Oxygen metabolism","Cytochrome c oxidase, caa3-type","K02274",
  "Oxygen metabolism","Cytochrome c oxidase, caa3-type","K02275",
  "Oxygen metabolism","Cytochrome c oxidase, cbb3-type","K00404",
  "Oxygen metabolism","Cytochrome c oxidase, cbb3-type","K00405",
  "Oxygen metabolism","Cytochrome c oxidase, cbb3-type","K00406",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, bo type","K02297",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, bo type","K02300",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, bo type","K02298",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, bo type","K02299",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, bd type","K00425",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, bd type","K00426",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, aa3 type, QoxABCD","K02826",
  "Oxygen metabolism","Cytochrome (quinone) oxidase, aa3 type, QoxABCD","K02827",
  
  "Iron cycling","Iron oxidation","Cyc1.hmm",
  "Iron cycling","Iron oxidation","Cyc2_repCluster1.hmm",
  "Iron cycling","Iron oxidation","Cyc2_repCluster2.hmm",
  "Iron cycling","Iron oxidation","Cyc2_repCluster3.hmm",
  "Iron cycling","Iron oxidation","FoxA.hmm",
  "Iron cycling","Iron oxidation","FoxB.hmm",
  "Iron cycling","Iron oxidation","FoxC.hmm",
  "Iron cycling","Iron oxidation","FoxE.hmm",
  "Iron cycling","Iron oxidation","FoxY.hmm",
  "Iron cycling","Iron oxidation","FoxZ.hmm",
  "Iron cycling","Iron oxidation","MtoA.hmm",
  "Iron cycling","Iron oxidation","sulfocyanin.hmm",
  
  "Iron cycling","Iron reduction","CymA.hmm",
  "Iron cycling","Iron reduction","OmcF.hmm",
  "Iron cycling","Iron reduction","OmcS.hmm",
  "Iron cycling","Iron reduction","OmcZ.hmm",
  "Iron cycling","Iron reduction","FmnA.hmm",
  "Iron cycling","Iron reduction","DmkA.hmm",
  "Iron cycling","Iron reduction","FmnB.hmm",
  "Iron cycling","Iron reduction","PplA.hmm",
  "Iron cycling","Iron reduction","Ndh2.hmm",
  "Iron cycling","Iron reduction","EetA.hmm",
  "Iron cycling","Iron reduction","EetB.hmm",
  "Iron cycling","Iron reduction","DmkB.hmm",
  "Iron cycling","Iron reduction","DFE_0448.hmm",
  "Iron cycling","Iron reduction","DFE_0449.hmm",
  "Iron cycling","Iron reduction","DFE_0450.hmm",
  "Iron cycling","Iron reduction","DFE_0451.hmm",
  "Iron cycling","Iron reduction","DFE_0461.hmm",
  "Iron cycling","Iron reduction","DFE_0462.hmm",
  "Iron cycling","Iron reduction","DFE_0463.hmm",
  "Iron cycling","Iron reduction","DFE_0464.hmm",
  "Iron cycling","Iron reduction","DFE_0465.hmm",
  
  "Iron cycling","Metal (Iron/Manganese) reduction","MtrA.hmm",
  "Iron cycling","Metal (Iron/Manganese) reduction","MtrB_TIGR03509.hmm",
  "Iron cycling","Metal (Iron/Manganese) reduction","MtrC_TIGR03507.hmm",
  
  "Manganese cycling","Manganese oxidation","K06324",
  "Manganese cycling","Manganese oxidation","K22348",
  "Manganese cycling","Manganese oxidation","K22349",
  "Manganese cycling","Manganese oxidation","K22350",
  
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K00531",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K02586",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K02591",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K22896",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K22897",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K22898",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K02588",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K02567",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K02568",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K00370",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K00371",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K15876",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K03385",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K04015",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K00362",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K00363",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K00368",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K15864",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K04561",
  "Nitrogen cycling","Nitrogen cycling – red. pathways","K02305",
  
  "Nitrogen cycling","Nitrogen cycling – ox. pathways","K00370",
  "Nitrogen cycling","Nitrogen cycling – ox. pathways","K00371",
  "Nitrogen cycling","Nitrogen cycling – ox. pathways","K10944",
  "Nitrogen cycling","Nitrogen cycling – ox. pathways","K10945",
  "Nitrogen cycling","Nitrogen cycling – ox. pathways","K10946",
  
  "Sulfur cycling","Sulfur cycling – ox. pathways","K17229",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K17218",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K11180",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K11181",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K17725",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K16952",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K17224",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K17226",
  "Sulfur cycling","Sulfur cycling – ox. pathways","K17225",
  
  "Sulfur cycling","Sulfur cycling – red. pathways","K17219",
  "Sulfur cycling","Sulfur cycling – red. pathways","K17220",
  "Sulfur cycling","Sulfur cycling – red. pathways","K17221",
  "Sulfur cycling","Sulfur cycling – red. pathways","K16950",
  "Sulfur cycling","Sulfur cycling – red. pathways","K16951",
  "Sulfur cycling","Sulfur cycling – red. pathways","K00385",
  "Sulfur cycling","Sulfur cycling – red. pathways","K00394",
  "Sulfur cycling","Sulfur cycling – red. pathways","K00958",
  
  "Sulfur cycling","Choline-sulfatase","K01133"
)

# explode markers
marker_df <- marker_tbl %>%
  dplyr::mutate(Marker_ID_raw = str_split(Marker_ID_raw, "\\s*,\\s*")) %>%
  unnest(Marker_ID_raw) %>%
  dplyr::mutate(
    Marker_ID_raw = str_trim(Marker_ID_raw),
    is_hmm = str_detect(Marker_ID_raw, "\\.hmm$"),
    Marker_ID = dplyr::if_else(is_hmm, str_remove(Marker_ID_raw, "\\.hmm$"), Marker_ID_raw)
  ) %>%
  dplyr::select(Macro, Function, Marker_ID, is_hmm) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    Macro = factor(Macro, levels = macro_levels)
  )

# enforce category order respecting macro order + first appearance within macro
ordered_functions <- marker_tbl %>%
  dplyr::distinct(Macro, Function) %>%
  dplyr::mutate(Macro = factor(Macro, levels = macro_levels)) %>%
  dplyr::arrange(Macro) %>%
  dplyr::pull(Function) %>%
  unique()

marker_df <- marker_df %>%
  dplyr::mutate(Function = factor(Function, levels = ordered_functions))

# marker counts per category (for SUM/n_markers)
marker_counts <- marker_df %>%
  dplyr::distinct(Function, Marker_ID) %>%
  dplyr::count(Function, name = "n_markers")

needed_kegg <- marker_df %>% dplyr::filter(!is_hmm) %>% dplyr::pull(Marker_ID) %>% unique()
needed_hmm  <- marker_df %>% dplyr::filter(is_hmm)  %>% dplyr::pull(Marker_ID) %>% unique()

############################################################
### 4) Filter to needed markers + map Marker → Category + Macro
###     + build Depth_group + Depth_layer_plot
############################################################

df_fun <- df_MT %>%
  dplyr::mutate(source = dplyr::if_else(str_detect(Marker_ID, "^K\\d+"), "KEGG", "METAB")) %>%
  dplyr::filter(
    (source == "KEGG"  & Marker_ID %in% needed_kegg) |
      (source == "METAB" & Marker_ID %in% needed_hmm)
  ) %>%
  dplyr::left_join(marker_df %>% dplyr::select(Macro, Function, Marker_ID), by = "Marker_ID") %>%
  dplyr::filter(!is.na(Function), !is.na(Macro)) %>%
  dplyr::mutate(
    Macro    = factor(Macro, levels = macro_levels),
    Function = factor(Function, levels = ordered_functions),
    
    Depth_group = dplyr::if_else(as.character(Depth_layer) %in% shallow_levels, "Shallow", "Deep"),
    Depth_group = factor(Depth_group, levels = c("Shallow","Deep")),
    
    Depth_layer_plot = factor(as.character(Depth_layer), levels = y_levels)
  )

# depth combos really present AFTER marker-filtering
avail_PD <- df_fun %>%
  dplyr::distinct(Peatland, Depth_group, Depth_layer_plot) %>%
  dplyr::mutate(Peatland = factor(Peatland, levels = peatland_levels))

############################################################
### 5) Aggregate per cell (Peatland × Depth × Category)
###     SUM across samples + markers, then / n_markers
############################################################

plot_df <- df_fun %>%
  dplyr::group_by(Peatland, Depth_group, Depth_layer_plot, Macro, Function) %>%
  dplyr::summarise(value_sum = sum(Coverage_per_cell, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(marker_counts, by = "Function") %>%
  dplyr::mutate(value = value_sum / n_markers) %>%
  dplyr::select(-value_sum, -n_markers)

############################################################
### 6) Complete grid + FIXED activity levels
############################################################

thr_low  <- 0.01
thr_high <- 0.30

lab_zero <- "0"
lab_low  <- "0–0.01"
lab_med  <- "0.01–0.3"
lab_high <- "> 0.3"

activity_levels <- c(
  lab_zero,
  lab_low,
  lab_med,
  lab_high
)

plot_MT4 <- avail_PD %>%
  tidyr::crossing(
    Macro = factor(
      macro_levels,
      levels = macro_levels
    )
  ) %>%
  dplyr::left_join(
    marker_tbl %>%
      dplyr::distinct(
        .data$Macro,
        .data$Function
      ) %>%
      dplyr::mutate(
        Macro = factor(
          .data$Macro,
          levels = macro_levels
        ),
        Function = factor(
          .data$Function,
          levels = ordered_functions
        )
      ),
    by = "Macro"
  ) %>%
  dplyr::left_join(
    plot_df %>%
      dplyr::select(
        .data$Peatland,
        .data$Depth_group,
        .data$Depth_layer_plot,
        .data$Macro,
        .data$Function,
        .data$value
      ),
    by = c(
      "Peatland",
      "Depth_group",
      "Depth_layer_plot",
      "Macro",
      "Function"
    )
  ) %>%
  dplyr::mutate(
    value = replace_na(
      .data$value,
      0
    ),
    Activity_level = dplyr::case_when(
      .data$value == 0 ~ lab_zero,
      .data$value > 0 &
        .data$value <= thr_low ~ lab_low,
      .data$value > thr_low &
        .data$value <= thr_high ~ lab_med,
      .data$value > thr_high ~ lab_high,
      TRUE ~ NA_character_
    ),
    Activity_level = factor(
      .data$Activity_level,
      levels = activity_levels
    )
  )


############################################################
### 7) One common white-to-black scale
############################################################

activity_palette <- c(
  "0"        = "#F3F5F2",
  "0–0.01"   = "#D7DED5",
  "0.01–0.3" = "#A8B5A7",
  "> 0.3"    = "#6F7E70"
)
############################################################

readr::write_tsv(
  plot_MT4 %>%
    dplyr::transmute(
      Peatland = as.character(.data$Peatland),
      Depth_group = as.character(.data$Depth_group),
      Depth_layer = as.character(.data$Depth_layer_plot),
      Macro = as.character(.data$Macro),
      Function = as.character(.data$Function),
      MT_coverage_per_cell_category_value = .data$value,
      Activity_level = as.character(.data$Activity_level),
      Detectable_expression = .data$value > 0
    ),
  out_values
)

############################################################
### 9) Plot
###
### Macro-categories are separated by:
###   - separate column facets
###   - narrow white gaps
###   - black panel borders
###
### Fill represents only MT_coverage_per_cell activity.
############################################################

p_big <- ggplot(
  plot_MT4,
  aes(
    x = Function,
    y = Depth_layer_plot
  )
) +
  geom_tile(
    aes(
      fill = Activity_level
    ),
    colour = "grey15",
    linewidth = 0.28,
    width = 0.94,
    height = 0.94
  ) +
  facet_grid(
    rows = vars(
      Peatland,
      Depth_group
    ),
    cols = vars(
      Macro
    ),
    switch = "y",
    scales = "free",
    space = "free",
    drop = TRUE,
    labeller = labeller(
      Depth_group = function(x) {
        rep(
          "",
          length(x)
        )
      },
      Macro = label_wrap_gen(
        width = 18
      )
    )
  ) +
  scale_fill_manual(
    values = activity_palette,
    breaks = activity_levels,
    limits = activity_levels,
    name = "MT coverage per cell\n(contig level)",
    drop = FALSE,
    na.value = "white"
  ) +
  scale_x_discrete(
    labels = function(x) {
      stringr::str_wrap(
        x,
        width = 18
      )
    },
    position = "top",
    drop = TRUE
  ) +
  scale_y_discrete(
    drop = TRUE,
    labels = function(x) {
      stringr::str_replace_all(
        x,
        "\\s*cm",
        ""
      )
    }
  ) +
  theme_bw(
    base_size = 12,
    base_family = "Helvetica"
  ) +
  theme(
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    panel.spacing.x = grid::unit(
      3.0,
      "pt"
    ),
    
    panel.spacing.y = grid::unit(
      gap_pt,
      "pt"
    ),
    
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    
    strip.placement = "outside",
    
    strip.background.x = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.35
    ),
    
    strip.background.y = element_blank(),
    
    strip.text.x = element_text(
      face = "bold",
      size = 8,
      colour = "black",
      margin = margin(
        t = 3,
        r = 2,
        b = 3,
        l = 2
      )
    ),
    
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      colour = "black"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key = element_rect(
      colour = "grey15",
      linewidth = 0.3
    ),
    
    legend.key.height = grid::unit(
      13,
      "pt"
    ),
    
    legend.key.width = grid::unit(
      13,
      "pt"
    ),
    
    plot.margin = grid::unit(
      c(
        10,
        15,
        10,
        10
      ),
      "pt"
    )
  ) +
  labs(
    x = "Metabolic function",
    y = "Depth layer"
  )

p_clear <- p_big +
  theme(
    legend.position = "none"
  )

print(p_big)

############################################################
### 10) Save outputs
############################################################

ggsave(
  filename = out_svg,
  plot = p_big,
  device = svglite::svglite,
  width = 14,
  height = 8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = out_svg_no_legend,
  plot = p_clear,
  device = svglite::svglite,
  width = 12.5,
  height = 8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

if (capabilities("cairo")) {
  
  ggsave(
    filename = out_pdf,
    plot = p_big,
    device = grDevices::cairo_pdf,
    width = 14,
    height = 8,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  
} else {
  
  warning(
    "PDF was not saved because Cairo support is unavailable."
  )
}

cat("\nSaved heatmap with legend:\n")
cat(out_svg, "\n")

cat("\nSaved heatmap without legend:\n")
cat(out_svg_no_legend, "\n")

if (capabilities("cairo")) {
  cat("\nSaved heatmap PDF:\n")
  cat(out_pdf, "\n")
}

cat("\nSaved heatmap value table:\n")
cat(out_values, "\n")

cat(
  paste0(
    "\nDONE. Figure 2A uses MT_coverage_per_cell ",
    "and one common white-to-black activity scale.\n"
  )
)
