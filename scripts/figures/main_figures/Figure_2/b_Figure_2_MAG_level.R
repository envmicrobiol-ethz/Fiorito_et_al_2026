#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Figure 2B: dominant-phylum bubbles for metatranscriptomic activity
# across metabolic functions, peatlands and depth layers.
#
# Marker-level MT_coverage_per_cell values are summed within each MAG and
# normalized by the number of markers assigned to the corresponding function.
# Values are then summed by phylum. The dominant globally abundant displayed
# phylum is selected for each heatmap cell, while bubble size represents one
# of three fixed activity classes.
#
# INPUT
# 1. MAG taxonomy TSV containing genome and Phylum columns.
# 2. MAG-level KEGG expression-profile TSV.
# 3. MAG-level METABOLIC expression-profile TSV.
# 4. Output directory.
#
# OUTPUT
# Figure 2B as SVG and PDF, versions without legends, and the plotted values TSV.
#
# USAGE
# Rscript Figure_2B.R \
#   MAG_taxonomy.tsv \
#   mags.KEGG.expression_profile.tsv \
#   mags.METABOLIC.expression_profile.tsv \
#   figure_output_directory

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript Figure_2B.R",
      "<MAG_taxonomy.tsv>",
      "<KEGG_expression.tsv>",
      "<METABOLIC_expression.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

tax_file <- args[[1]]
kegg_file <- args[[2]]
metab_file <- args[[3]]
output_dir <- args[[4]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg_bubble <- file.path(
  output_dir,
  "Figure_2B_dominant_phylum_bubbles.svg"
)

out_svg_bubble_no_legend <- file.path(
  output_dir,
  "Figure_2B_dominant_phylum_bubbles_no_legend.svg"
)

out_pdf_bubble <- file.path(
  output_dir,
  "Figure_2B_dominant_phylum_bubbles.pdf"
)

out_pdf_bubble_no_legend <- file.path(
  output_dir,
  "Figure_2B_dominant_phylum_bubbles_no_legend.pdf"
)

out_tsv <- file.path(
  output_dir,
  "Figure_2B_dominant_phylum_bubble_values.tsv"
)


############################################################
### 0b) SETTINGS
############################################################

peatland_levels <- c("LM","NR","HV","BM")

depth_levels <- c(
  "0–10 cm","10–30 cm","30–60 cm",
  "200–225 cm","265–280 cm","325–350 cm","525–540 cm","525–550 cm","650–675 cm"
)

# ===== Y AXIS / GAP CONTROL (same logic as your other script) =====
shallow_levels <- c("0–10 cm","10–30 cm","30–60 cm")

# order within shallow (0–10 on top)
shallow_y_levels <- c("30–60 cm","10–30 cm","0–10 cm")

# order within deep (200–225 on top of deep block; deepest at bottom)
deep_y_levels <- c("650–675 cm","525–550 cm","525–540 cm","325–350 cm","265–280 cm","200–225 cm")

# global y levels (deep block below + shallow block above)
y_levels <- c(deep_y_levels, shallow_y_levels)

# gap between shallow and deep + between peatlands (pt)
gap_pt <- 1.2

# alpha per i 3 livelli
alpha_low  <- 0
alpha_med  <- 0.45
alpha_high <- 1

# top N phyla globali
topN_phyla <- 15

# thresholds fissi
thr_low  <- 0.01
thr_high <- 0.30

############################################################
### helper: SUM che non trasforma "tutto NA" in 0
############################################################
safe_sum <- function(x){
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

fmt <- function(x) formatC(x, format = "fg", digits = 3)

############################################################

############################################################
### 1) TAXONOMY
############################################################

tax <- readr::read_tsv(tax_file, show_col_types = FALSE) %>%
  dplyr::rename(MAG_name = genome) %>%
  dplyr::select(MAG_name, Phylum) %>%
  dplyr::filter(!is.na(MAG_name), !is.na(Phylum))

############################################################
### 2) CATEGORIES / MARKERS (KOs + HMMs)
############################################################

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

marker_tbl <- tibble::tribble(
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


marker_df <- marker_tbl %>%
  dplyr::mutate(
    Marker_ID_raw = stringr::str_split(
      .data$Marker_ID_raw,
      "\\s*,\\s*"
    )
  ) %>%
  tidyr::unnest(cols = dplyr::all_of("Marker_ID_raw")) %>%
  dplyr::mutate(
    Marker_ID_raw = stringr::str_trim(.data$Marker_ID_raw),
    is_hmm = stringr::str_detect(.data$Marker_ID_raw, "\\.hmm$"),
    Marker_ID = dplyr::if_else(
      .data$is_hmm,
      stringr::str_remove(.data$Marker_ID_raw, "\\.hmm$"),
      .data$Marker_ID_raw
    ),
    Macro = factor(.data$Macro, levels = macro_levels)
  ) %>%
  dplyr::select(
    .data$Macro,
    .data$Function,
    .data$Marker_ID,
    .data$is_hmm
  ) %>%
  dplyr::distinct()

# Function order is identical to the contig-level heatmap:
# macro order first, then first appearance within each macro-category.
ordered_functions <- marker_tbl %>%
  dplyr::distinct(.data$Macro, .data$Function) %>%
  dplyr::mutate(
    Macro = factor(.data$Macro, levels = macro_levels)
  ) %>%
  dplyr::arrange(.data$Macro) %>%
  dplyr::pull(.data$Function) %>%
  unique()

marker_df <- marker_df %>%
  dplyr::mutate(
    Function = factor(.data$Function, levels = ordered_functions)
  )

marker_counts <- marker_df %>%
  dplyr::distinct(.data$Function, .data$Marker_ID) %>%
  dplyr::count(.data$Function, name = "n_markers")

needed_kegg <- marker_df %>%
  dplyr::filter(!.data$is_hmm) %>%
  dplyr::pull(.data$Marker_ID) %>%
  unique()

needed_hmm <- marker_df %>%
  dplyr::filter(.data$is_hmm) %>%
  dplyr::pull(.data$Marker_ID) %>%
  unique()


############################################################
### 3) HELPERS: peatland + depth parsing
############################################################

clean_depth <- function(x){
  x %>%
    stringr::str_replace_all("–", "-") %>%
    stringr::str_replace_all("—", "-") %>%
    stringr::str_replace_all(" ", "") %>%
    stringr::str_to_lower()
}

depth_map <- tibble::tribble(
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


############################################################
### 4) LOAD EXPRESSION TABLES (KEGG + METABOLIC) and MERGE
############################################################

kegg <- readr::read_tsv(kegg_file, show_col_types = FALSE) %>%
  dplyr::transmute(
    Marker_ID = as.character(KEGG_ko),
    MetaT_sample = as.character(MetaT_sample),
    MAG_name = as.character(MAG_name),
    MT_coverage_per_cell = suppressWarnings(as.numeric(MT_coverage_per_cell))
  ) %>%
  dplyr::filter(Marker_ID %in% needed_kegg) %>%
  dplyr::filter(!is.na(MAG_name), !is.na(MetaT_sample), !is.na(Marker_ID))

metab <- readr::read_tsv(metab_file, show_col_types = FALSE) %>%
  dplyr::transmute(
    Marker_ID = as.character(METABOLIC_hmm),
    MetaT_sample = as.character(MetaT_sample),
    MAG_name = as.character(MAG_name),
    MT_coverage_per_cell = suppressWarnings(as.numeric(MT_coverage_per_cell))
  ) %>%
  dplyr::filter(Marker_ID %in% needed_hmm) %>%
  dplyr::filter(!is.na(MAG_name), !is.na(MetaT_sample), !is.na(Marker_ID))

df <- dplyr::bind_rows(kegg, metab) %>%
  dplyr::filter(!is.na(.data$MT_coverage_per_cell))


############################################################
############################################################
### 5) JOIN taxonomy + parse Peatland/Depth + map Marker→Function
############################################################

df <- df %>%
  dplyr::left_join(tax, by = "MAG_name") %>%
  dplyr::filter(!is.na(Phylum)) %>%
  dplyr::mutate(
    Peatland = dplyr::case_when(
      stringr::str_detect(MetaT_sample, stringr::regex("bj(ö|o)rsmossen|bjorsmossen", ignore_case = TRUE)) ~ "BM",
      stringr::str_detect(MetaT_sample, stringr::regex("norra[_-]?romyren|norraromyren|romyren", ignore_case = TRUE)) ~ "NR",
      stringr::str_detect(MetaT_sample, stringr::regex("havsj(ö|o)mossen|havsjomossen", ignore_case = TRUE)) ~ "HV",
      stringr::str_detect(MetaT_sample, stringr::regex("lungsmossen|lungs?mossen", ignore_case = TRUE)) ~ "LM",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Peatland)) %>%
  dplyr::mutate(
    Depth_raw = clean_depth(
      dplyr::coalesce(
        stringr::str_extract(MetaT_sample, "\\d+-\\d+cm|\\d+cm"),
        stringr::str_extract(MAG_name,     "\\d+-\\d+cm|\\d+cm")
      )
    )
  ) %>%
  dplyr::left_join(depth_map, by = "Depth_raw") %>%

  dplyr::mutate(
    Depth_layer = factor(
      .data$Depth_layer,
      levels = depth_levels
    ),
    Peatland = factor(
      .data$Peatland,
      levels = peatland_levels
    )
  ) %>%

  dplyr::left_join(
    marker_df %>%
      dplyr::select(.data$Macro, .data$Function, .data$Marker_ID),
    by = "Marker_ID"
  ) %>%
  dplyr::filter(
    !is.na(.data$Macro),
    !is.na(.data$Function),
    !is.na(.data$Depth_layer)
  ) %>%
  dplyr::mutate(
    Macro = factor(.data$Macro, levels = macro_levels),
    Function = factor(.data$Function, levels = ordered_functions)
  )

# ===== Depth_group + Depth_layer_plot (needed for faceting + exact Y order) =====
df <- df %>%
  dplyr::mutate(
    Depth_group = dplyr::if_else(as.character(Depth_layer) %in% shallow_levels, "Shallow", "Deep"),
    Depth_group = factor(Depth_group, levels = c("Shallow","Deep")),
    Depth_layer_plot = factor(as.character(Depth_layer), levels = y_levels)
  )

# combos really present AFTER all filters
avail_PD <- df %>%
  dplyr::distinct(Peatland, Depth_group, Depth_layer_plot) %>%
  dplyr::mutate(Peatland = factor(Peatland, levels = peatland_levels))

############################################################
### 6) CATEGORY EXPRESSION PER MAG PER SAMPLE
###     markers within category: SUM, poi NORMALIZZA per n_markers (SUM/n_markers)
############################################################

mag_fun <- df %>%
  dplyr::group_by(MAG_name, MetaT_sample, Peatland, Depth_layer, Depth_group, Depth_layer_plot, Macro, Function, Phylum) %>%
  dplyr::summarise(expr_sum = safe_sum(.data$MT_coverage_per_cell), .groups = "drop") %>%
  dplyr::left_join(marker_counts, by = "Function") %>%
  dplyr::mutate(expr = expr_sum / n_markers) %>%
  dplyr::select(-expr_sum, -n_markers)

############################################################
### 7) PHYLUM EXPRESSION PER SAMPLE
############################################################

phy_fun_sample <- mag_fun %>%
  dplyr::group_by(Peatland, Depth_layer, Depth_group, Depth_layer_plot, MetaT_sample, Macro, Function, Phylum) %>%
  dplyr::summarise(expr = safe_sum(expr), .groups = "drop")

phy_fun <- phy_fun_sample %>%
  dplyr::group_by(Peatland, Depth_layer, Depth_group, Depth_layer_plot, Macro, Function, Phylum) %>%
  dplyr::summarise(expr = safe_sum(expr), .groups = "drop") %>%
  dplyr::mutate(expr = dplyr::if_else(is.finite(expr), expr, NA_real_))

############################################################
### 8) TOP N PHYLA GLOBALLY + WINNER per cell
############################################################

top_phyla <- phy_fun %>%
  dplyr::filter(!is.na(expr), expr > 0) %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(total = sum(expr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::slice_head(n = topN_phyla) %>%
  dplyr::pull(Phylum)

phy_fun_top <- phy_fun %>%
  dplyr::filter(Phylum %in% top_phyla, !is.na(expr), expr > 0)

# ===== winner keeps Depth_group + Depth_layer_plot so plotting join is stable =====
winner <- phy_fun_top %>%
  dplyr::group_by(Peatland, Depth_group, Depth_layer_plot, Macro, Function) %>%
  dplyr::slice_max(order_by = expr, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

############################################################
### 9) COMPLETE GRID (only existing Peatland×Depth) + LEVELS
############################################################

function_layout <- marker_tbl %>%
  dplyr::distinct(.data$Macro, .data$Function) %>%
  dplyr::mutate(
    Macro = factor(.data$Macro, levels = macro_levels),
    Function = factor(.data$Function, levels = ordered_functions)
  )

grid <- avail_PD %>%
  tidyr::crossing(
    Macro = factor(macro_levels, levels = macro_levels)
  ) %>%
  dplyr::left_join(
    function_layout,
    by = "Macro"
  )

plot_df <- grid %>%
  dplyr::left_join(
    winner,
    by = c(
      "Peatland",
      "Depth_group",
      "Depth_layer_plot",
      "Macro",
      "Function"
    )
  ) %>%
  dplyr::mutate(
    Macro = factor(.data$Macro, levels = macro_levels),
    Function = factor(.data$Function, levels = ordered_functions)
  )

lab_low  <- paste0("≤ ", fmt(thr_low))
lab_med  <- paste0(fmt(thr_low), " – ", fmt(thr_high))
lab_high <- paste0("> ", fmt(thr_high))

plot_df <- plot_df %>%
  dplyr::mutate(
    Level = dplyr::case_when(
      is.na(expr) | expr <= thr_low  ~ lab_low,
      expr <= thr_high              ~ lab_med,
      TRUE                          ~ lab_high
    ),
    Level = factor(Level, levels = c(lab_low, lab_med, lab_high))
  )

############################################################
### 10) COLORS FOR TOP PHYLA (manual palette + fallback)
###     + legend order EXACTLY as your vectors (bact then arch)
############################################################

bact_colors <- c(
  "p__Acidobacteriota"    = "#B2182B",
  "p__Actinomycetota"     = "#FFF04F",
  "p__Bacteroidota"       = "#FE8CA1",
  "p__Chloroflexota"      = "darkorange",
  "p__Desulfobacterota"   = "#C77526",
  "p__Desulfobacterota_G" = "#DD1C77",
  "p__FCPU426"            = "violet",
  "p__Planctomycetota"    = "#901394",
  "p__Pseudomonadota"     = "#8C510A",
  "p__Verrucomicrobiota"  = "#D8B365"
)

arch_colors <- c(
  "p__Thermoproteota" = "#08306B",
  "p__Halobacteriota" = "#2171B5"
)

# Base palette. Archaeal and bacterial colours are retained from the original figure.
phylum_colors_manual <- c(bact_colors, arch_colors)

# Add fallback colours only for top phyla not present in the manual palette.
missing_top_phyla <- setdiff(top_phyla, names(phylum_colors_manual))

if (length(missing_top_phyla) > 0) {
  phylum_colors_manual <- c(
    phylum_colors_manual,
    stats::setNames(
      scales::hue_pal()(length(missing_top_phyla)),
      missing_top_phyla
    )
  )
}

# Show only the actual global top-N phyla in the legend.
legend_phyla_order <- c(
  intersect(names(bact_colors), top_phyla),
  intersect(names(arch_colors), top_phyla),
  setdiff(top_phyla, c(names(bact_colors), names(arch_colors)))
)

phylum_colors <- phylum_colors_manual[legend_phyla_order]

plot_df <- plot_df %>%
  dplyr::mutate(
    Phylum = factor(.data$Phylum, levels = legend_phyla_order)
  )

############################################################
### 11) COMMON PLOT GEOMETRY
###
### This geometry matches the contig-level heatmap:
###   rows = Peatland × Depth_group
###   columns = Macro
###   free x-space = separate metabolic sectors
###
### Use the no-legend PDF/SVG for exact panel alignment.
############################################################

common_facet <- ggplot2::facet_grid(
  rows = ggplot2::vars(
    Peatland,
    Depth_group
  ),
  cols = ggplot2::vars(
    Macro
  ),
  switch = "y",
  scales = "free",
  space = "free",
  drop = TRUE,
  labeller = ggplot2::labeller(
    Depth_group = function(x) {
      rep("", length(x))
    },
    Macro = ggplot2::label_wrap_gen(
      width = 18
    )
  )
)

common_theme <- ggplot2::theme_bw(
  base_size = 12,
  base_family = "Arial"
) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    
    panel.border = ggplot2::element_rect(
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
    
    axis.text.x = ggplot2::element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8.5,
      colour = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 9,
      colour = "black"
    ),
    
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    
    strip.placement = "outside",
    
    strip.background.x = ggplot2::element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.35
    ),
    
    strip.background.y = ggplot2::element_blank(),
    
    strip.text.x = ggplot2::element_text(
      face = "bold",
      size = 8,
      colour = "black",
      margin = ggplot2::margin(
        t = 3,
        r = 2,
        b = 3,
        l = 2
      )
    ),
    
    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold",
      colour = "black"
    ),
    
    legend.position = "right",
    legend.box = "vertical",
    
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = ggplot2::element_text(
      size = 9
    ),
    
    plot.margin = grid::unit(
      c(10, 15, 10, 10),
      "pt"
    )
  )

common_x_scale <- ggplot2::scale_x_discrete(
  labels = function(x) {
    stringr::str_wrap(
      x,
      width = 18
    )
  },
  position = "top",
  drop = TRUE
)

common_y_scale <- ggplot2::scale_y_discrete(
  drop = TRUE,
  labels = function(x) {
    stringr::str_replace_all(
      x,
      "\\s*cm",
      ""
    )
  }
)


############################################################
### 13) BUBBLE PLOT
############################################################

plot_df2 <- plot_df %>%
  dplyr::mutate(
    Level = dplyr::case_when(
      is.na(.data$expr) ~ NA_character_,
      .data$expr <= thr_low ~ lab_low,
      .data$expr <= thr_high ~ lab_med,
      TRUE ~ lab_high
    ),
    Level = factor(
      .data$Level,
      levels = c(lab_low, lab_med, lab_high)
    )
  )

size_vals <- stats::setNames(
  c(1.8, 3.4, 5.2),
  c(lab_low, lab_med, lab_high)
)

p_bubble <- ggplot2::ggplot(
  plot_df2,
  ggplot2::aes(
    x = .data$Function,
    y = .data$Depth_layer_plot
  )
) +
  ggplot2::geom_tile(
    fill = "white",
    colour = "grey55",
    linewidth = 0.25,
    width = 0.94,
    height = 0.94
  ) +
  ggplot2::geom_point(
    data = plot_df2 %>%
      dplyr::filter(
        !is.na(.data$Level),
        !is.na(.data$Phylum)
      ),
    ggplot2::aes(
      fill = .data$Phylum,
      size = .data$Level
    ),
    shape = 21,
    colour = "grey20",
    stroke = 0.25
  ) +
  common_facet +
  ggplot2::scale_fill_manual(
    values = phylum_colors,
    breaks = legend_phyla_order,
    limits = legend_phyla_order,
    drop = FALSE,
    name = "Phylum (top global)",
    labels = function(x) {
      stringr::str_remove(x, "^p__")
    },
    na.value = "white",
    na.translate = FALSE
  ) +
  ggplot2::scale_size_manual(
    values = size_vals,
    name = "MT coverage per cell\n(sum/n_markers)",
    drop = FALSE
  ) +
  common_x_scale +
  common_y_scale +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        size = 4
      )
    ),
    size = ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        fill = "grey70"
      )
    )
  ) +
  common_theme +
  ggplot2::labs(
    x = "Metabolic function",
    y = "Depth layer"
  )

p_bubble_no_legend <- p_bubble +
  ggplot2::theme(
    legend.position = "none"
  )

print(p_bubble)

############################################################
### 14) EXPORT VALUES
############################################################

bubble_values <- plot_df2 %>%
  dplyr::transmute(
    Peatland = as.character(.data$Peatland),
    Depth_group = as.character(.data$Depth_group),
    Depth_layer = as.character(.data$Depth_layer_plot),
    Macro = as.character(.data$Macro),
    Function = as.character(.data$Function),
    Dominant_top15_phylum = as.character(.data$Phylum),
    MT_coverage_per_cell = .data$expr,
    MT_coverage_per_cell_bin = as.character(.data$Level),
    Bubble_shown = !is.na(.data$Phylum) & !is.na(.data$expr)
  )

readr::write_tsv(
  bubble_values,
  out_tsv
)

############################################################
### 15) SAVE OUTPUTS
###
### Dimensions match the neutral contig-level heatmap:
###   with legend:    14 × 8 inches
###   without legend: 12.5 × 8 inches
###
### For exact vertical stacking/alignment, use the two
### no-legend files, both exported at 12.5 × 8 inches.
############################################################



ggplot2::ggsave(
  filename = out_svg_bubble,
  plot = p_bubble,
  device = svglite::svglite,
  width = 14,
  height = 8,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)

ggplot2::ggsave(
  filename = out_svg_bubble_no_legend,
  plot = p_bubble_no_legend,
  device = svglite::svglite,
  width = 12.5,
  height = 8,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)

if (capabilities("cairo")) {
  
  ggplot2::ggsave(
    filename = out_pdf_bubble,
    plot = p_bubble,
    device = grDevices::cairo_pdf,
    width = 14,
    height = 8,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = out_pdf_bubble_no_legend,
    plot = p_bubble_no_legend,
    device = grDevices::cairo_pdf,
    width = 12.5,
    height = 8,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )
  
} else {
  
  warning(
    "PDF files were not saved because Cairo support is unavailable."
  )
}



cat("\nSaved bubble plot with legend:\n")
cat(out_svg_bubble, "\n")

cat("\nSaved bubble plot without legend:\n")
cat(out_svg_bubble_no_legend, "\n")

if (capabilities("cairo")) {
  cat("\nSaved bubble PDF with legend:\n")
  cat(out_pdf_bubble, "\n")
  
  cat("\nSaved bubble PDF without legend:\n")
  cat(out_pdf_bubble_no_legend, "\n")
}

cat("\nSaved bubble values table:\n")
cat(out_tsv, "\n")

cat(
  paste0(
    "\nDONE. MT_coverage_per_cell was used throughout; ",
    "metabolic macro-categories are separated exactly as in the ",
    "contig-level heatmap.\n"
  )
)
