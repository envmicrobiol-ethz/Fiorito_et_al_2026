#!/usr/bin/env Rscript

# DESCRIPTION
# Generates Extended Data Figure 3, showing the dominant phylum associated with expression
# of terminal reductases and hydrogenases across peatlands and
# depth layers.
#
# MT_coverage_per_cell values are summed for each MAG, sample and gene and
# divided by the number of markers assigned to that gene. Values are then
# summed by phylum, and the true dominant phylum is selected before reducing
# the legend to eight explicitly displayed phyla plus Other phyla.
#
# INPUT
# 1. MAG taxonomy TSV containing genome and Phylum.
# 2. MAG-level KEGG expression-profile TSV.
# 3. MAG-level METABOLIC expression-profile TSV.
# 4. MAG-level PFAM expression-profile TSV.
# 5. Output directory.
#
# OUTPUT
# Extended Data Figure 3 as SVG and PDF, plotted values TSV and phylum-ranking TSV.
#
# USAGE
# Rscript ED_Figure_3.R \
#   MAG_taxonomy.tsv \
#   MAGs.KEGG.expression_profile.tsv \
#   MAGs.METABOLIC.expression_profile.tsv \
#   MAGs.PFAM.expression_profile.tsv \
#   output_directory

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

if (length(args) != 5) {
  stop(
    paste(
      "Usage: Rscript ED_Figure_3.R",
      "<MAG_taxonomy.tsv>",
      "<KEGG_expression.tsv>",
      "<METABOLIC_expression.tsv>",
      "<PFAM_expression.tsv>",
      "<output_directory>"
    ),
    call. = FALSE
  )
}

tax_file <- args[[1]]
kegg_file <- args[[2]]
metab_file <- args[[3]]
pfam_file <- args[[4]]
output_dir <- args[[5]]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_svg <- file.path(
  output_dir,
  "ED_Figure_3_terminal_reductases_EET_hydrogenases.svg"
)

out_pdf <- file.path(
  output_dir,
  "ED_Figure_3_terminal_reductases_EET_hydrogenases.pdf"
)

out_values <- file.path(
  output_dir,
  "ED_Figure_3_terminal_reductases_EET_hydrogenases_values.tsv"
)

out_ranking <- file.path(
  output_dir,
  "ED_Figure_3_phylum_expression_ranking.tsv"
)

plot_font <- "Helvetica"

peatland_levels <- c("LM", "NR", "HV", "BM")

depth_levels <- c(
  "0–10 cm", "10–30 cm", "30–60 cm",
  "200–225 cm", "265–280 cm", "325–350 cm",
  "525–540 cm", "525–550 cm", "650–675 cm"
)

shallow_levels <- c(
  "0–10 cm",
  "10–30 cm",
  "30–60 cm"
)

shallow_y_levels <- c(
  "30–60 cm",
  "10–30 cm",
  "0–10 cm"
)

deep_y_levels <- c(
  "650–675 cm",
  "525–550 cm",
  "525–540 cm",
  "325–350 cm",
  "265–280 cm",
  "200–225 cm"
)

y_levels <- c(
  deep_y_levels,
  shallow_y_levels
)

gap_pt <- 1.2

thr_low <- 0.01
thr_high <- 0.30

selected_phyla <- c(
  "p__Halobacteriota",
  "p__Thermoproteota",
  "p__Acidobacteriota",
  "p__Actinomycetota",
  "p__Verrucomicrobiota",
  "p__Chloroflexota",
  "p__Desulfobacterota",
  "p__Pseudomonadota"
)

phylum_levels <- c(
  selected_phyla,
  "Other phyla"
)

phylum_colors <- c(
  "p__Halobacteriota"    = "#2171B5",
  "p__Thermoproteota"    = "#08306B",
  "p__Acidobacteriota"   = "#B2182B",
  "p__Actinomycetota"    = "#FFF04F",
  "p__Verrucomicrobiota" = "#D8B365",
  "p__Chloroflexota"     = "#FF8C00",
  "p__Desulfobacterota"  = "#C77526",
  "p__Pseudomonadota"    = "#8C510A",
  "Other phyla"          = "grey60"
)

safe_sum <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  sum(
    x,
    na.rm = TRUE
  )
}

fmt <- function(x) {
  formatC(
    x,
    format = "fg",
    digits = 3
  )
}

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
  "200-225cm",  "200–225 cm",
  "265-280cm",  "265–280 cm",
  "325-350cm",  "325–350 cm",
  "525-540cm",  "525–540 cm",
  "525-550cm",  "525–550 cm",
  "650-675cm",  "650–675 cm"
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


############################################################
### 2) TAXONOMY
############################################################

tax <- readr::read_tsv(
  tax_file,
  show_col_types = FALSE
) |>
  dplyr::rename(
    MAG_name = .data$genome
  ) |>
  dplyr::select(
    .data$MAG_name,
    .data$Phylum
  ) |>
  dplyr::filter(
    !is.na(.data$MAG_name),
    !is.na(.data$Phylum)
  )

############################################################
### 3) GENE CATALOG
############################################################

gene_tbl <- tribble(
  ~Gene,        ~Description,                                                      ~Marker_ID,
  "coxA",       "cytochrome c oxidase subunit I",                                  "K02274",
  "coxB",       "cytochrome c oxidase subunit II",                                 "K02275",
  "ccoN",       "cytochrome c oxidase cbb3-type subunit I",                         "K00404",
  "ccoO",       "cytochrome c oxidase cbb3-type subunit II",                        "K00405",
  "ccoP",       "cytochrome c oxidase cbb3-type subunit III",                       "K00406",
  "cyoA",       "cytochrome o ubiquinol oxidase subunit II",                        "K02297",
  "cyoD",       "cytochrome o ubiquinol oxidase subunit IV",                        "K02300",
  "cyoB",       "cytochrome o ubiquinol oxidase subunit I",                         "K02298",
  "cyoC",       "cytochrome o ubiquinol oxidase subunit III",                       "K02299",
  "cydA",       "cytochrome bd ubiquinol oxidase subunit I",                        "K00425",
  "cydB",       "cytochrome bd ubiquinol oxidase subunit II",                       "K00426",
  "bcrC",       "benzoyl-coa reductase subunit C",                                  "K04112",
  "bcrB",       "benzoyl-coa reductase subunit B",                                  "K04113",
  "bcrA",       "benzoyl-coa reductase subunit A",                                  "K04114",
  "bcrD",       "benzoyl-coa reductase subunit D",                                  "K04115",
  "nxrA",       "nitrate reductase / nitrite oxidoreductase, alpha subunit",        "K00370",
  "nxrB",       "nitrate reductase / nitrite oxidoreductase, beta subunit",         "K00371",
  "napA",       "periplasmic nitrate reductase napa",                               "K02567",
  "napB",       "cytochrome c-type protein napb",                                   "K02568",
  "narG",       "nitrate reductase / nitrite oxidoreductase, alpha subunit",        "K00370",
  "narH",       "nitrate reductase / nitrite oxidoreductase, beta subunit",         "K00371",
  "nrfH",       "cytochrome c nitrite reductase small subunit",                      "K15876",
  "nrfA",       "nitrite reductase (cytochrome c-552)",                             "K03385",
  "nrfD",       "cytochrome c nitrite reductase, nrfd subunit",                      "K04015",
  "nirB",       "nitrite reductase (nadh) large subunit",                           "K00362",
  "nirD",       "nitrite reductase (nadh) small subunit",                           "K00363",
  "nirK",       "nitrite reductase (no-forming)",                                   "K00368",
  "nirS",       "nitrite reductase (no-forming) / hydroxylamine reductase",         "K15864",
  "norB",       "nitric oxide reductase subunit B",                                 "K04561",
  "norC",       "nitric oxide reductase subunit C",                                 "K02305",
  "nosD",       "nitrous oxidase accessory protein",                                "PF05048",
  "nosZ",       "nitrous-oxide reductase",                                          "K00376",
  "dsrD",       "dissimilatory sulfite reductase delta subunit",                    "PF08679",
  "dsrA",       "dissimilatory sulfite reductase alpha subunit",                    "K11180",
  "dsrB",       "dissimilatory sulfite reductase beta subunit",                     "K11181",
  "sreA",       "sulfur reductase molybdopterin subunit",                           "K17219",
  "sreB",       "sulfur reductase fes subunit",                                     "K17220",
  "sreC",       "sulfur reductase membrane anchor",                                 "K17221",
  "asrA",       "anaerobic sulfite reductase subunit A",                            "K16950",
  "asrB",       "anaerobic sulfite reductase subunit B",                            "K16951",
  "asrC",       "anaerobic sulfite reductase subunit C",                            "K00385",
  "aprA",       "adenylylsulfate reductase, subunit A",                             "K00394",
  "dmsA",       "anaerobic dimethyl sulfoxide reductase subunit A",                 "K07306",
  "phsA",       "thiosulfate reductase / polysulfide reductase chain A",            "K08352",
  "ClrA",       "chlorate reductase subunit alpha",                                 "K17050",
  "ClrB",       "chlorate reductase subunit gamma",                                 "K17051",
  "ClrC",       "chlorate reductase subunit beta",                                  "K17052",
  "arrA",       "dissimilatory arsenate reductase",                                 "K28466",
  "arsC (grx)", "arsenate reductase (glutaredoxin)",                                "K00537",
  "arsC (trx)", "arsenate reductase (thioredoxin)",                                 "K03741",
  "ygfM",       "putative selenate reductase fad-binding subunit",                  "K12529",
  "xdhD",       "putative selenate reductase molybdopterin-binding subunit",        "K12528",
  "YgfK",       "putative selenate reductase",                                      "K12527",
  "frdA",       "fumarate reductase subunit A",                                     "K00244",
  "frdB",       "fumarate reductase subunit B",                                     "K00245",
  "frdC",       "fumarate reductase subunit C",                                     "K00246",
  "frdD",       "fumarate reductase subunit D",                                     "K00247"
) |>
  dplyr::mutate(
    Description = dplyr::if_else(
      is.na(Description) | Description == "",
      "",
      str_replace(Description, "^[[:alpha:]]", toupper)
    )
  )

new_metab_genes_raw <- c(
  "CymA.hmm","OmcF.hmm","OmcS.hmm","OmcZ.hmm",
  "FmnA.hmm","DmkA.hmm","FmnB.hmm","PplA.hmm","Ndh2.hmm",
  "EetA.hmm","EetB.hmm","DmkB.hmm",
  "DFE_0448.hmm","DFE_0449.hmm","DFE_0450.hmm","DFE_0451.hmm",
  "DFE_0461.hmm","DFE_0462.hmm","DFE_0463.hmm","DFE_0465.hmm",
  "MtrA.hmm","MtrB_TIGR03509.hmm","MtrC_TIGR03507.hmm",
  "fefe-group-a13.hmm","fefe-group-a2.hmm","fefe-group-a4.hmm",
  "fefe-group-b.hmm","fefe-group-c1.hmm","fefe-group-c2.hmm","fefe-group-c3.hmm",
  "nife-group-1.hmm","nife-group-2ade.hmm","nife-group-2bc.hmm",
  "nife-group-3abd.hmm","nife-group-3c.hmm","nife-group-4a-g.hmm","nife-group-4hi.hmm"
)

new_metab_genes <- stringr::str_remove(new_metab_genes_raw, "\\.hmm$")

new_metab_tbl <- tibble(
  Gene        = new_metab_genes,
  Description = "",
  Marker_ID   = new_metab_genes
)

gene_tbl <- dplyr::bind_rows(gene_tbl, new_metab_tbl) |>
  dplyr::distinct(.data$Gene, .data$Marker_ID, .keep_all = TRUE)

gene_order <- gene_tbl$Gene

marker_df <- gene_tbl |>
  dplyr::transmute(
    Marker_ID = as.character(.data$Marker_ID),
    Gene = factor(.data$Gene, levels = gene_order),
    Description = .data$Description
  )

marker_counts <- marker_df |>
  dplyr::distinct(
    .data$Gene,
    .data$Marker_ID
  ) |>
  dplyr::count(
    .data$Gene,
    name = "n_markers"
  )

needed_kegg <- marker_df |>
  dplyr::filter(
    stringr::str_starts(.data$Marker_ID, "K")
  ) |>
  dplyr::pull(.data$Marker_ID) |>
  unique()

needed_pfam <- marker_df |>
  dplyr::filter(
    stringr::str_starts(.data$Marker_ID, "PF")
  ) |>
  dplyr::pull(.data$Marker_ID) |>
  unique()

needed_metab <- marker_df |>
  dplyr::filter(
    !stringr::str_starts(.data$Marker_ID, "K"),
    !stringr::str_starts(.data$Marker_ID, "PF")
  ) |>
  dplyr::pull(.data$Marker_ID) |>
  unique()

############################################################


kegg <- readr::read_tsv(
  kegg_file,
  show_col_types = FALSE
) %>%
  dplyr::rename(
    Marker_ID = .data$KEGG_ko
  ) %>%
  dplyr::filter(
    .data$Marker_ID %in% needed_kegg
  ) %>%
  dplyr::select(
    .data$MAG_name,
    .data$Marker_ID,
    .data$MetaT_sample,
    .data$MT_coverage_per_cell
  )

pfam <- readr::read_tsv(
  pfam_file,
  show_col_types = FALSE
) %>%
  dplyr::rename(
    Marker_ID = .data$PFAM_accession
  ) %>%
  dplyr::filter(
    .data$Marker_ID %in% needed_pfam
  ) %>%
  dplyr::select(
    .data$MAG_name,
    .data$Marker_ID,
    .data$MetaT_sample,
    .data$MT_coverage_per_cell
  )

metab <- readr::read_tsv(
  metab_file,
  show_col_types = FALSE
) %>%
  dplyr::rename(
    Marker_ID = .data$METABOLIC_hmm
  ) %>%
  dplyr::mutate(
    Marker_ID = as.character(.data$Marker_ID),
    Marker_ID = stringr::str_remove(
      .data$Marker_ID,
      "\\.hmm$"
    )
  ) %>%
  dplyr::filter(
    .data$Marker_ID %in% needed_metab
  ) %>%
  dplyr::select(
    .data$MAG_name,
    .data$Marker_ID,
    .data$MetaT_sample,
    .data$MT_coverage_per_cell
  )

df <- dplyr::bind_rows(
  KEGG = kegg,
  PFAM = pfam,
  METABOLIC = metab,
  .id = "Source"
) %>%
  dplyr::mutate(
    MT_coverage_per_cell = suppressWarnings(
      as.numeric(.data$MT_coverage_per_cell)
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$MAG_name),
    !is.na(.data$MetaT_sample),
    !is.na(.data$Marker_ID)
  )



df <- df %>%
  dplyr::left_join(
    tax,
    by = "MAG_name"
  ) %>%
  dplyr::filter(
    !is.na(.data$Phylum)
  ) %>%
  dplyr::mutate(
    Peatland = extract_peatland(
      .data$MetaT_sample
    ),
    Depth_raw = clean_depth(
      dplyr::coalesce(
        stringr::str_extract(
          .data$MetaT_sample,
          "\\d+-\\d+cm|\\d+cm"
        ),
        stringr::str_extract(
          .data$MAG_name,
          "\\d+-\\d+cm|\\d+cm"
        )
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$Peatland)
  ) %>%
  dplyr::left_join(
    depth_map,
    by = "Depth_raw"
  ) %>%
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
    marker_df,
    by = "Marker_ID",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(
    !is.na(.data$Gene),
    !is.na(.data$Depth_layer)
  ) %>%
  dplyr::mutate(
    Gene = factor(
      .data$Gene,
      levels = gene_order
    ),
    Depth_group = dplyr::if_else(
      as.character(.data$Depth_layer) %in% shallow_levels,
      "Shallow",
      "Deep"
    ),
    Depth_group = factor(
      .data$Depth_group,
      levels = c("Shallow", "Deep")
    ),
    Depth_layer_plot = factor(
      as.character(.data$Depth_layer),
      levels = y_levels
    )
  )


avail_PD <- df |>
  dplyr::distinct(
    .data$Peatland,
    .data$Depth_group,
    .data$Depth_layer_plot
  ) |>
  dplyr::mutate(
    Peatland = factor(
      .data$Peatland,
      levels = peatland_levels
    )
  )

############################################################
### 6) GENE EXPRESSION PER MAG PER SAMPLE
###    SUM(MT_coverage_per_cell) / number of markers per gene
############################################################

mag_gene <- df |>
  dplyr::group_by(
    .data$MAG_name,
    .data$MetaT_sample,
    .data$Peatland,
    .data$Depth_group,
    .data$Depth_layer_plot,
    .data$Phylum,
    .data$Gene,
    .data$Description
  ) |>
  dplyr::summarise(
    expr_sum = safe_sum(.data$MT_coverage_per_cell),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    marker_counts,
    by = "Gene"
  ) |>
  dplyr::mutate(
    expr = .data$expr_sum / .data$n_markers
  ) |>
  dplyr::select(
    -dplyr::all_of(c("expr_sum", "n_markers"))
  )

############################################################
### 7) PHYLUM EXPRESSION
############################################################

phy_gene_sample <- mag_gene |>
  dplyr::group_by(
    .data$Peatland,
    .data$Depth_group,
    .data$Depth_layer_plot,
    .data$MetaT_sample,
    .data$Gene,
    .data$Description,
    .data$Phylum
  ) |>
  dplyr::summarise(
    expr = safe_sum(.data$expr),
    .groups = "drop"
  )

phy_gene <- phy_gene_sample |>
  dplyr::group_by(
    .data$Peatland,
    .data$Depth_group,
    .data$Depth_layer_plot,
    .data$Gene,
    .data$Description,
    .data$Phylum
  ) |>
  dplyr::summarise(
    expr = safe_sum(.data$expr),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    expr = dplyr::if_else(
      is.finite(.data$expr),
      .data$expr,
      NA_real_
    )
  )

############################################################
### 8) DROP GENES WITHOUT DETECTABLE EXPRESSION
############################################################

genes_keep <- phy_gene |>
  dplyr::group_by(
    .data$Gene
  ) |>
  dplyr::summarise(
    any_positive = any(
      !is.na(.data$expr) & .data$expr > 0
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(
    .data$any_positive
  ) |>
  dplyr::pull(
    .data$Gene
  ) |>
  as.character()

gene_order_keep <- gene_order[
  gene_order %in% genes_keep
]

phy_gene_keep <- phy_gene |>
  dplyr::filter(
    as.character(.data$Gene) %in% gene_order_keep
  )

############################################################
### 9) PHYLUM RANKING AND TRUE WINNER PER CELL
###
### Winner is selected across ALL phyla first.
### Only afterwards is taxonomy simplified to the selected eight
### phyla plus "Other phyla".
############################################################

phylum_ranking <- phy_gene_keep |>
  dplyr::filter(
    !is.na(.data$expr),
    .data$expr > 0
  ) |>
  dplyr::group_by(
    .data$Phylum
  ) |>
  dplyr::summarise(
    total_MT_coverage_per_cell = sum(
      .data$expr,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$total_MT_coverage_per_cell)
  ) |>
  dplyr::mutate(
    rank = dplyr::row_number(),
    displayed_explicitly = .data$Phylum %in% selected_phyla
  ) |>
  dplyr::select(
    .data$rank,
    .data$Phylum,
    .data$total_MT_coverage_per_cell,
    .data$displayed_explicitly
  )

readr::write_tsv(
  phylum_ranking,
  out_ranking
)


winner_all <- phy_gene_keep |>
  dplyr::filter(
    !is.na(.data$expr),
    .data$expr > 0
  ) |>
  dplyr::group_by(
    .data$Peatland,
    .data$Depth_group,
    .data$Depth_layer_plot,
    .data$Gene,
    .data$Description
  ) |>
  dplyr::slice_max(
    order_by = .data$expr,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::rename(
    Phylum_winner_raw = .data$Phylum
  ) |>
  dplyr::select(
    .data$Peatland,
    .data$Depth_group,
    .data$Depth_layer_plot,
    .data$Gene,
    .data$Description,
    .data$expr,
    .data$Phylum_winner_raw
  )

############################################################
### 10) COMPLETE GRID AND EXPRESSION BINS
############################################################

gene_layout <- gene_tbl |>
  dplyr::filter(
    .data$Gene %in% gene_order_keep
  ) |>
  dplyr::distinct(
    .data$Gene,
    .data$Description
  ) |>
  dplyr::mutate(
    Gene = factor(
      .data$Gene,
      levels = gene_order_keep
    )
  ) |>
  dplyr::arrange(
    .data$Gene
  )

grid <- avail_PD |>
  tidyr::crossing(
    Gene = factor(
      gene_order_keep,
      levels = gene_order_keep
    )
  ) |>
  dplyr::left_join(
    gene_layout,
    by = "Gene"
  )

lab_low <- paste0("≤ ", fmt(thr_low))
lab_medium <- paste0(fmt(thr_low), "–", fmt(thr_high))
lab_high <- paste0("> ", fmt(thr_high))

plot_df <- grid |>
  dplyr::left_join(
    winner_all,
    by = c(
      "Peatland",
      "Depth_group",
      "Depth_layer_plot",
      "Gene",
      "Description"
    )
  ) |>
  dplyr::mutate(
    Gene = factor(
      .data$Gene,
      levels = gene_order_keep
    ),
    x_pos = as.numeric(.data$Gene),
    Level = dplyr::case_when(
      is.na(.data$expr) ~ NA_character_,
      .data$expr <= thr_low ~ lab_low,
      .data$expr <= thr_high ~ lab_medium,
      TRUE ~ lab_high
    ),
    Level = factor(
      .data$Level,
      levels = c(
        lab_low,
        lab_medium,
        lab_high
      )
    ),
    Phylum_plot = dplyr::case_when(
      is.na(.data$Phylum_winner_raw) ~ NA_character_,
      .data$Phylum_winner_raw %in% selected_phyla ~ .data$Phylum_winner_raw,
      TRUE ~ "Other phyla"
    ),
    Phylum_plot = factor(
      .data$Phylum_plot,
      levels = phylum_levels
    )
  )


bubble_points <- plot_df |>
  dplyr::filter(
    !is.na(.data$Level),
    !is.na(.data$Phylum_plot)
  )

cat("\nBubble rows retained for plotting: ", nrow(bubble_points), "\n", sep = "")

size_values <- stats::setNames(
  c(1.9, 3.6, 5.4),
  c(lab_low, lab_medium, lab_high)
)

readr::write_tsv(
  plot_df |>
    dplyr::mutate(
      Gene = as.character(.data$Gene),
      Peatland = as.character(.data$Peatland),
      Depth_group = as.character(.data$Depth_group),
      Depth_layer_plot = as.character(.data$Depth_layer_plot),
      Phylum_plot = as.character(.data$Phylum_plot),
      Level = as.character(.data$Level)
    ),
  out_values
)

############################################################
### 11) LABELS
############################################################

description_map <- gene_layout |>
  dplyr::mutate(
    Description = dplyr::if_else(
      is.na(.data$Description),
      "",
      .data$Description
    )
  ) |>
  dplyr::arrange(
    .data$Gene
  )

x_breaks <- seq_along(gene_order_keep)
top_labels <- description_map$Description

############################################################
### 12) BUBBLE PLOT
############################################################

p_bubble <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = .data$x_pos,
    y = .data$Depth_layer_plot
  )
) +
  ggplot2::geom_tile(
    fill = "white",
    colour = "grey82",
    linewidth = 0.25,
    width = 0.92,
    height = 0.92
  ) +
  ggplot2::geom_point(
    data = bubble_points,
    ggplot2::aes(
      fill = .data$Phylum_plot,
      size = .data$Level
    ),
    shape = 21,
    colour = "grey15",
    stroke = 0.30
  ) +
  ggplot2::facet_grid(
    Peatland + Depth_group ~ .,
    switch = "y",
    scales = "free_y",
    space = "free_y",
    labeller = ggplot2::labeller(
      Depth_group = function(x) {
        rep("", length(x))
      }
    )
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = gene_order_keep,
    expand = ggplot2::expansion(mult = c(0, 0)),
    sec.axis = ggplot2::dup_axis(
      breaks = x_breaks,
      labels = top_labels,
      name = NULL
    )
  ) +
  ggplot2::scale_fill_manual(
    values = phylum_colors,
    breaks = phylum_levels,
    limits = phylum_levels,
    labels = function(x) {
      stringr::str_remove(x, "^p__")
    },
    drop = FALSE,
    name = "Dominant phylum"
  ) +
  ggplot2::scale_size_manual(
    values = size_values,
    breaks = c(lab_low, lab_medium, lab_high),
    limits = c(lab_low, lab_medium, lab_high),
    drop = FALSE,
    name = "MT coverage per cell"
  ) +
  ggplot2::scale_y_discrete(
    drop = TRUE,
    labels = function(x) {
      stringr::str_replace_all(x, "\\s*cm", "")
    },
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      order = 1,
      ncol = 5,
      byrow = TRUE,
      override.aes = list(
        shape = 21,
        size = 4,
        colour = "grey15"
      )
    ),
    size = ggplot2::guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "grey70",
        colour = "grey15"
      )
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Depth layer"
  ) +
  ggplot2::theme_bw(
    base_size = 12,
    base_family = plot_font
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.title.x = ggplot2::element_blank(),
    axis.text.x.bottom = ggplot2::element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.x.top = ggplot2::element_text(
      angle = 60,
      hjust = 0,
      vjust = 0,
      size = 8
    ),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold"
    ),
    panel.spacing.y = grid::unit(
      gap_pt,
      "pt"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.title.align = 0.5,
    plot.margin = ggplot2::margin(
      t = 10,
      r = 15,
      b = 8,
      l = 10,
      unit = "pt"
    )
  )

print(p_bubble)

############################################################
### 13) SAVE OUTPUTS
############################################################

ggplot2::ggsave(
  filename = out_svg,
  plot = p_bubble,
  device = svglite::svglite,
  width = 20,
  height = 8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

if (capabilities("cairo")) {
  ggplot2::ggsave(
    filename = out_pdf,
    plot = p_bubble,
    device = grDevices::cairo_pdf,
    width = 16,
    height = 9,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
} else {
  warning(
    "PDF was not saved because this R installation lacks Cairo support."
  )
}

cat("\nSaved SVG:\n", out_svg, "\n", sep = "")

if (capabilities("cairo")) {
  cat("\nSaved PDF:\n", out_pdf, "\n", sep = "")
}

cat("\nSaved plotted values:\n", out_values, "\n", sep = "")
cat("\nSaved phylum ranking:\n", out_ranking, "\n", sep = "")
cat("\nDONE. Metric used throughout: MT_coverage_per_cell.\n")
