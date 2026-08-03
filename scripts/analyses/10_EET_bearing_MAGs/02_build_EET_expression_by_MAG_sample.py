#!/usr/bin/env python3

# DESCRIPTION
# Maps curated EET proteins to MAGs and builds EET gene-, component- and
# machinery-level expression tables for each MAG and metatranscriptomic
# sample.
#
# INPUT
# 1. Wide EET-gene expression TSV.
# 2. Analysis directory containing outputs from script 09.
# 3. Directory containing dereplicated MAG FASTA files.
#
# OUTPUT
# Protein-to-MAG mapping tables, long EET expression tables, machinery-level
# summaries and QC files.
#
# USAGE
# python3 10_build_EET_expression_by_MAG_sample.py \
#   EET_gene_expression_wide.tsv \
#   analysis_output_directory \
#   dereplicated_MAG_FASTA_directory

from collections import defaultdict
from pathlib import Path
import argparse
import re
import sys

import pandas as pd


OUT_DIR = None
EET_EXPRESSION_FILE = None
EET_MEMBERSHIP_FILE = None
ALTERNATIVE_MODULE_FILE = None
DEREPLICATED_GENOMES_DIR = None

FILL_MISSING_EXPRESSION_WITH_ZERO = True


###############################################################################
# NAME NORMALISATION
###############################################################################

def remove_fasta_suffix(name):
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", str(name), flags=re.I)


def canonical_mag_name(value):
    return remove_fasta_suffix(Path(str(value)).name.strip())


def canonical_sample_name(value):
    text = str(value).strip()
    text = re.sub(r"^MT_coverage_per_cell_", "", text)
    text = re.sub(r"_2024-06_iMT$", "", text)
    return text


def protein_to_contig_id(protein_id):
    """
    Prodigal IDs are expected to be:
        <contig_ID>_<gene_number>
    Example:
        NODE_1703_length_55991_cov_1.550211_34
        -> NODE_1703_length_55991_cov_1.550211
    """
    protein_id = str(protein_id).strip()
    match = re.match(r"^(.+)_([0-9]+)$", protein_id)

    if not match:
        return None

    return match.group(1)


###############################################################################
# INPUT TABLES
###############################################################################

def read_membership():
    if not EET_MEMBERSHIP_FILE.is_file():
        raise FileNotFoundError(EET_MEMBERSHIP_FILE)

    membership = pd.read_csv(EET_MEMBERSHIP_FILE, sep="\t", dtype=str)

    required = {"MAG_name", "EET_groups"}
    missing = required - set(membership.columns)
    if missing:
        raise RuntimeError(
            f"Missing columns in {EET_MEMBERSHIP_FILE}: {sorted(missing)}"
        )

    membership = membership[["MAG_name", "EET_groups"]].drop_duplicates().copy()
    membership = membership.rename(columns={"MAG_name": "source_MAG_name"})
    membership["source_MAG_name"] = membership["source_MAG_name"].map(str).str.strip()

    # Use MAG names exactly as supplied in the membership table.
    membership["MAG_name"] = membership["source_MAG_name"].map(canonical_mag_name)
    membership["canonical_MAG_name"] = membership["MAG_name"]

    duplicated_canonical = membership[
        membership["canonical_MAG_name"].duplicated(keep=False)
    ].sort_values("canonical_MAG_name")

    if not duplicated_canonical.empty:
        duplicated_canonical.to_csv(
            OUT_DIR / "QC_duplicated_canonical_MAG_names.tsv",
            sep="\t",
            index=False,
        )
        raise RuntimeError(
            "Duplicated MAG names were found in the membership table. "
            "See QC_duplicated_canonical_MAG_names.tsv"
        )


    return membership


def read_reference_mt_samples():
    if not ALTERNATIVE_MODULE_FILE.is_file():
        raise FileNotFoundError(ALTERNATIVE_MODULE_FILE)

    samples = pd.read_csv(
        ALTERNATIVE_MODULE_FILE,
        sep="\t",
        usecols=["MetaT_sample"],
        dtype=str,
    )["MetaT_sample"].dropna().drop_duplicates()

    reference = pd.DataFrame({"source_MetaT_sample": sorted(samples)})
    reference["canonical_sample"] = reference["source_MetaT_sample"].map(
        canonical_sample_name
    )

    # Reconstruct the complete MetaT sample name without changing it.
    reference["MetaT_sample"] = (
        reference["canonical_sample"] + "_2024-06_iMT"
    )

    duplicated = reference[
        reference["canonical_sample"].duplicated(keep=False)
    ].sort_values("canonical_sample")

    if not duplicated.empty:
        duplicated.to_csv(
            OUT_DIR / "QC_duplicated_canonical_MetaT_samples.tsv",
            sep="\t",
            index=False,
        )
        raise RuntimeError(
            "Multiple MetaT samples collapse to the same canonical sample. "
            "See QC_duplicated_canonical_MetaT_samples.tsv"
        )


    return reference


###############################################################################
# FASTA-BASED PROTEIN -> MAG MAPPING
###############################################################################

def find_eet_mag_fasta_files(membership):
    """
    Resolve each MAG by exact filename construction.

    MAG names are used exactly as supplied in the membership table.
    No aliases, depth corrections or approximate matching are applied.
    """
    rows = []
    missing = []

    for record in membership.itertuples(index=False):
        fasta_path = DEREPLICATED_GENOMES_DIR / f"{record.MAG_name}.fa"

        if not fasta_path.is_file():
            missing.append(
                {
                    "MAG_name": record.MAG_name,
                    "path_tested": str(fasta_path),
                }
            )
            continue

        rows.append(
            {
                "source_MAG_name": record.source_MAG_name,
                "MAG_name": record.MAG_name,
                "EET_groups": record.EET_groups,
                "canonical_MAG_name": record.canonical_MAG_name,
                "matched_FASTA_name": record.MAG_name,
                "FASTA_file": str(fasta_path),
            }
        )

    missing_df = pd.DataFrame(
        missing,
        columns=[
            "MAG_name",
            "path_tested",
        ],
    )
    missing_df.to_csv(
        OUT_DIR / "QC_EET_MAGs_without_FASTA.tsv",
        sep="	",
        index=False,
    )

    if missing:
        raise RuntimeError(
            f"Could not resolve FASTA files for {len(missing)} EET MAGs. "
            "See QC_EET_MAGs_without_FASTA.tsv"
        )

    fasta_table = pd.DataFrame(rows)

    if fasta_table["FASTA_file"].duplicated().any():
        duplicated = fasta_table[
            fasta_table["FASTA_file"].duplicated(keep=False)
        ].sort_values("FASTA_file")
        duplicated.to_csv(
            OUT_DIR / "QC_duplicate_FASTA_assignments.tsv",
            sep="	",
            index=False,
        )
        raise RuntimeError(
            "The same FASTA was assigned to more than one EET MAG. "
            "See QC_duplicate_FASTA_assignments.tsv"
        )

    fasta_table.to_csv(
        OUT_DIR / "03_EET_MAG_FASTA_files.tsv",
        sep="	",
        index=False,
    )

    return fasta_table


def build_contig_to_mag_index(fasta_table):
    contig_to_mags = defaultdict(set)
    n_contigs_by_mag = {}

    for row in fasta_table.itertuples(index=False):
        fasta_path = Path(row.FASTA_file)
        n_contigs = 0

        with fasta_path.open() as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue

                contig_id = line[1:].strip().split()[0]
                contig_to_mags[contig_id].add(row.MAG_name)
                n_contigs += 1

        n_contigs_by_mag[row.MAG_name] = n_contigs

    contig_counts = pd.DataFrame(
        sorted(n_contigs_by_mag.items()),
        columns=["MAG_name", "n_contigs"],
    )
    contig_counts.to_csv(
        OUT_DIR / "QC_EET_MAG_contig_counts.tsv",
        sep="\t",
        index=False,
    )

    return contig_to_mags


def parse_group_set(value):
    return {
        item.strip()
        for item in str(value).split(";")
        if item.strip()
    }


def map_proteins_to_mags(eet_wide, membership, contig_to_mags):
    membership_groups = {
        row.MAG_name: parse_group_set(row.EET_groups)
        for row in membership.itertuples(index=False)
    }

    unique_proteins = (
        eet_wide[["Protein_ID", "EET_type", "EET_group"]]
        .drop_duplicates()
        .copy()
    )
    unique_proteins["contig_ID"] = unique_proteins["Protein_ID"].map(
        protein_to_contig_id
    )

    mapping_rows = []
    unresolved_rows = []

    for row in unique_proteins.itertuples(index=False):
        if row.contig_ID is None:
            unresolved_rows.append(
                {
                    "Protein_ID": row.Protein_ID,
                    "EET_type": row.EET_type,
                    "EET_group": row.EET_group,
                    "contig_ID": "",
                    "candidate_MAGs_before_group_filter": "",
                    "candidate_MAGs_after_group_filter": "",
                    "problem": "Protein_ID does not end with _<integer>",
                }
            )
            continue

        candidates_before = sorted(contig_to_mags.get(row.contig_ID, set()))

        candidates_after = [
            mag
            for mag in candidates_before
            if row.EET_group in membership_groups.get(mag, set())
        ]

        if len(candidates_after) == 1:
            mapping_rows.append(
                {
                    "Protein_ID": row.Protein_ID,
                    "EET_type": row.EET_type,
                    "EET_group": row.EET_group,
                    "contig_ID": row.contig_ID,
                    "MAG_name": candidates_after[0],
                    "n_candidate_MAGs_before_group_filter": len(candidates_before),
                }
            )
        else:
            if len(candidates_before) == 0:
                problem = "Contig not found in the EET MAG FASTAs"
            elif len(candidates_after) == 0:
                problem = (
                    "Contig found, but no candidate MAG belongs to this EET group"
                )
            else:
                problem = "Contig remains ambiguous after EET-group filtering"

            unresolved_rows.append(
                {
                    "Protein_ID": row.Protein_ID,
                    "EET_type": row.EET_type,
                    "EET_group": row.EET_group,
                    "contig_ID": row.contig_ID,
                    "candidate_MAGs_before_group_filter": ";".join(
                        candidates_before
                    ),
                    "candidate_MAGs_after_group_filter": ";".join(
                        candidates_after
                    ),
                    "problem": problem,
                }
            )

    mapping = pd.DataFrame(mapping_rows)
    unresolved = pd.DataFrame(
        unresolved_rows,
        columns=[
            "Protein_ID",
            "EET_type",
            "EET_group",
            "contig_ID",
            "candidate_MAGs_before_group_filter",
            "candidate_MAGs_after_group_filter",
            "problem",
        ],
    )

    mapping.to_csv(
        OUT_DIR / "04_EET_Protein_ID_to_MAG.tsv",
        sep="\t",
        index=False,
    )
    unresolved.to_csv(
        OUT_DIR / "QC_unresolved_EET_Protein_ID_to_MAG.tsv",
        sep="\t",
        index=False,
    )

    if not unresolved.empty:
        raise RuntimeError(
            f"{len(unresolved)} EET proteins could not be assigned uniquely "
            "to a MAG. See QC_unresolved_EET_Protein_ID_to_MAG.tsv"
        )

    # Protein_ID values such as NODE_<...>_<gene_number> are generated
    # independently within different assemblies and therefore are not
    # guaranteed to be globally unique across all MAGs. The valid unique key
    # for this table is Protein_ID + EET_type + EET_group.
    duplicate_protein_ids = mapping[
        mapping["Protein_ID"].duplicated(keep=False)
    ].sort_values(
        ["Protein_ID", "EET_group", "EET_type", "MAG_name"]
    )

    duplicate_protein_ids.to_csv(
        OUT_DIR / "QC_reused_EET_Protein_IDs_across_contexts.tsv",
        sep="\t",
        index=False,
    )

    composite_key = ["Protein_ID", "EET_type", "EET_group"]
    duplicated_composite = mapping[
        mapping.duplicated(subset=composite_key, keep=False)
    ].sort_values(composite_key + ["MAG_name"])

    duplicated_composite.to_csv(
        OUT_DIR / "QC_duplicated_EET_composite_mapping_key.tsv",
        sep="\t",
        index=False,
    )

    if not duplicated_composite.empty:
        raise RuntimeError(
            "At least one Protein_ID + EET_type + EET_group combination "
            "has multiple mappings. See "
            "QC_duplicated_EET_composite_mapping_key.tsv"
        )

    mapped_mags = set(mapping["MAG_name"])
    expected_mags = set(membership["MAG_name"])
    missing_mags = sorted(expected_mags - mapped_mags)

    Path(OUT_DIR / "QC_EET_MAGs_without_mapped_EET_proteins.txt").write_text(
        "\n".join(missing_mags) + ("\n" if missing_mags else "")
    )

    # The EET expression table contains a curated subset of EET genes rather
    # than genes from all 120 EET-bearing MAGs. Therefore, downstream
    # expression analyses must use only the MAGs represented in this curated
    # EET expression subset. MAGs without mapped EET proteins are retained in
    # the QC file and are not treated as biological zeros.
    if len(mapped_mags) == 0:
        raise RuntimeError(
            "No EET proteins could be mapped to EET MAGs."
        )

    included_mags = (
        membership[membership["MAG_name"].isin(mapped_mags)]
        [["MAG_name", "EET_groups"]]
        .drop_duplicates()
        .sort_values(["EET_groups", "MAG_name"])
    )
    included_mags.to_csv(
        OUT_DIR / "04b_EET_MAGs_in_expression_subset.tsv",
        sep="\t",
        index=False,
    )

    excluded_mags = (
        membership[~membership["MAG_name"].isin(mapped_mags)]
        [["MAG_name", "EET_groups"]]
        .drop_duplicates()
        .sort_values(["EET_groups", "MAG_name"])
    )
    excluded_mags.to_csv(
        OUT_DIR / "04c_EET_MAGs_not_in_expression_subset.tsv",
        sep="\t",
        index=False,
    )

    return mapping


###############################################################################
# EET EXPRESSION RESHAPING AND SUMMARISATION
###############################################################################

def read_and_validate_eet_wide():
    if not EET_EXPRESSION_FILE.is_file():
        raise FileNotFoundError(EET_EXPRESSION_FILE)

    eet = pd.read_csv(
        EET_EXPRESSION_FILE,
        sep="\t",
        dtype={"Protein_ID": str, "EET_type": str},
        low_memory=False,
    )

    required = {"Protein_ID", "EET_type"}
    missing = required - set(eet.columns)
    if missing:
        raise RuntimeError(
            f"Missing columns in {EET_EXPRESSION_FILE}: {sorted(missing)}"
        )

    sample_columns = [
        column
        for column in eet.columns
        if column.startswith("MT_coverage_per_cell_")
    ]


    invalid_types = eet[
        ~eet["EET_type"].str.match(r"^EET[1-4]_.+", na=False)
    ][["Protein_ID", "EET_type"]].drop_duplicates()

    invalid_types.to_csv(
        OUT_DIR / "QC_invalid_EET_type_values.tsv",
        sep="\t",
        index=False,
    )

    if not invalid_types.empty:
        raise RuntimeError(
            "Invalid EET_type values found. See QC_invalid_EET_type_values.tsv"
        )

    eet["EET_group"] = eet["EET_type"].str.extract(
        r"^(EET[1-4])_", expand=False
    )
    eet["EET_component"] = eet["EET_type"].str.replace(
        r"^EET[1-4]_", "", regex=True
    )

    duplicated = eet[
        eet.duplicated(subset=["Protein_ID", "EET_type"], keep=False)
    ].sort_values(["Protein_ID", "EET_type"])

    duplicated.to_csv(
        OUT_DIR / "QC_duplicated_Protein_ID_EET_type_rows.tsv",
        sep="\t",
        index=False,
    )

    if not duplicated.empty:
        raise RuntimeError(
            "Duplicated Protein_ID + EET_type rows found. "
            "See QC_duplicated_Protein_ID_EET_type_rows.tsv"
        )

    return eet, sample_columns


def map_sample_columns(sample_columns, reference_samples):
    reference_map = {
        row.canonical_sample: (
            row.source_MetaT_sample,
            row.MetaT_sample,
        )
        for row in reference_samples.itertuples(index=False)
    }

    rows = []

    for source_column in sample_columns:
        source_label = re.sub(
            r"^MT_coverage_per_cell_", "", source_column
        )
        canonical = canonical_sample_name(source_label)
        matched = reference_map.get(canonical)

        source_reference_sample = matched[0] if matched else None
        corrected_metat_sample = matched[1] if matched else None

        rows.append(
            {
                "source_expression_column": source_column,
                "source_sample_label": source_label,
                "canonical_sample": canonical,
                "source_reference_MetaT_sample": source_reference_sample,
                "MetaT_sample": corrected_metat_sample,
            }
        )

    mapping = pd.DataFrame(rows)
    mapping.to_csv(
        OUT_DIR / "05_EET_sample_name_mapping.tsv",
        sep="\t",
        index=False,
    )

    unmatched = mapping[mapping["MetaT_sample"].isna()]

    if not unmatched.empty:
        unmatched.to_csv(
            OUT_DIR / "QC_unmatched_EET_sample_columns.tsv",
            sep="\t",
            index=False,
        )
        raise RuntimeError(
            f"{len(unmatched)} EET expression columns could not be matched to "
            "the 20 MetaT samples. See QC_unmatched_EET_sample_columns.tsv"
        )

    if mapping["MetaT_sample"].duplicated().any():
        duplicated = mapping[
            mapping["MetaT_sample"].duplicated(keep=False)
        ].sort_values("MetaT_sample")
        duplicated.to_csv(
            OUT_DIR / "QC_duplicated_EET_sample_mapping.tsv",
            sep="\t",
            index=False,
        )
        raise RuntimeError(
            "More than one EET expression column maps to the same MetaT sample. "
            "See QC_duplicated_EET_sample_mapping.tsv"
        )


    return mapping


def reshape_eet_expression(eet_wide, sample_columns, protein_mapping, sample_mapping):
    long = eet_wide.melt(
        id_vars=[
            "Protein_ID",
            "EET_type",
            "EET_group",
            "EET_component",
        ],
        value_vars=sample_columns,
        var_name="source_expression_column",
        value_name="MT_coverage_per_cell",
    )

    long = long.merge(
        sample_mapping[
            [
                "source_expression_column",
                "source_sample_label",
                "source_reference_MetaT_sample",
                "MetaT_sample",
            ]
        ],
        on="source_expression_column",
        how="left",
        validate="many_to_one",
    )

    long = long.merge(
        protein_mapping[
            [
                "Protein_ID",
                "EET_type",
                "EET_group",
                "contig_ID",
                "MAG_name",
            ]
        ],
        on=["Protein_ID", "EET_type", "EET_group"],
        how="left",
        validate="many_to_one",
    )

    missing_before_fill = long["MT_coverage_per_cell"].isna().sum()
    long["MT_coverage_per_cell"] = pd.to_numeric(
        long["MT_coverage_per_cell"],
        errors="coerce",
    )
    missing_after_numeric_conversion = long["MT_coverage_per_cell"].isna().sum()

    if missing_after_numeric_conversion:
        if FILL_MISSING_EXPRESSION_WITH_ZERO:
            long["MT_coverage_per_cell"] = (
                long["MT_coverage_per_cell"].fillna(0.0)
            )
        else:
            long[long["MT_coverage_per_cell"].isna()].to_csv(
                OUT_DIR / "QC_missing_EET_expression_values.tsv",
                sep="\t",
                index=False,
            )
            raise RuntimeError(
                "Missing or non-numeric EET expression values found. "
                "See QC_missing_EET_expression_values.tsv"
            )

    long["EET_gene_is_expressed"] = long["MT_coverage_per_cell"] > 0
    long["EET_feature_ID"] = (
        long["MAG_name"].astype(str)
        + "|"
        + long["Protein_ID"].astype(str)
        + "|"
        + long["EET_type"].astype(str)
    )

    long = long[
        [
            "EET_feature_ID",
            "MAG_name",
            "Protein_ID",
            "contig_ID",
            "EET_group",
            "EET_component",
            "EET_type",
            "source_sample_label",
            "source_reference_MetaT_sample",
            "MetaT_sample",
            "MT_coverage_per_cell",
            "EET_gene_is_expressed",
        ]
    ].sort_values(
        [
            "MAG_name",
            "EET_group",
            "Protein_ID",
            "MetaT_sample",
        ]
    )

    long.to_csv(
        OUT_DIR / "06_EET_gene_expression_by_MAG_sample.tsv",
        sep="\t",
        index=False,
    )

    return long, missing_before_fill, missing_after_numeric_conversion


def summarise_eet_expression(gene_expression):
    component_keys = [
        "MAG_name",
        "EET_group",
        "EET_component",
        "EET_type",
        "MetaT_sample",
    ]

    component = (
        gene_expression
        .groupby(component_keys, as_index=False, dropna=False)
        .agg(
            n_EET_proteins=("EET_feature_ID", "nunique"),
            n_EET_proteins_expressed=(
                "EET_gene_is_expressed",
                "sum",
            ),
            sum_MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                "sum",
            ),
            mean_MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                "mean",
            ),
            median_MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                "median",
            ),
        )
    )

    component["EET_component_is_expressed"] = (
        component["sum_MT_coverage_per_cell"] > 0
    )

    component.to_csv(
        OUT_DIR / "07_EET_component_expression_by_MAG_sample.tsv",
        sep="\t",
        index=False,
    )

    machinery_gene = (
        gene_expression
        .groupby(
            ["MAG_name", "EET_group", "MetaT_sample"],
            as_index=False,
            dropna=False,
        )
        .agg(
            n_EET_proteins=("EET_feature_ID", "nunique"),
            n_EET_proteins_expressed=("EET_gene_is_expressed", "sum"),
            sum_gene_MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                "sum",
            ),
            mean_gene_MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                "mean",
            ),
            median_gene_MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                "median",
            ),
        )
    )

    machinery_component = (
        component
        .groupby(
            ["MAG_name", "EET_group", "MetaT_sample"],
            as_index=False,
            dropna=False,
        )
        .agg(
            n_EET_components=("EET_component", "nunique"),
            n_EET_components_expressed=(
                "EET_component_is_expressed",
                "sum",
            ),
            mean_component_sum_MT_coverage_per_cell=(
                "sum_MT_coverage_per_cell",
                "mean",
            ),
            median_component_sum_MT_coverage_per_cell=(
                "sum_MT_coverage_per_cell",
                "median",
            ),
            mean_component_mean_MT_coverage_per_cell=(
                "mean_MT_coverage_per_cell",
                "mean",
            ),
        )
    )

    machinery = machinery_gene.merge(
        machinery_component,
        on=["MAG_name", "EET_group", "MetaT_sample"],
        how="inner",
        validate="one_to_one",
    )

    machinery["EET_machinery_is_expressed"] = (
        machinery["sum_gene_MT_coverage_per_cell"] > 0
    )
    machinery["all_EET_components_expressed"] = (
        machinery["n_EET_components_expressed"]
        == machinery["n_EET_components"]
    )

    machinery.to_csv(
        OUT_DIR / "08_EET_machinery_expression_by_MAG_sample.tsv",
        sep="\t",
        index=False,
    )

    return component, machinery



def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Map curated EET proteins to MAGs and summarize EET "
            "expression by MAG and sample."
        )
    )
    parser.add_argument(
        "eet_expression_wide_tsv",
        type=Path,
    )
    parser.add_argument(
        "analysis_directory",
        type=Path,
    )
    parser.add_argument(
        "dereplicated_genomes_directory",
        type=Path,
    )
    return parser.parse_args()


###############################################################################
# MAIN
###############################################################################

def main():
    global OUT_DIR
    global EET_EXPRESSION_FILE
    global EET_MEMBERSHIP_FILE
    global ALTERNATIVE_MODULE_FILE
    global DEREPLICATED_GENOMES_DIR

    args = parse_args()

    OUT_DIR = args.analysis_directory
    EET_EXPRESSION_FILE = args.eet_expression_wide_tsv
    EET_MEMBERSHIP_FILE = OUT_DIR / "00_EET_MAG_membership.tsv"
    ALTERNATIVE_MODULE_FILE = (
        OUT_DIR / "02_selected_modules_expression_by_MAG_sample.tsv"
    )
    DEREPLICATED_GENOMES_DIR = args.dereplicated_genomes_directory

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    membership = read_membership()
    reference_samples = read_reference_mt_samples()

    eet_wide, sample_columns = read_and_validate_eet_wide()
    sample_mapping = map_sample_columns(sample_columns, reference_samples)

    fasta_table = find_eet_mag_fasta_files(membership)
    contig_to_mags = build_contig_to_mag_index(fasta_table)
    protein_mapping = map_proteins_to_mags(
        eet_wide,
        membership,
        contig_to_mags,
    )

    gene_expression, n_missing_before, n_missing_after = reshape_eet_expression(
        eet_wide,
        sample_columns,
        protein_mapping,
        sample_mapping,
    )

    component_expression, machinery_expression = summarise_eet_expression(
        gene_expression
    )

    mapped_mags = gene_expression["MAG_name"].nunique()
    mapped_samples = gene_expression["MetaT_sample"].nunique()
    mapped_proteins = gene_expression["EET_feature_ID"].nunique()

    expected_gene_rows = mapped_proteins * mapped_samples
    if len(gene_expression) != expected_gene_rows:
        raise RuntimeError(
            f"Expected {expected_gene_rows} gene-level rows "
            f"({mapped_proteins} proteins x {mapped_samples} samples), "
            f"found {len(gene_expression)}."
        )

    total_membership_mags = membership["MAG_name"].nunique()
    reference_sample_count = reference_samples["MetaT_sample"].nunique()

    qc_lines = [
        f"Total EET-bearing MAGs in membership files: {total_membership_mags}",
        f"EET MAGs represented in curated EET expression table: {mapped_mags}",
        (
            "EET MAGs excluded because no curated EET gene was present: "
            f"{total_membership_mags - mapped_mags}"
        ),
        f"MetaT samples in reference table: {reference_sample_count}",
        f"MetaT samples mapped: {mapped_samples}",
        f"Unique EET proteins: {mapped_proteins}",
        f"EET gene-level rows: {len(gene_expression)}",
        f"EET component-level rows: {len(component_expression)}",
        f"EET machinery-level rows: {len(machinery_expression)}",
        f"Missing wide-table values before numeric conversion: {n_missing_before}",
        f"Missing/non-numeric values after numeric conversion: {n_missing_after}",
        (
            "Missing/non-numeric expression values filled with zero: "
            f"{FILL_MISSING_EXPRESSION_WITH_ZERO}"
        ),
        (
            "EET MAG-sample-group rows with detectable machinery expression: "
            f"{int(machinery_expression['EET_machinery_is_expressed'].sum())}"
        ),
        (
            "EET MAG-sample-group rows with all machinery components expressed: "
            f"{int(machinery_expression['all_EET_components_expressed'].sum())}"
        ),
    ]

    # The curated EET-gene table is not expected to include all 120
    # EET-bearing MAGs. The final table therefore contains only the mapped
    # subset and this subset must be used in all correlations and models.
    if mapped_mags == 0:
        raise RuntimeError(
            "The final EET expression table contains no mapped MAGs."
        )

    if mapped_samples != reference_sample_count:
        raise RuntimeError(
            "The final EET expression table does not contain the same "
            "number of MetaT samples as the reference table."
        )

    (OUT_DIR / "QC_EET_expression_summary.txt").write_text(
        "\n".join(qc_lines) + "\n"
    )

    print("\n[COMPLETED]")
    for line in qc_lines:
        print(f"  {line}")

    print(f"\nOutputs written to:\n  {OUT_DIR}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"\n[ERROR] {error}", file=sys.stderr)
        sys.exit(1)