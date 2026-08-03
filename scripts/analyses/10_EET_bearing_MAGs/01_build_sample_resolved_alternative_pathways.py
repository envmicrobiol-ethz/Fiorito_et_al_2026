#!/usr/bin/env python3

# DESCRIPTION
# Builds sample-resolved expression tables for selected alternative
# respiration pathways and hydrogenases in EET-bearing MAGs.
#
# INPUT
# 1. Directory containing EET1-EET4 genome-name files.
# 2. MAG-level KEGG expression-profile TSV.
# 3. MAG-level PFAM expression-profile TSV.
# 4. MAG-level METABOLIC HMM expression-profile TSV.
# 5. Output directory.
#
# OUTPUT
# EET MAG membership, marker-level expression, module-level expression
# and QC files.
#
# USAGE
# python3 01_build_sample_resolved_alternative_pathways.py \
#   EET_membership_directory \
#   MAG_KEGG_expression.tsv \
#   MAG_PFAM_expression.tsv \
#   MAG_METABOLIC_HMM_expression.tsv \
#   output_directory

from pathlib import Path
import argparse
import re
import sys

import pandas as pd


DRAM_DIR = None
KEGG_FILE = None
PFAM_FILE = None
HMM_FILE = None
OUT_DIR = None

CHUNK_SIZE = 500_000


###############################################################################
# SELECTED MARKERS
###############################################################################

MARKERS = [
    ("Oxygen metabolism", "coxA", "K02274"),
    ("Oxygen metabolism", "coxB", "K02275"),
    ("Oxygen metabolism", "ccoN", "K00404"),
    ("Oxygen metabolism", "ccoO", "K00405"),
    ("Oxygen metabolism", "ccoP", "K00406"),
    ("Oxygen metabolism", "cyoA", "K02297"),
    ("Oxygen metabolism", "cyoD", "K02300"),
    ("Oxygen metabolism", "cyoB", "K02298"),
    ("Oxygen metabolism", "cyoC", "K02299"),
    ("Oxygen metabolism", "cydA", "K00425"),
    ("Oxygen metabolism", "cydB", "K00426"),

    ("Aromatics degradation", "bcrC", "K04112"),
    ("Aromatics degradation", "bcrB", "K04113"),
    ("Aromatics degradation", "bcrA", "K04114"),
    ("Aromatics degradation", "bcrD", "K04115"),

    ("Nitrogen cycling", "nxrA", "K00370"),
    ("Nitrogen cycling", "nxrB", "K00371"),
    ("Nitrogen cycling", "napA", "K02567"),
    ("Nitrogen cycling", "napB", "K02568"),
    ("Nitrogen cycling", "narG", "K00370"),
    ("Nitrogen cycling", "narH", "K00371"),
    ("Nitrogen cycling", "nrfH", "K15876"),
    ("Nitrogen cycling", "nrfA", "K03385"),
    ("Nitrogen cycling", "nirB", "K00362"),
    ("Nitrogen cycling", "nirK", "K00368"),
    ("Nitrogen cycling", "nirS", "K15864"),
    ("Nitrogen cycling", "norB", "K04561"),
    ("Nitrogen cycling", "norC", "K02305"),
    ("Nitrogen cycling", "nosD", "PF05048"),
    ("Nitrogen cycling", "nosZ", "K00376"),

    ("Sulfur cycling", "dsrD", "PF08679"),
    ("Sulfur cycling", "dsrA", "K11180"),
    ("Sulfur cycling", "dsrB", "K11181"),
    ("Sulfur cycling", "asrA", "K16950"),
    ("Sulfur cycling", "asrB", "K16951"),
    ("Sulfur cycling", "asrC", "K00385"),
    ("Sulfur cycling", "aprA", "K00394"),
    ("Sulfur cycling", "dmsA", "K07306"),
    ("Sulfur cycling", "phsA", "K08352"),

    ("Chlorate reduction", "ClrB", "K17051"),

    ("Arsenate reduction", "arsC (grx)", "K00537"),
    ("Arsenate reduction", "arsC (trx)", "K03741"),

    ("Selenate reduction", "ygfM", "K12529"),
    ("Selenate reduction", "xdhD", "K12528"),
    ("Selenate reduction", "YgfK", "K12527"),

    ("Fumarate reduction", "frdA", "K00244"),
    ("Fumarate reduction", "frdB", "K00245"),
    ("Fumarate reduction", "frdC", "K00246"),
    ("Fumarate reduction", "frdD", "K00247"),

    ("Hydrogenases", "fefe-group-a13", "fefe-group-a13.hmm"),
    ("Hydrogenases", "fefe-group-a2", "fefe-group-a2.hmm"),
    ("Hydrogenases", "fefe-group-a4", "fefe-group-a4.hmm"),
    ("Hydrogenases", "fefe-group-b", "fefe-group-b.hmm"),
    ("Hydrogenases", "fefe-group-c2", "fefe-group-c2.hmm"),
    ("Hydrogenases", "fefe-group-c3", "fefe-group-c3.hmm"),
    ("Hydrogenases", "nife-group-1", "nife-group-1.hmm"),
    ("Hydrogenases", "nife-group-2bc", "nife-group-2bc.hmm"),
    ("Hydrogenases", "nife-group-3abd", "nife-group-3abd.hmm"),
    ("Hydrogenases", "nife-group-3c", "nife-group-3c.hmm"),
    ("Hydrogenases", "nife-group-4a-g", "nife-group-4a-g.hmm"),
]


MODULE_ROLE = {
    "Oxygen metabolism": ("alternative_terminal_respiration", True),
    "Aromatics degradation": ("electron_donor_metabolism", False),
    "Nitrogen cycling": ("potential_alternative_terminal_respiration", True),
    "Sulfur cycling": ("potential_alternative_terminal_respiration", True),
    "Chlorate reduction": ("potential_alternative_terminal_respiration", True),
    "Arsenate reduction": ("detoxification_not_terminal_respiration", False),
    "Selenate reduction": ("potential_alternative_terminal_respiration", True),
    "Fumarate reduction": ("alternative_terminal_respiration", True),
    "Hydrogenases": ("hydrogen_metabolism", False),
}


###############################################################################
# FUNCTIONS
###############################################################################

def clean_mag_name(value):
    name = str(value).strip().split("\t")[0]
    name = Path(name).name
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", name, flags=re.I)


def read_eet_mag_lists():
    rows = []
    pattern = re.compile(r"^EET([1-4])_genome_names\.(txt|tsv)$", re.I)

    for path in DRAM_DIR.rglob("EET*_genome_names.*"):
        match = pattern.match(path.name)
        if not match:
            continue

        eet_group = f"EET{match.group(1)}"

        with path.open() as handle:
            names = set()

            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                name = clean_mag_name(line.split("\t")[0])
                if name.lower() in {"mag_name", "genome", "genome_name"}:
                    continue

                names.add(name)

        print(f"[INPUT] {eet_group}: {len(names)} names from {path}")

        for name in names:
            rows.append((name, eet_group))

    if not rows:
        raise RuntimeError(f"No EET genome-name files found under {DRAM_DIR}")

    membership_long = pd.DataFrame(rows, columns=["MAG_name", "EET_group"])
    membership_long = membership_long.drop_duplicates()

    membership = (
        membership_long
        .groupby("MAG_name", as_index=False)
        .agg(
            EET_groups=("EET_group", lambda x: ";".join(sorted(set(x)))),
            n_EET_groups=("EET_group", "nunique"),
        )
        .sort_values("MAG_name")
    )

    print(f"[OK] Unique EET MAGs: {len(membership)}")
    return membership


def marker_source(accession):
    if re.fullmatch(r"K\d{5}", accession, re.I):
        return "KEGG"
    if re.fullmatch(r"PF\d{5}(\.\d+)?", accession, re.I):
        return "PFAM"
    return "METABOLIC_HMM"


def normalize_accession(value, source):
    text = str(value).strip()

    if source == "KEGG":
        match = re.search(r"K\d{5}", text, re.I)
        return match.group(0).upper() if match else None

    if source == "PFAM":
        match = re.search(r"PF\d{5}", text, re.I)
        return match.group(0).upper() if match else None

    return re.sub(r"\.hmm$", "", Path(text).name, flags=re.I).lower()


def build_marker_table():
    rows = []

    for module, gene, accession in MARKERS:
        source = marker_source(accession)
        role, include = MODULE_ROLE[module]

        rows.append(
            {
                "metabolic_module": module,
                "gene_label": gene,
                "feature_source": source,
                "feature_accession": normalize_accession(accession, source),
                "analysis_class": role,
                "include_in_competing_respiration_test": include,
            }
        )

    table = pd.DataFrame(rows)

    # K00370 cannot distinguish nxrA from narG, and K00371 cannot distinguish
    # nxrB from narH. Collapse them so their expression is not counted twice.
    table = (
        table
        .groupby(
            [
                "metabolic_module",
                "feature_source",
                "feature_accession",
                "analysis_class",
                "include_in_competing_respiration_test",
            ],
            as_index=False,
            sort=False,
        )
        .agg(
            gene_label=(
                "gene_label",
                lambda x: "/".join(dict.fromkeys(x)),
            )
        )
    )

    return table


def find_hmm_file():
    path = Path(HMM_FILE)

    if not path.is_file():
        raise FileNotFoundError(path)

    return path


def normalize_series(series, source):
    values = series.astype("string").str.strip()

    if source == "KEGG":
        return values.str.extract(r"(K\d{5})", flags=re.I, expand=False).str.upper()

    if source == "PFAM":
        return values.str.extract(r"(PF\d{5})", flags=re.I, expand=False).str.upper()

    return values.str.replace(r"\.hmm$", "", regex=True, flags=re.I).str.lower()


def filter_profile(
    path,
    source,
    id_column,
    sample_column,
    marker_table,
    eet_mags,
    matched_mags,
    all_samples,
    mg_column=None,
):
    source_markers = marker_table[marker_table["feature_source"] == source].copy()
    targets = set(source_markers["feature_accession"])

    usecols = ["MAG_name", id_column, sample_column, "MT_coverage_per_cell"]
    if mg_column:
        usecols.append(mg_column)
    if source != "METABOLIC_HMM":
        usecols.append("MetaG_sample")

    output = []

    print(f"[READ] {source}: {path}")

    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        chunksize=CHUNK_SIZE,
        low_memory=False,
    ):
        chunk["MAG_name"] = chunk["MAG_name"].map(clean_mag_name)
        chunk = chunk[chunk["MAG_name"].isin(eet_mags)].copy()

        if chunk.empty:
            continue

        matched_mags.update(chunk["MAG_name"].unique())
        all_samples.update(chunk[sample_column].dropna().astype(str).unique())

        chunk["feature_accession"] = normalize_series(chunk[id_column], source)
        chunk = chunk[chunk["feature_accession"].isin(targets)].copy()

        if chunk.empty:
            continue

        chunk = chunk.merge(
            source_markers,
            on=["feature_source", "feature_accession"]
            if "feature_source" in chunk.columns
            else ["feature_accession"],
            how="inner",
        )

        chunk["feature_source"] = source
        chunk["MetaT_sample"] = chunk[sample_column]

        if source == "METABOLIC_HMM":
            chunk["MetaG_sample"] = chunk["MetaT_sample"].str.replace(
                r"_iMT$", "_iMG", regex=True
            )
            chunk["MG_coverage_per_cell"] = pd.NA
        else:
            chunk["MG_coverage_per_cell"] = chunk[mg_column]

        output.append(
            chunk[
                [
                    "MAG_name",
                    "metabolic_module",
                    "gene_label",
                    "feature_source",
                    "feature_accession",
                    "analysis_class",
                    "include_in_competing_respiration_test",
                    "MetaG_sample",
                    "MetaT_sample",
                    "MG_coverage_per_cell",
                    "MT_coverage_per_cell",
                ]
            ]
        )

    return pd.concat(output, ignore_index=True) if output else pd.DataFrame()



def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build sample-resolved expression tables for selected "
            "alternative pathways in EET-bearing MAGs."
        )
    )
    parser.add_argument(
        "eet_membership_directory",
        type=Path,
    )
    parser.add_argument(
        "kegg_expression_tsv",
        type=Path,
    )
    parser.add_argument(
        "pfam_expression_tsv",
        type=Path,
    )
    parser.add_argument(
        "metabolic_hmm_expression_tsv",
        type=Path,
    )
    parser.add_argument(
        "output_directory",
        type=Path,
    )
    return parser.parse_args()


###############################################################################
# MAIN
###############################################################################

def main():
    global DRAM_DIR
    global KEGG_FILE
    global PFAM_FILE
    global HMM_FILE
    global OUT_DIR

    args = parse_args()

    DRAM_DIR = args.eet_membership_directory
    KEGG_FILE = args.kegg_expression_tsv
    PFAM_FILE = args.pfam_expression_tsv
    HMM_FILE = args.metabolic_hmm_expression_tsv
    OUT_DIR = args.output_directory

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    membership = read_eet_mag_lists()
    membership.to_csv(
        OUT_DIR / "00_EET_MAG_membership.tsv",
        sep="\t",
        index=False,
    )

    eet_mags = set(membership["MAG_name"])
    group_map = membership.set_index("MAG_name")["EET_groups"].to_dict()
    markers = build_marker_table()

    # Add feature_source before merging in filter_profile.
    markers["feature_source"] = markers["feature_source"].astype(str)

    hmm_file = find_hmm_file()
    print(f"[INPUT] HMM table: {hmm_file}")

    matched_mags = set()
    all_samples = set()

    frames = []

    for path, source, id_col, sample_col, mg_col in [
        (KEGG_FILE, "KEGG", "KEGG_ko", "MetaT_sample", "MG_coverage_per_cell"),
        (PFAM_FILE, "PFAM", "PFAM_accession", "MetaT_sample", "MG_coverage_per_cell"),
        (hmm_file, "METABOLIC_HMM", "METABOLIC_hmm", "Sample", None),
    ]:
        if not Path(path).is_file():
            raise FileNotFoundError(path)

        frame = filter_profile(
            path=Path(path),
            source=source,
            id_column=id_col,
            sample_column=sample_col,
            marker_table=markers,
            eet_mags=eet_mags,
            matched_mags=matched_mags,
            all_samples=all_samples,
            mg_column=mg_col,
        )

        if not frame.empty:
            frames.append(frame)

    missing_from_profiles = sorted(eet_mags - matched_mags)
    (OUT_DIR / "QC_EET_MAGs_not_found_in_profiles.txt").write_text(
        "\n".join(missing_from_profiles) + ("\n" if missing_from_profiles else "")
    )


    expression = pd.concat(frames, ignore_index=True)
    expression["EET_groups"] = expression["MAG_name"].map(group_map)

    for col in ["MG_coverage_per_cell", "MT_coverage_per_cell"]:
        expression[col] = pd.to_numeric(expression[col], errors="coerce")

    keys = [
        "MAG_name",
        "EET_groups",
        "metabolic_module",
        "gene_label",
        "feature_source",
        "feature_accession",
        "analysis_class",
        "include_in_competing_respiration_test",
        "MetaG_sample",
        "MetaT_sample",
    ]

    marker_expression = (
        expression
        .groupby(keys, as_index=False, dropna=False)
        .agg(
            MG_coverage_per_cell=(
                "MG_coverage_per_cell",
                lambda x: x.sum(min_count=1),
            ),
            MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                lambda x: x.sum(min_count=1),
            ),
            n_profile_rows=("MT_coverage_per_cell", "size"),
        )
    )

    marker_expression["marker_is_expressed"] = (
        marker_expression["MT_coverage_per_cell"].fillna(0) > 0
    )

    marker_expression.to_csv(
        OUT_DIR / "01_selected_markers_expression_by_MAG_sample.tsv",
        sep="\t",
        index=False,
    )

    marker_expression["marker_id"] = (
        marker_expression["feature_source"]
        + ":"
        + marker_expression["feature_accession"]
    )
    marker_expression["expressed_marker_id"] = marker_expression["marker_id"].where(
        marker_expression["marker_is_expressed"]
    )

    module_keys = [
        "MAG_name",
        "EET_groups",
        "MetaG_sample",
        "MetaT_sample",
        "metabolic_module",
        "analysis_class",
        "include_in_competing_respiration_test",
    ]

    module_expression = (
        marker_expression
        .groupby(module_keys, as_index=False, dropna=False)
        .agg(
            n_markers_detected=("marker_id", "nunique"),
            n_markers_expressed=("expressed_marker_id", "nunique"),
            sum_MT_coverage_per_cell=(
                "MT_coverage_per_cell",
                lambda x: x.sum(min_count=1),
            ),
            mean_MT_coverage_per_cell=("MT_coverage_per_cell", "mean"),
            median_MT_coverage_per_cell=("MT_coverage_per_cell", "median"),
        )
    )

    module_expression["module_is_expressed"] = (
        module_expression["sum_MT_coverage_per_cell"].fillna(0) > 0
    )

    module_expression.to_csv(
        OUT_DIR / "02_selected_modules_expression_by_MAG_sample.tsv",
        sep="\t",
        index=False,
    )

    mags_with_markers = set(marker_expression["MAG_name"])
    no_markers = sorted(eet_mags - mags_with_markers)

    (OUT_DIR / "QC_EET_MAGs_without_selected_markers.txt").write_text(
        "\n".join(no_markers) + ("\n" if no_markers else "")
    )

    qc = [
        f"Unique EET MAGs in membership files: {len(eet_mags)}",
        f"EET MAGs matched in expression profiles: {len(matched_mags)}",
        f"MetaT samples observed: {len(all_samples)}",
        f"EET MAGs with at least one selected marker: {len(mags_with_markers)}",
        f"EET MAGs without selected markers: {len(no_markers)}",
        f"Marker-level rows: {len(marker_expression)}",
        f"Module-level rows: {len(module_expression)}",
    ]

    (OUT_DIR / "QC_summary.txt").write_text("\n".join(qc) + "\n")

    print("\n[COMPLETED]")
    for line in qc:
        print(f"  {line}")

    print(f"\nOutputs: {OUT_DIR}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"\n[ERROR] {error}", file=sys.stderr)
        sys.exit(1)