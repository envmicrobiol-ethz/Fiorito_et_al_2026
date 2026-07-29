#!/usr/bin/env python3

# DESCRIPTION
# Assigns cluster-level credibility scores to candidate MHC/EET protein
# families using motif content, subcellular localization, PFAM support,
# heme density, and metatranscriptomic expression.
#
# INPUT
# A member-level TSV or CSV table containing MMseqs2 cluster assignments and,
# where available, motif counts, protein lengths, PSORTb localization, PFAM
# annotations, and expression values.
#
# OUTPUT
# A cluster-level TSV table with component metrics, credibility scores, and
# Tier 1-3 assignments, plus a JSON file recording the scoring rules.
#
# USAGE
# python3 04_score_MHC_EET_credibility.py \
#   --input MHC_cluster_members.tsv \
#   --output MHC_EET_credibility_scores.tsv

"""
Cluster-level credibility scoring for multiheme cytochrome (MHC) / extracellular
electron transfer (EET) candidates.

The script scores candidate protein clusters using:
  - multiheme motif content,
  - protein localization,
  - PFAM support for c-type cytochromes / EET-like proteins,
  - heme density,
  - expression support,
  - penalties for likely false-positive PFAMs or cytoplasmic localization.

Input:
  A member-level table where each row is a protein/member assigned to a cluster.
  The table must contain a cluster identifier column, by default: SeqC_ID.

Output:
  A cluster-level table with credibility scores and tiers.

Example:
  python3 04_score_MHC_EET_credibility.py \
    --input MHC_cluster_members.tsv \
    --output MHC_credibility_scores.tsv
    
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# 1) DEFAULT SCORING PARAMETERS
# ============================================================

CORE_CYTC_PF = {
    "PF00034",  # Cytochrome c-like domain
    "PF02335",  # Cytochrome c552
    "PF13435",  # Cytochrome c554 / c'-type
    "PF14537",  # Tetrahaem cytochrome
    "PF11783",  # Octaheme c-type cytochrome
    "PF02085",  # Class III cytochrome c
    "PF14522",  # Cytochrome c7-like
    "PF07635",  # Planctomycete-type cytochrome c
    "PF22678",  # NrfB-like c-type cytochrome domain
    "PF09698",  # Geobacter / double-motif cytochrome signature
    "PF09699",  # Geobacter / double-motif cytochrome signature
    "PF22112",  # OmcA-like N-terminus
    "PF22113",  # MtrC/MtrF-like domains II/IV
}

OM_EET_PF = {
    "PF22112",  # OmcA-like N-terminus
    "PF22113",  # MtrC/MtrF-like domains II/IV
    "PF03264",  # NapC/NirT cytochrome c, N-terminal
}

NEGATIVE_PF = {
    # Fe-S / redox scaffolds
    "PF00037", "PF10588", "PF12801", "PF12838", "PF13183", "PF13510", "PF22117",

    # Molybdopterin oxidoreductase-related
    "PF00384", "PF04879",

    # Radical SAM
    "PF04055",

    # ABC/UvrA
    "PF00005", "PF17755", "PF17760",

    # tRNA synthetase-like
    "PF00133", "PF08264",
}

STRONG_ENVELOPE_LOCALIZATIONS = {
    "Periplasmic",
    "OuterMembrane",
    "Extracellular",
}

WEAK_ENVELOPE_LOCALIZATIONS = {
    "CytoplasmicMembrane",
}


# ============================================================
# 2) HELPER FUNCTIONS
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score MHC/EET candidate clusters using heme motifs, localization, "
            "PFAM evidence, heme density and expression."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input member-level table. TSV or CSV. Must contain SeqC_ID by default.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output cluster-level scoring table. Recommended extension: .tsv",
    )

    parser.add_argument(
        "--metadata-output",
        default=None,
        help=(
            "Optional JSON file recording selected columns and scoring parameters. "
            "Default: <output>.metadata.json"
        ),
    )

    parser.add_argument(
        "--sep",
        default="auto",
        choices=["auto", "tab", "comma"],
        help="Input separator. Default: auto from file extension.",
    )

    parser.add_argument(
        "--cluster-col",
        default="SeqC_ID",
        help="Cluster identifier column. Default: SeqC_ID.",
    )

    parser.add_argument(
        "--member-col",
        default="Member",
        help="Protein/member identifier column. Default: Member.",
    )

    parser.add_argument(
        "--heme-col",
        default=None,
        help=(
            "Column with the heme motif count used for scoring. "
            "If omitted, the script uses sum_heme_motifs_member when present. "
            "Otherwise, it uses the maximum count across the overlapping "
            "CxxCH, CxxCnH/CxxC_H and CXnCH motif categories."
        ),
    )

    parser.add_argument(
        "--length-col",
        default=None,
        help=(
            "Protein length column. If omitted, the script searches for "
            "seq_length_x, seq_length, seq_length_x_y, seq_length_x_x."
        ),
    )

    parser.add_argument(
        "--expression-col",
        default=None,
        help=(
            "Expression column. If omitted, the script searches for "
            "SUM_TPM_member, then SUM_TPM. If none is found, expression is set to 0."
        ),
    )

    parser.add_argument(
        "--localization-col",
        default=None,
        help=(
            "Localization column. If omitted, the script searches for "
            "PSORT_Localization_member, then PSORT_Localization_representative."
        ),
    )

    return parser.parse_args()


def infer_separator(path: Path, sep_arg: str) -> str:
    if sep_arg == "tab":
        return "\t"
    if sep_arg == "comma":
        return ","

    lower_name = path.name.lower()

    if lower_name.endswith(".csv") or lower_name.endswith(".csv.gz"):
        return ","

    return "\t"


def safe_numeric(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def pick_first_existing_column(
    df: pd.DataFrame,
    candidates: list[str],
    user_choice: str | None = None,
) -> str | None:
    if user_choice:
        if user_choice not in df.columns:
            raise ValueError(f"Requested column not found: {user_choice}")
        return user_choice

    for col in candidates:
        if col in df.columns:
            return col

    return None


def clean_localization(series: pd.Series) -> pd.Series:
    cleaned = series.where(series.notna(), "Unknown").astype(str).str.strip()
    cleaned = cleaned.replace({"": "Unknown", "nan": "Unknown", "None": "Unknown"})
    return cleaned


def extract_pfam_accessions(row: pd.Series, pfam_cols: list[str]) -> set[str]:
    """
    Extract PFAM accessions from pfam_* columns.

    Expected format:
      PF14537 | Tetrahaem cytochrome ...
    """
    pfams = set()

    for col in pfam_cols:
        value = row.get(col, "")

        if pd.notna(value) and str(value).strip():
            accession = str(value).split("|")[0].strip()

            if accession:
                pfams.add(accession)

    return pfams


def make_metadata(
    args: argparse.Namespace,
    input_rows: int,
    output_rows: int,
    selected_columns: dict,
) -> dict:
    return {
        "script": "04_score_MHC_EET_credibility.py",
        "description": "Cluster-level credibility scoring for MHC/EET candidates.",
        "input_file": str(args.input),
        "output_file": str(args.output),
        "input_rows": input_rows,
        "output_rows": output_rows,
        "selected_columns": selected_columns,
        "scoring_rules": {
            "multiheme": {
                "+2": "p_heme6 >= 0.5",
                "+1": "median_hemes >= 6",
                "+1_additional": "p_heme10 >= 0.2",
            },
            "localization": {
                "+2": "p_env_strong >= 0.4 and p_cyt <= 0.2",
                "+1": "p_env_strong < 0.4 and p_env_weak >= 0.5 and p_cyt <= 0.2",
                "-2": "p_cyt >= 0.5",
            },
            "pfam": {
                "+2": "p_core >= 0.5",
                "+1": "p_om >= 0.2",
                "-2": "p_negative >= 0.5",
            },
            "density": {
                "+1": "median_density >= 2 motifs per 100 aa",
                "+1_additional": "median_density >= 4 motifs per 100 aa",
            },
            "expression": {
                "+1": "max_TPM_member > 0 and >= 75th percentile of cluster max_TPM_member",
            },
            "tiers": {
                "Tier1_credible": "score >= 6",
                "Tier2_putative": "score 3 to 5",
                "Tier3_uncertain": "score <= 2",
            },
        },
        "pfam_sets": {
            "core_cytochrome_c_or_multiheme_pfams": sorted(CORE_CYTC_PF),
            "outer_membrane_eet_pfams": sorted(OM_EET_PF),
            "negative_pfams": sorted(NEGATIVE_PF),
        },
        "localization_sets": {
            "strong_envelope": sorted(STRONG_ENVELOPE_LOCALIZATIONS),
            "weak_envelope": sorted(WEAK_ENVELOPE_LOCALIZATIONS),
        },
    }


# ============================================================
# 3) MAIN SCORING FUNCTION
# ============================================================

def score_clusters(df: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, dict]:
    d = df.copy()

    cluster_col = args.cluster_col

    if cluster_col not in d.columns:
        raise ValueError(
            f"Missing required cluster column: {cluster_col}. "
            "Use --cluster-col if your cluster column has a different name."
        )

    member_col = args.member_col if args.member_col in d.columns else None

    length_col = pick_first_existing_column(
        d,
        ["seq_length_x", "seq_length", "seq_length_x_y", "seq_length_x_x"],
        user_choice=args.length_col,
    )

    expression_col = pick_first_existing_column(
        d,
        ["SUM_TPM_member", "SUM_TPM"],
        user_choice=args.expression_col,
    )

    localization_col = pick_first_existing_column(
        d,
        ["PSORT_Localization_member", "PSORT_Localization_representative"],
        user_choice=args.localization_col,
    )

    # -----------------------------
    # Heme motif column
    # -----------------------------

    if args.heme_col:
        if args.heme_col not in d.columns:
            raise ValueError(f"Requested heme motif column not found: {args.heme_col}")
        heme_col = args.heme_col

    elif "sum_heme_motifs_member" in d.columns:
        heme_col = "sum_heme_motifs_member"

    else:
        # The three motif definitions overlap. In particular, a canonical
        # CXXCH motif can also match the broader CXXCnH and CXnCH patterns.
        # Therefore, summing the three columns can count the same motif more
        # than once. The fallback uses the largest same-type motif count.
        motif_aliases = {
            "CxxCH": ["num_CxxCH_member"],
            "CxxCnH": ["num_CxxCnH_member", "num_CxxC_H_member"],
            "CXnCH": ["num_CXnCH_member"],
        }

        selected_motif_cols = []

        for aliases in motif_aliases.values():
            selected = next((col for col in aliases if col in d.columns), None)
            if selected is not None:
                selected_motif_cols.append(selected)

        if not selected_motif_cols:
            print(
                "[WARNING] No heme motif columns found. Heme counts will be set to 0.",
                file=sys.stderr,
            )
            d["max_same_type_heme_motifs_member_fallback"] = 0
        else:
            motif_counts = pd.concat(
                [
                    safe_numeric(d[col]).fillna(0).rename(col)
                    for col in selected_motif_cols
                ],
                axis=1,
            )

            d["max_same_type_heme_motifs_member_fallback"] = motif_counts.max(axis=1)

            print(
                "[INFO] Using the maximum same-type motif count across: "
                + ", ".join(selected_motif_cols),
                file=sys.stderr,
            )

        heme_col = "max_same_type_heme_motifs_member_fallback"

    # -----------------------------
    # Expression
    # -----------------------------

    expression_available = expression_col is not None

    if expression_col is None:
        print(
            "[WARNING] No expression column found. Expression values will be set to 0 "
            "and no expression bonus will be assigned.",
            file=sys.stderr,
        )
        d["__expression__"] = 0
        expression_col = "__expression__"

    # -----------------------------
    # PFAM columns
    # -----------------------------

    pfam_cols = [col for col in d.columns if col.startswith("pfam_")]

    if not pfam_cols:
        print(
            "[WARNING] No pfam_* columns found. PFAM support will be set to False.",
            file=sys.stderr,
        )

    # -----------------------------
    # Numeric cleaning
    # -----------------------------

    d[heme_col] = safe_numeric(d[heme_col]).fillna(0)

    if length_col:
        d[length_col] = safe_numeric(d[length_col])
    else:
        print(
            "[WARNING] No sequence length column found. Heme density bonuses cannot be assigned.",
            file=sys.stderr,
        )

    d[expression_col] = safe_numeric(d[expression_col]).fillna(0)

    # -----------------------------
    # Member-level flags
    # -----------------------------

    d["heme6"] = d[heme_col] >= 6
    d["heme10"] = d[heme_col] >= 10

    if length_col:
        d = d[d[length_col].fillna(0) > 0].copy()

        if d.empty:
            raise ValueError(
                "No rows remain after filtering for positive sequence length."
            )

        d["heme_density_100aa"] = 100 * d[heme_col] / d[length_col]

    else:
        d["heme_density_100aa"] = np.nan

    d["dense4"] = d["heme_density_100aa"] >= 4

    if localization_col:
        d["loc_clean"] = clean_localization(d[localization_col])
    else:
        d["loc_clean"] = "Unknown"

    d["is_env_strong"] = d["loc_clean"].isin(STRONG_ENVELOPE_LOCALIZATIONS)
    d["is_env_weak"] = d["loc_clean"].isin(WEAK_ENVELOPE_LOCALIZATIONS)
    d["is_cytoplasmic"] = d["loc_clean"].eq("Cytoplasmic")
    d["is_unknown_loc"] = d["loc_clean"].eq("Unknown")

    if pfam_cols:
        d["_pfam_set"] = d.apply(
            lambda row: extract_pfam_accessions(row, pfam_cols),
            axis=1,
        )

        d["has_core_cytc"] = d["_pfam_set"].apply(
            lambda pfams: len(pfams & CORE_CYTC_PF) > 0
        )
        d["has_om_eet"] = d["_pfam_set"].apply(
            lambda pfams: len(pfams & OM_EET_PF) > 0
        )
        d["has_negative"] = d["_pfam_set"].apply(
            lambda pfams: len(pfams & NEGATIVE_PF) > 0
        )

    else:
        d["has_core_cytc"] = False
        d["has_om_eet"] = False
        d["has_negative"] = False

    # -----------------------------
    # Cluster-level aggregation
    # -----------------------------

    if member_col:
        n_members_agg = (member_col, "nunique")
    else:
        d["__row_count__"] = 1
        n_members_agg = ("__row_count__", "sum")

    cluster_scores = (
        d.groupby(cluster_col, dropna=False)
        .agg(
            n_members=n_members_agg,

            p_heme6=("heme6", "mean"),
            p_heme10=("heme10", "mean"),
            median_hemes=(heme_col, "median"),
            max_hemes=(heme_col, "max"),
            median_density=("heme_density_100aa", "median"),
            p_dense4=("dense4", "mean"),

            p_env_strong=("is_env_strong", "mean"),
            p_env_weak=("is_env_weak", "mean"),
            p_cyt=("is_cytoplasmic", "mean"),
            p_unknown=("is_unknown_loc", "mean"),

            p_core=("has_core_cytc", "mean"),
            p_om=("has_om_eet", "mean"),
            p_negative=("has_negative", "mean"),

            SUM_TPM_cluster=(expression_col, "sum"),
            mean_TPM_member=(expression_col, "mean"),
            median_TPM_member=(expression_col, "median"),
            max_TPM_member=(expression_col, "max"),
        )
        .reset_index()
        .rename(columns={cluster_col: "SeqC_ID"})
    )

    # -----------------------------
    # Scoring
    # -----------------------------

    cluster_scores["score"] = 0

    cluster_scores.loc[cluster_scores["p_heme6"] >= 0.5, "score"] += 2
    cluster_scores.loc[cluster_scores["median_hemes"] >= 6, "score"] += 1
    cluster_scores.loc[cluster_scores["p_heme10"] >= 0.2, "score"] += 1

    cluster_scores.loc[
        (cluster_scores["p_env_strong"] >= 0.4)
        & (cluster_scores["p_cyt"] <= 0.2),
        "score",
    ] += 2

    cluster_scores.loc[
        (cluster_scores["p_env_strong"] < 0.4)
        & (cluster_scores["p_env_weak"] >= 0.5)
        & (cluster_scores["p_cyt"] <= 0.2),
        "score",
    ] += 1

    cluster_scores.loc[cluster_scores["p_core"] >= 0.5, "score"] += 2
    cluster_scores.loc[cluster_scores["p_om"] >= 0.2, "score"] += 1

    cluster_scores.loc[cluster_scores["median_density"] >= 2, "score"] += 1
    cluster_scores.loc[cluster_scores["median_density"] >= 4, "score"] += 1

    q75_max_tpm = cluster_scores["max_TPM_member"].quantile(0.75)

    if expression_available:
        cluster_scores.loc[
            (cluster_scores["max_TPM_member"] > 0)
            & (cluster_scores["max_TPM_member"] >= q75_max_tpm),
            "score",
        ] += 1

    cluster_scores.loc[cluster_scores["p_negative"] >= 0.5, "score"] -= 2
    cluster_scores.loc[cluster_scores["p_cyt"] >= 0.5, "score"] -= 2

    # -----------------------------
    # Tiering
    # -----------------------------

    cluster_scores["tier"] = np.select(
        [
            cluster_scores["score"] >= 6,
            cluster_scores["score"].between(3, 5, inclusive="both"),
            cluster_scores["score"] <= 2,
        ],
        [
            "Tier1_credible",
            "Tier2_putative",
            "Tier3_uncertain",
        ],
        default="Tier3_uncertain",
    )

    tier_rank = {
        "Tier1_credible": 1,
        "Tier2_putative": 2,
        "Tier3_uncertain": 3,
    }

    cluster_scores["tier_rank"] = cluster_scores["tier"].map(tier_rank)

    cluster_scores = cluster_scores.sort_values(
        ["tier_rank", "score", "max_TPM_member", "n_members"],
        ascending=[True, False, False, False],
    ).drop(columns=["tier_rank"])

    selected_columns = {
        "cluster_col": cluster_col,
        "member_col": member_col,
        "heme_col": heme_col,
        "length_col": length_col,
        "expression_col": expression_col,
        "expression_available": expression_available,
        "expression_q75_max_tpm": float(q75_max_tpm),
        "localization_col": localization_col,
        "pfam_cols": pfam_cols,
    }

    return cluster_scores, selected_columns


# ============================================================
# 4) RUN
# ============================================================

def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    sep = infer_separator(input_path, args.sep)

    df = pd.read_csv(
        input_path,
        sep=sep,
        dtype=str,
        low_memory=False,
    )

    print(f"Input file: {input_path}")
    print(f"Input rows: {len(df)}")
    print(f"Input columns: {len(df.columns)}")

    cluster_scores, selected_columns = score_clusters(df, args)

    cluster_scores.to_csv(output_path, sep="\t", index=False)

    if args.metadata_output:
        metadata_path = Path(args.metadata_output)
    else:
        metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")

    metadata = make_metadata(
        args=args,
        input_rows=len(df),
        output_rows=len(cluster_scores),
        selected_columns=selected_columns,
    )

    with open(metadata_path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    print("\nScoring completed.")
    print(f"Clusters scored: {len(cluster_scores)}")
    print(f"Output table: {output_path}")
    print(f"Metadata JSON: {metadata_path}")

    print("\nTier counts:")
    print(cluster_scores["tier"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
