#!/usr/bin/env python3

# DESCRIPTION
# Calls credible external-dataset matches to conserved EET gene cassettes from
# TBLASTN results using an explicit query-gene-to-cassette map.
#
# Default protein-level filters:
#   E-value <= 1e-10
#   amino-acid identity >= 50%
#   query coverage >= 90%
#
# TBLASTN HSPs on the same target strand are clustered into target loci when
# they overlap by at least 50% of the shorter interval. A credible cassette
# match requires, on the same target contig:
#   1. one essential EET query gene;
#   2. one different query gene from the same explicit query cassette; and
#   3. the two genes must map to distinct target loci.
#
# INPUT
# 1. Raw TBLASTN table from 02_run_TBLASTN_external_MAGs.sh.
# 2. Explicit query map from
#    03_build_explicit_EET_query_cassette_map.py.
# 3. Contig-to-MAG mapping from:
#    ../01_BLASTn_external_MAGs/02_build_contig_to_MAG_mapping.py
#
# OUTPUT
# Filtered HSPs, target-locus associations, cassette candidates, credible
# cassette calls, a target-MAG list for GTDB-Tk and a summary table.
#
# USAGE
# python3 04_call_credible_EET_cassettes_from_TBLASTN.py \
#   --tblastn-tsv SPRUCE_EET_proteins_TBLASTN.tsv \
#   --dataset SPRUCE \
#   --query-cassette-map EET_query_gene_to_cassette.tsv \
#   --contig-to-mag SPRUCE_contig_to_MAG.tsv \
#   --output-dir SPRUCE_TBLASTN_credible_calls

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]

NUMERIC_COLUMNS = [
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Call credible EET cassette matches from TBLASTN using an "
            "explicit query-gene-to-cassette map."
        )
    )
    parser.add_argument("--tblastn-tsv", required=True, type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--query-cassette-map", required=True, type=Path)
    parser.add_argument("--contig-to-mag", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-pident", type=float, default=50.0)
    parser.add_argument("--min-qcov", type=float, default=90.0)
    parser.add_argument("--max-evalue", type=float, default=1e-10)
    parser.add_argument(
        "--same-locus-overlap",
        type=float,
        default=0.50,
        help=(
            "Minimum overlap as a fraction of the shorter target interval "
            "for two HSPs to be assigned to the same target locus. "
            "Default: 0.50."
        ),
    )
    parser.add_argument(
        "--allow-unmapped-contigs",
        action="store_true",
        help=(
            "Drop filtered hits whose target contig is absent from the "
            "contig-to-MAG mapping. By default, these cause an error."
        ),
    )
    return parser.parse_args()


def require_file(path: Path, label: str, allow_empty: bool = False) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} is missing: {path}")

    if not allow_empty and path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is empty: {path}")


def parse_boolean(value: object) -> bool:
    text = str(value).strip().casefold()

    if text in {"true", "1", "yes", "y"}:
        return True

    if text in {"false", "0", "no", "n"}:
        return False

    raise ValueError(f"Cannot interpret Boolean value: {value!r}")


def normalize_target_contig(sseqid: object) -> str:
    token = str(sseqid).strip().split()[0]

    accession = re.search(r"([A-Za-z0-9_]+\.[0-9]+)", token)
    return accession.group(1) if accession else token.strip("|")


def overlap_fraction_of_shorter(
    start_a: int,
    end_a: int,
    start_b: int,
    end_b: int,
) -> float:
    overlap = max(
        0,
        min(end_a, end_b) - max(start_a, start_b) + 1,
    )

    if overlap == 0:
        return 0.0

    length_a = end_a - start_a + 1
    length_b = end_b - start_b + 1

    return overlap / min(length_a, length_b)


def cluster_hits_into_target_loci(
    group: pd.DataFrame,
    overlap_threshold: float,
) -> pd.DataFrame:
    rows = (
        group.sort_values(
            ["bitscore", "pident", "qcov", "evalue"],
            ascending=[False, False, False, True],
        )
        .to_dict("records")
    )

    loci: list[list[dict[str, object]]] = []

    for row in rows:
        matching_loci: list[int] = []

        for locus_index, locus in enumerate(loci):
            if any(
                overlap_fraction_of_shorter(
                    int(row["target_start"]),
                    int(row["target_end"]),
                    int(member["target_start"]),
                    int(member["target_end"]),
                )
                >= overlap_threshold
                for member in locus
            ):
                matching_loci.append(locus_index)

        if not matching_loci:
            loci.append([row])
            continue

        first_index = matching_loci[0]
        loci[first_index].append(row)

        for other_index in reversed(matching_loci[1:]):
            loci[first_index].extend(loci[other_index])
            del loci[other_index]

    output_tables: list[pd.DataFrame] = []

    loci.sort(
        key=lambda locus: min(
            int(row["target_start"])
            for row in locus
        )
    )

    for locus_number, locus in enumerate(loci, start=1):
        locus_start = min(
            int(row["target_start"])
            for row in locus
        )
        locus_end = max(
            int(row["target_end"])
            for row in locus
        )
        strand = str(locus[0]["target_strand"])
        target_contig = str(locus[0]["target_contig"])

        locus_id = (
            f"{target_contig}:{locus_start}-{locus_end}:"
            f"{strand}:L{locus_number:03d}"
        )

        locus_table = (
            pd.DataFrame(locus)
            .sort_values(
                ["bitscore", "pident", "qcov", "evalue"],
                ascending=[False, False, False, True],
            )
            .drop_duplicates(
                subset=["qseqid"],
                keep="first",
            )
        )

        locus_table["target_locus_id"] = locus_id
        locus_table["target_locus_start"] = locus_start
        locus_table["target_locus_end"] = locus_end
        output_tables.append(locus_table)

    if not output_tables:
        return pd.DataFrame()

    return pd.concat(output_tables, ignore_index=True)


def find_credible_witness_pair(
    group: pd.DataFrame,
) -> dict[str, object]:
    associations = (
        group.sort_values(
            ["bitscore", "pident", "qcov", "evalue"],
            ascending=[False, False, False, True],
        )
        .drop_duplicates(
            subset=["query_gene", "target_locus_id"]
        )
    )

    essential_hits = associations[
        associations["is_essential"]
    ]

    for _, essential_hit in essential_hits.iterrows():
        additional_hits = associations[
            (
                associations["query_gene"]
                != essential_hit["query_gene"]
            )
            & (
                associations["target_locus_id"]
                != essential_hit["target_locus_id"]
            )
        ]

        if additional_hits.empty:
            continue

        additional_hit = additional_hits.iloc[0]

        return {
            "credible_match": True,
            "witness_essential_gene": essential_hit["query_gene"],
            "witness_essential_locus": essential_hit["target_locus_id"],
            "witness_additional_gene": additional_hit["query_gene"],
            "witness_additional_locus": additional_hit["target_locus_id"],
        }

    return {
        "credible_match": False,
        "witness_essential_gene": "",
        "witness_essential_locus": "",
        "witness_additional_gene": "",
        "witness_additional_locus": "",
    }


def load_contig_mapping(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path, sep="\t", dtype=str)

    if {"target_contig", "target_MAG"}.issubset(mapping.columns):
        selected = mapping.copy()
    elif {
        "contig_accession",
        "target_mag_file",
    }.issubset(mapping.columns):
        selected = mapping.rename(
            columns={
                "contig_accession": "target_contig",
                "target_mag_file": "target_MAG",
            }
        )
    else:
        raise ValueError(
            "Contig-to-MAG table must contain either target_contig and "
            "target_MAG, or the legacy columns contig_accession and "
            "target_mag_file."
        )

    ambiguous = (
        selected.groupby("target_contig")["target_MAG"]
        .nunique()
    )
    ambiguous = ambiguous[ambiguous > 1]

    if not ambiguous.empty:
        raise ValueError(
            "One or more target contigs map to multiple MAGs."
        )

    return selected.drop_duplicates(
        subset=["target_contig"]
    )


def write_table(
    path: Path,
    rows: pd.DataFrame,
    columns: list[str],
) -> None:
    if rows.empty:
        pd.DataFrame(columns=columns).to_csv(
            path,
            sep="\t",
            index=False,
        )
    else:
        rows.to_csv(path, sep="\t", index=False)


def main() -> None:
    args = parse_args()

    if not args.dataset.strip():
        raise ValueError("--dataset cannot be empty.")

    if not 0 <= args.same_locus_overlap <= 1:
        raise ValueError("--same-locus-overlap must be between 0 and 1.")

    tblastn_path = args.tblastn_tsv.resolve()
    query_map_path = args.query_cassette_map.resolve()
    contig_map_path = args.contig_to_mag.resolve()
    output_dir = args.output_dir.resolve()

    require_file(tblastn_path, "TBLASTN table", allow_empty=True)
    require_file(query_map_path, "Query cassette map")
    require_file(contig_map_path, "Contig-to-MAG mapping")

    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory exists and is not empty: {output_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    query_map = pd.read_csv(
        query_map_path,
        sep="\t",
        dtype=str,
    )

    required_query_columns = {
        "qseqid",
        "query_MAG",
        "query_gene",
        "query_contig",
        "gene_number",
        "cassette_id",
        "is_essential",
    }

    missing_query_columns = (
        required_query_columns - set(query_map.columns)
    )

    if missing_query_columns:
        raise ValueError(
            "Query map is missing columns: "
            + ", ".join(sorted(missing_query_columns))
        )

    query_map["is_essential"] = (
        query_map["is_essential"].map(parse_boolean)
    )

    if query_map["qseqid"].duplicated().any():
        duplicates = (
            query_map.loc[
                query_map["qseqid"].duplicated(keep=False),
                "qseqid",
            ]
            .drop_duplicates()
            .head(10)
            .tolist()
        )
        raise ValueError(
            "Duplicate qseqid values in query map: "
            + ", ".join(duplicates)
        )

    contig_map = load_contig_mapping(contig_map_path)

    try:
        blast = pd.read_csv(
            tblastn_path,
            sep="\t",
            header=None,
            names=BLAST_COLUMNS,
            low_memory=False,
        )
    except EmptyDataError:
        blast = pd.DataFrame(columns=BLAST_COLUMNS)

    filtered_path = (
        output_dir
        / f"{args.dataset}__tblastn_filtered_parsed.tsv"
    )
    locus_path = (
        output_dir
        / f"{args.dataset}__tblastn_target_locus_associations.tsv"
    )
    candidate_path = (
        output_dir
        / f"{args.dataset}__all_cassette_candidates.tsv"
    )
    calls_path = (
        output_dir
        / f"{args.dataset}__credible_EET_cassette_calls.tsv"
    )
    mag_list_path = (
        output_dir
        / f"{args.dataset}__target_MAGs_for_GTDBTk.list"
    )
    summary_path = (
        output_dir
        / f"{args.dataset}__summary_counts.tsv"
    )

    if blast.empty:
        write_table(filtered_path, blast, BLAST_COLUMNS)
        write_table(locus_path, pd.DataFrame(), [])
        write_table(candidate_path, pd.DataFrame(), [])
        write_table(calls_path, pd.DataFrame(), [])
        mag_list_path.write_text("", encoding="utf-8")

        pd.DataFrame(
            [
                {
                    "dataset": args.dataset,
                    "n_raw_tblastn_HSPs": 0,
                    "n_filtered_HSPs": 0,
                    "n_explicit_query_cassettes": (
                        query_map["cassette_id"].nunique()
                    ),
                    "n_target_loci": 0,
                    "n_credible_cassette_calls": 0,
                    "n_target_MAGs": 0,
                    "n_target_contigs": 0,
                    "min_pident": args.min_pident,
                    "min_qcov": args.min_qcov,
                    "max_evalue": args.max_evalue,
                    "same_locus_overlap": args.same_locus_overlap,
                }
            ]
        ).to_csv(summary_path, sep="\t", index=False)

        print(f"No TBLASTN hits were reported for {args.dataset}.")
        return

    for column in NUMERIC_COLUMNS:
        blast[column] = pd.to_numeric(
            blast[column],
            errors="raise",
        )

    blast["qseqid"] = (
        blast["qseqid"]
        .astype(str)
        .str.split()
        .str[0]
    )
    blast["dataset"] = args.dataset

    blast = blast.merge(
        query_map,
        on="qseqid",
        how="left",
        validate="many_to_one",
    )

    if blast["cassette_id"].isna().any():
        missing_ids = (
            blast.loc[
                blast["cassette_id"].isna(),
                "qseqid",
            ]
            .drop_duplicates()
            .head(20)
            .tolist()
        )
        raise ValueError(
            "TBLASTN query IDs absent from the explicit cassette map: "
            + ", ".join(missing_ids)
        )

    blast["query_aligned_span_aa"] = (
        blast["qend"] - blast["qstart"]
    ).abs() + 1

    blast["qcov"] = (
        blast["query_aligned_span_aa"]
        / blast["qlen"]
        * 100.0
    )

    blast["target_contig"] = (
        blast["sseqid"].map(normalize_target_contig)
    )
    blast["target_start"] = (
        blast[["sstart", "send"]]
        .min(axis=1)
        .astype(int)
    )
    blast["target_end"] = (
        blast[["sstart", "send"]]
        .max(axis=1)
        .astype(int)
    )
    blast["target_strand"] = (
        blast["sstart"] <= blast["send"]
    ).map({True: "+", False: "-"})

    blast = blast.merge(
        contig_map,
        on="target_contig",
        how="left",
        validate="many_to_one",
    )

    filtered = blast[
        (blast["pident"] >= args.min_pident)
        & (blast["qcov"] >= args.min_qcov)
        & (blast["evalue"] <= args.max_evalue)
    ].copy()

    filtered.to_csv(
        filtered_path,
        sep="\t",
        index=False,
    )

    unmapped = filtered["target_MAG"].isna()

    if unmapped.any() and not args.allow_unmapped_contigs:
        examples = (
            filtered.loc[unmapped, "target_contig"]
            .drop_duplicates()
            .head(10)
            .tolist()
        )
        raise ValueError(
            f"{int(unmapped.sum())} filtered hits could not be mapped "
            "to a target MAG. Examples: "
            + ", ".join(examples)
        )

    mapped = filtered.loc[~unmapped].copy()

    locus_group_columns = [
        "cassette_id",
        "query_MAG",
        "query_contig",
        "dataset",
        "target_MAG",
        "target_contig",
        "target_strand",
    ]

    locus_tables: list[pd.DataFrame] = []

    for _, group in mapped.groupby(
        locus_group_columns,
        dropna=False,
        sort=True,
    ):
        clustered = cluster_hits_into_target_loci(
            group,
            overlap_threshold=args.same_locus_overlap,
        )

        if not clustered.empty:
            locus_tables.append(clustered)

    locus_hits = (
        pd.concat(locus_tables, ignore_index=True)
        if locus_tables
        else pd.DataFrame()
    )

    write_table(
        locus_path,
        locus_hits,
        list(filtered.columns)
        + [
            "target_locus_id",
            "target_locus_start",
            "target_locus_end",
        ],
    )

    call_group_columns = [
        "cassette_id",
        "query_MAG",
        "query_contig",
        "dataset",
        "target_MAG",
        "target_contig",
    ]

    candidate_rows: list[dict[str, object]] = []

    if not locus_hits.empty:
        for keys, group in locus_hits.groupby(
            call_group_columns,
            dropna=False,
            sort=True,
        ):
            (
                cassette_id,
                query_mag,
                query_contig,
                dataset,
                target_mag,
                target_contig,
            ) = keys

            associations = group.drop_duplicates(
                subset=["query_gene", "target_locus_id"]
            )

            query_genes = sorted(
                associations["query_gene"]
                .astype(str)
                .unique()
            )
            essential_genes = sorted(
                associations.loc[
                    associations["is_essential"],
                    "query_gene",
                ]
                .astype(str)
                .unique()
            )

            witness = find_credible_witness_pair(associations)

            candidate_rows.append(
                {
                    "cassette_id": cassette_id,
                    "query_MAG": query_mag,
                    "query_contig": query_contig,
                    "dataset": dataset,
                    "target_MAG": target_mag,
                    "target_contig": target_contig,
                    "n_total_genes": len(query_genes),
                    "n_essential_genes": len(essential_genes),
                    "n_target_loci": (
                        associations["target_locus_id"].nunique()
                    ),
                    "query_genes": ",".join(query_genes),
                    "essential_genes": ",".join(essential_genes),
                    **witness,
                    "target_region_start": int(
                        associations["target_locus_start"].min()
                    ),
                    "target_region_end": int(
                        associations["target_locus_end"].max()
                    ),
                    "minimum_pident": associations["pident"].min(),
                    "minimum_qcov": associations["qcov"].min(),
                    "best_evalue": associations["evalue"].min(),
                    "best_bitscore": associations["bitscore"].max(),
                }
            )

    candidates = pd.DataFrame(candidate_rows)

    candidate_columns = [
        "cassette_id",
        "query_MAG",
        "query_contig",
        "dataset",
        "target_MAG",
        "target_contig",
        "n_total_genes",
        "n_essential_genes",
        "n_target_loci",
        "query_genes",
        "essential_genes",
        "credible_match",
        "witness_essential_gene",
        "witness_essential_locus",
        "witness_additional_gene",
        "witness_additional_locus",
        "target_region_start",
        "target_region_end",
        "minimum_pident",
        "minimum_qcov",
        "best_evalue",
        "best_bitscore",
    ]

    write_table(candidate_path, candidates, candidate_columns)

    calls = (
        candidates.loc[candidates["credible_match"]].copy()
        if not candidates.empty
        else pd.DataFrame(columns=candidate_columns)
    )

    write_table(calls_path, calls, candidate_columns)

    target_mags = sorted(
        calls["target_MAG"].dropna().astype(str).unique()
    )
    mag_list_path.write_text(
        "\n".join(target_mags)
        + ("\n" if target_mags else ""),
        encoding="utf-8",
    )

    summary = pd.DataFrame(
        [
            {
                "dataset": args.dataset,
                "n_raw_tblastn_HSPs": len(blast),
                "n_filtered_HSPs": len(filtered),
                "n_filtered_hits_without_MAG_mapping": int(
                    unmapped.sum()
                ),
                "n_explicit_query_cassettes": (
                    query_map["cassette_id"].nunique()
                ),
                "n_target_loci": (
                    locus_hits["target_locus_id"].nunique()
                    if not locus_hits.empty
                    else 0
                ),
                "n_credible_cassette_calls": len(calls),
                "n_distinct_query_cassettes_found": (
                    calls["cassette_id"].nunique()
                    if not calls.empty
                    else 0
                ),
                "n_target_MAGs": len(target_mags),
                "n_target_contigs": (
                    calls["target_contig"].nunique()
                    if not calls.empty
                    else 0
                ),
                "min_pident": args.min_pident,
                "min_qcov": args.min_qcov,
                "max_evalue": args.max_evalue,
                "same_locus_overlap": args.same_locus_overlap,
            }
        ]
    )

    summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
    )

    print("Credible EET cassette calling from TBLASTN completed successfully.")
    print(f"Dataset: {args.dataset}")
    print(f"Raw TBLASTN HSPs: {len(blast)}")
    print(f"Filtered HSPs: {len(filtered)}")
    print(f"Explicit query cassettes: {query_map['cassette_id'].nunique()}")
    print(f"Credible cassette calls: {len(calls)}")
    print(f"Target MAGs: {len(target_mags)}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
