#!/usr/bin/env python3

# DESCRIPTION
# Builds an explicit query-gene to EET-cassette map from the translated query
# FASTA. Genes on the same query MAG and contig are sorted by their terminal
# Prodigal gene number and split into consecutive runs.
#
# By default, a new cassette begins when consecutive query gene numbers differ
# by more than one.
#
# INPUT
# 1. Protein query FASTA with identifiers formatted as:
#      query_MAG|queryGene
#    where queryGene ends in a numeric Prodigal suffix.
# 2. A text file containing one essential EET query-gene identifier per line.
#
# OUTPUT
# A query-gene-to-cassette TSV and a cassette summary TSV.
#
# USAGE
# python3 03_build_explicit_EET_query_cassette_map.py \
#   --query-fasta EET_query_proteins.faa \
#   --essential-gene-ids EET_essential_gene_ids.txt \
#   --output-map EET_query_gene_to_cassette.tsv \
#   --summary-tsv EET_query_cassettes_summary.tsv \
#   --max-gene-gap 1

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build an explicit query-gene to EET-cassette map from "
            "consecutive Prodigal gene numbers."
        )
    )
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument("--essential-gene-ids", required=True, type=Path)
    parser.add_argument("--output-map", required=True, type=Path)
    parser.add_argument("--summary-tsv", required=True, type=Path)
    parser.add_argument(
        "--max-gene-gap",
        type=int,
        default=1,
        help=(
            "Maximum difference between consecutive Prodigal gene numbers "
            "within one query cassette. Default: 1."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing outputs.",
    )
    return parser.parse_args()


def require_nonempty(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is missing or empty: {path}")


def fasta_ids(path: Path) -> Iterator[str]:
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith(">"):
                yield raw_line[1:].strip().split()[0]


def parse_query_id(
    query_id: str,
) -> tuple[str, str, str, int]:
    if "|" not in query_id:
        raise ValueError(
            f"Query ID does not contain query_MAG|queryGene: {query_id}"
        )

    query_mag, query_gene = query_id.split("|", 1)

    if not query_mag or not query_gene:
        raise ValueError(f"Malformed query identifier: {query_id}")

    match = re.fullmatch(r"(.+)_([0-9]+)", query_gene)

    if match is None:
        raise ValueError(
            f"Query gene does not end in a numeric Prodigal suffix: "
            f"{query_gene}"
        )

    return (
        query_mag,
        query_gene,
        match.group(1),
        int(match.group(2)),
    )


def read_essential_ids(path: Path) -> set[str]:
    require_nonempty(path, "Essential EET gene list")

    identifiers: set[str] = set()

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            token = raw_line.strip()

            if not token or token.startswith("#"):
                continue

            token = token.split()[0].lstrip(">")
            identifiers.add(token)

            if "|" in token:
                identifiers.add(token.split("|", 1)[1])

    if not identifiers:
        raise ValueError("Essential EET gene list contains no identifiers.")

    return identifiers


def main() -> None:
    args = parse_args()

    if args.max_gene_gap < 1:
        raise ValueError("--max-gene-gap must be at least 1.")

    query_fasta = args.query_fasta.resolve()
    essential_path = args.essential_gene_ids.resolve()
    output_map = args.output_map.resolve()
    summary_tsv = args.summary_tsv.resolve()

    require_nonempty(query_fasta, "Protein query FASTA")
    essential = read_essential_ids(essential_path)

    existing = [
        path
        for path in (output_map, summary_tsv)
        if path.exists()
    ]

    if existing and not args.overwrite:
        raise FileExistsError(
            "Output already exists. Use --overwrite to replace it: "
            + ", ".join(str(path) for path in existing)
        )

    grouped: dict[
        tuple[str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    seen_ids: set[str] = set()

    for query_id in fasta_ids(query_fasta):
        if query_id in seen_ids:
            raise ValueError(f"Duplicate query FASTA ID: {query_id}")
        seen_ids.add(query_id)

        (
            query_mag,
            query_gene,
            query_contig,
            gene_number,
        ) = parse_query_id(query_id)

        grouped[(query_mag, query_contig)].append(
            {
                "qseqid": query_id,
                "query_MAG": query_mag,
                "query_gene": query_gene,
                "query_contig": query_contig,
                "gene_number": gene_number,
                "is_essential": (
                    query_id in essential
                    or query_gene in essential
                ),
            }
        )

    if not seen_ids:
        raise ValueError("No query FASTA identifiers were found.")

    mapping_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    def save_run(
        query_mag: str,
        query_contig: str,
        rows: list[dict[str, object]],
        run_number: int,
    ) -> None:
        if not rows:
            return

        cassette_id = (
            f"{query_mag}|{query_contig}|cassette_{run_number:03d}"
        )

        gene_numbers = [
            int(row["gene_number"])
            for row in rows
        ]

        for row in rows:
            mapping_rows.append(
                {
                    **row,
                    "cassette_id": cassette_id,
                }
            )

        summary_rows.append(
            {
                "cassette_id": cassette_id,
                "query_MAG": query_mag,
                "query_contig": query_contig,
                "first_gene_number": min(gene_numbers),
                "last_gene_number": max(gene_numbers),
                "n_query_genes": len(rows),
                "n_essential_genes": sum(
                    bool(row["is_essential"])
                    for row in rows
                ),
                "query_genes": ",".join(
                    str(row["query_gene"])
                    for row in rows
                ),
            }
        )

    for (query_mag, query_contig), genes in sorted(grouped.items()):
        genes = sorted(
            genes,
            key=lambda row: int(row["gene_number"]),
        )

        gene_numbers = [
            int(row["gene_number"])
            for row in genes
        ]

        if len(gene_numbers) != len(set(gene_numbers)):
            raise ValueError(
                f"Duplicate Prodigal gene numbers on "
                f"{query_mag}|{query_contig}: {gene_numbers}"
            )

        run_number = 1
        previous_gene_number: int | None = None
        run_rows: list[dict[str, object]] = []

        for row in genes:
            current_gene_number = int(row["gene_number"])

            if (
                previous_gene_number is not None
                and current_gene_number - previous_gene_number
                > args.max_gene_gap
            ):
                save_run(
                    query_mag,
                    query_contig,
                    run_rows,
                    run_number,
                )
                run_number += 1
                run_rows = []

            run_rows.append(row)
            previous_gene_number = current_gene_number

        save_run(
            query_mag,
            query_contig,
            run_rows,
            run_number,
        )

    mapping_rows.sort(
        key=lambda row: (
            str(row["query_MAG"]),
            str(row["query_contig"]),
            int(row["gene_number"]),
        )
    )
    summary_rows.sort(
        key=lambda row: (
            str(row["query_MAG"]),
            str(row["query_contig"]),
            int(row["first_gene_number"]),
        )
    )

    output_map.parent.mkdir(parents=True, exist_ok=True)
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    with output_map.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "qseqid",
                "query_MAG",
                "query_gene",
                "query_contig",
                "gene_number",
                "cassette_id",
                "is_essential",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(mapping_rows)

    with summary_tsv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cassette_id",
                "query_MAG",
                "query_contig",
                "first_gene_number",
                "last_gene_number",
                "n_query_genes",
                "n_essential_genes",
                "query_genes",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    multi_cassette_contigs = sum(
        count > 1
        for count in (
            sum(
                row["query_MAG"] == query_mag
                and row["query_contig"] == query_contig
                for row in summary_rows
            )
            for query_mag, query_contig in grouped
        )
    )

    print("Explicit EET query-cassette mapping completed successfully.")
    print(f"Query proteins: {len(mapping_rows)}")
    print(f"Query MAG-contig combinations: {len(grouped)}")
    print(f"Reconstructed query cassettes: {len(summary_rows)}")
    print(
        "Query contigs containing more than one cassette: "
        f"{multi_cassette_contigs}"
    )
    print(
        "Cassettes containing at least one essential gene: "
        f"{sum(int(row['n_essential_genes']) >= 1 for row in summary_rows)}"
    )
    print(f"Mapping: {output_map}")
    print(f"Summary: {summary_tsv}")


if __name__ == "__main__":
    main()
