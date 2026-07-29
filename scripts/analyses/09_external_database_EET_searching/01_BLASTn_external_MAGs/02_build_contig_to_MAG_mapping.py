#!/usr/bin/env python3

# DESCRIPTION
# Builds a mapping between external-dataset contig identifiers and the MAG
# FASTA files from which those contigs originated.
#
# INPUT
# A directory containing one nucleotide FASTA file per MAG.
#
# OUTPUT
# A TSV with columns:
#   target_contig
#   target_MAG
#   target_MAG_path
#
# USAGE
# python3 02_build_contig_to_MAG_mapping.py \
#   --genomes-dir external_MAG_directory \
#   --output contig_to_MAG.tsv

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path
from typing import Iterator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map target contig identifiers to external MAG files."
    )
    parser.add_argument(
        "--genomes-dir",
        required=True,
        type=Path,
        help="Directory containing one MAG nucleotide FASTA per file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output contig-to-MAG TSV.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def normalize_contig_identifier(header_or_id: str) -> str:
    value = header_or_id.strip()

    if value.startswith(">"):
        value = value[1:]

    token = value.split()[0]

    accession = re.search(r"([A-Za-z0-9_]+\.[0-9]+)", token)
    return accession.group(1) if accession else token


def fasta_paths(directory: Path) -> Iterator[Path]:
    patterns = (
        "*.fna",
        "*.fa",
        "*.fasta",
        "*.fna.gz",
        "*.fa.gz",
        "*.fasta.gz",
    )

    seen: set[Path] = set()

    for pattern in patterns:
        for path in sorted(directory.glob(pattern)):
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                yield resolved


def mag_name(path: Path) -> str:
    name = path.name

    for suffix in (
        ".fna.gz",
        ".fa.gz",
        ".fasta.gz",
        ".fna",
        ".fa",
        ".fasta",
    ):
        if name.endswith(suffix):
            return name[: -len(suffix)]

    return path.stem


def main() -> None:
    args = parse_args()

    genomes_dir = args.genomes_dir.resolve()
    output = args.output.resolve()

    if not genomes_dir.is_dir():
        raise FileNotFoundError(
            f"Genome directory does not exist: {genomes_dir}"
        )

    genomes = list(fasta_paths(genomes_dir))

    if not genomes:
        raise ValueError(
            f"No MAG nucleotide FASTA files found in: {genomes_dir}"
        )

    rows: list[dict[str, str]] = []
    contig_to_mag: dict[str, str] = {}

    for genome_path in genomes:
        target_mag = mag_name(genome_path)
        n_contigs = 0

        with open_text(genome_path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                if not raw_line.startswith(">"):
                    continue

                target_contig = normalize_contig_identifier(raw_line)

                if not target_contig:
                    raise ValueError(
                        f"Empty contig identifier in {genome_path}, "
                        f"line {line_number}"
                    )

                previous_mag = contig_to_mag.get(target_contig)

                if previous_mag is not None:
                    raise ValueError(
                        f"Contig identifier {target_contig!r} occurs in "
                        f"multiple MAG files: {previous_mag!r} and "
                        f"{target_mag!r}"
                    )

                contig_to_mag[target_contig] = target_mag

                rows.append(
                    {
                        "target_contig": target_contig,
                        "target_MAG": target_mag,
                        "target_MAG_path": str(genome_path),
                    }
                )
                n_contigs += 1

        if n_contigs == 0:
            raise ValueError(
                f"No FASTA records found in MAG file: {genome_path}"
            )

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "target_contig",
                "target_MAG",
                "target_MAG_path",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print("Contig-to-MAG mapping completed successfully.")
    print(f"MAG files: {len(genomes)}")
    print(f"Contigs mapped: {len(rows)}")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
