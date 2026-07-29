#!/usr/bin/env python3

# DESCRIPTION
# Translates nucleotide CDS sequences from conserved EET gene cassettes into
# protein sequences for TBLASTN searches. Translation follows bacterial
# genetic code 11, including recognition of alternative bacterial start
# codons as methionine.
#
# INPUT
# A nucleotide CDS FASTA. Query identifiers should be retained in the form:
#   query_MAG|queryGene
#
# OUTPUT
# A protein FASTA with unchanged query identifiers and a translation report.
#
# USAGE
# python3 01_translate_EET_CDS_to_proteins.py \
#   --input-fasta EET_query_genes.fna \
#   --output-fasta EET_query_proteins.faa \
#   --report-tsv EET_translation_report.tsv

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterator


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODONS = {
    "ATG", "GTG", "TTG", "CTG",
    "ATT", "ATC", "ATA",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Translate nucleotide EET cassette CDS sequences into proteins "
            "for TBLASTN searches."
        )
    )
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--report-tsv", required=True, type=Path)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output files.",
    )
    return parser.parse_args()


def require_nonempty(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is missing or empty: {path}")


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    sequence_parts: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:].strip()
                sequence_parts = []

                if not header:
                    raise ValueError(
                        f"Empty FASTA header at line {line_number} in {path}"
                    )
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence before first FASTA header at line "
                        f"{line_number} in {path}"
                    )
                sequence_parts.append(line)

    if header is not None:
        yield header, "".join(sequence_parts)


def translate_cds(
    nucleotide_sequence: str,
) -> tuple[str, int, int, bool]:
    sequence = re.sub(
        r"[\s-]+",
        "",
        nucleotide_sequence.upper().replace("U", "T"),
    )

    invalid = re.search(r"[^ACGTRYSWKMBDHVN]", sequence)
    if invalid:
        raise ValueError(
            f"Invalid nucleotide character: {invalid.group(0)!r}"
        )

    remainder = len(sequence) % 3
    trimmed_bases = remainder

    if remainder:
        sequence = sequence[:-remainder]

    amino_acids: list[str] = []

    for position in range(0, len(sequence), 3):
        codon = sequence[position : position + 3]
        amino_acids.append(CODON_TABLE.get(codon, "X"))

    if amino_acids and sequence[:3] in START_CODONS:
        amino_acids[0] = "M"

    protein = "".join(amino_acids)

    terminal_stop = protein.endswith("*")
    if terminal_stop:
        protein = protein[:-1]

    internal_stops = protein.count("*")
    protein = protein.replace("*", "X")

    return protein, trimmed_bases, internal_stops, terminal_stop


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(
        sequence[start : start + width]
        for start in range(0, len(sequence), width)
    )


def main() -> None:
    args = parse_args()

    input_fasta = args.input_fasta.resolve()
    output_fasta = args.output_fasta.resolve()
    report_tsv = args.report_tsv.resolve()

    require_nonempty(input_fasta, "Input nucleotide CDS FASTA")

    existing = [
        path
        for path in (output_fasta, report_tsv)
        if path.exists()
    ]

    if existing and not args.overwrite:
        raise FileExistsError(
            "Output already exists. Use --overwrite to replace it: "
            + ", ".join(str(path) for path in existing)
        )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    report_rows: list[dict[str, object]] = []
    seen_ids: set[str] = set()

    input_count = 0
    written_count = 0
    trimmed_count = 0
    internal_stop_count = 0
    empty_count = 0

    with output_fasta.open("w", encoding="utf-8") as output_handle:
        for full_header, nucleotide_sequence in fasta_records(input_fasta):
            input_count += 1
            sequence_id = full_header.split()[0]

            if sequence_id in seen_ids:
                raise ValueError(
                    f"Duplicate FASTA identifier: {sequence_id}"
                )
            seen_ids.add(sequence_id)

            (
                protein,
                trimmed_bases,
                internal_stops,
                terminal_stop,
            ) = translate_cds(nucleotide_sequence)

            if trimmed_bases:
                trimmed_count += 1

            if internal_stops:
                internal_stop_count += 1

            if protein:
                output_handle.write(f">{full_header}\n")
                output_handle.write(wrap(protein))
                output_handle.write("\n")
                written_count += 1
                status = "written"
            else:
                empty_count += 1
                status = "empty_translation"

            report_rows.append(
                {
                    "sequence_id": sequence_id,
                    "nucleotide_length": len(
                        re.sub(r"[\s-]+", "", nucleotide_sequence)
                    ),
                    "protein_length": len(protein),
                    "trimmed_terminal_bases": trimmed_bases,
                    "internal_stop_codons": internal_stops,
                    "terminal_stop_removed": int(terminal_stop),
                    "status": status,
                }
            )

    with report_tsv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as report_handle:
        writer = csv.DictWriter(
            report_handle,
            fieldnames=[
                "sequence_id",
                "nucleotide_length",
                "protein_length",
                "trimmed_terminal_bases",
                "internal_stop_codons",
                "terminal_stop_removed",
                "status",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(report_rows)

    if input_count == 0:
        raise ValueError("No nucleotide CDS records were found.")

    if empty_count:
        raise RuntimeError(
            f"{empty_count} CDS records produced empty translations. "
            f"Inspect: {report_tsv}"
        )

    print("EET CDS translation completed successfully.")
    print(f"Input CDS sequences: {input_count}")
    print(f"Protein sequences written: {written_count}")
    print(f"Sequences with terminal bases trimmed: {trimmed_count}")
    print(f"Sequences with internal stop codons replaced by X: {internal_stop_count}")
    print(f"Protein FASTA: {output_fasta}")
    print(f"Translation report: {report_tsv}")


if __name__ == "__main__":
    main()
