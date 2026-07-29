#!/usr/bin/env python3

# DESCRIPTION
# Identifies putative multiheme cytochromes from predicted protein FASTA files.
#
# INPUT
# A directory containing protein FASTA files generated from MAGs, reference
# genomes or metagenomic assemblies, for example using Prodigal.
#
# OUTPUT
# Tables and FASTA files containing all motif-bearing proteins and proteins
# with at least three motifs of the same type.
#
# USAGE
# python3 01_identify_putative_MHCs.py \
#   predicted_proteins_directory \
#   output_directory
#
# Optional recursive search:
# python3 01_identify_putative_MHCs.py \
#   predicted_proteins_directory \
#   output_directory \
#   --recursive

import argparse
import csv
import gzip
import re
from pathlib import Path
from typing import Iterator, TextIO


###############################################################################
# MOTIF DEFINITIONS
###############################################################################

# Matches are non-overlapping within each motif category.
# Motif categories are not mutually exclusive.
MOTIF_PATTERNS = {
    "CxxCH": re.compile(r"C[A-Z]{2}CH"),
    "CxxCnH": re.compile(r"C[A-Z]{2}C[A-Z]{0,70}?H"),
    "CXnCH": re.compile(r"C[A-Z]{1,70}?CH"),
}

DEFAULT_EXTENSIONS = (
    ".faa",
    ".faa.gz",
    ".fasta",
    ".fasta.gz",
)

OUTPUT_FIELDS = [
    "source_id",
    "source_file",
    "protein_id",
    "contig_id",
    "combined_id",
    "sequence_length",
    "num_CxxCH",
    "num_CxxCnH",
    "num_CXnCH",
    "motif_types_present",
    "qualifying_motif_types",
    "CxxCH_motifs",
    "CxxCnH_motifs",
    "CXnCH_motifs",
]


###############################################################################
# ARGUMENTS
###############################################################################

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Identify predicted proteins containing heme-binding motifs and "
            "retain candidate MHCs containing multiple motifs of the same type."
        )
    )

    parser.add_argument(
        "input_directory",
        type=Path,
        help="Directory containing predicted protein FASTA files.",
    )

    parser.add_argument(
        "output_directory",
        type=Path,
        help="Directory in which output files will be written.",
    )

    parser.add_argument(
        "--min-motifs",
        type=int,
        default=3,
        help=(
            "Minimum number of motifs of the same type required to classify "
            "a protein as a candidate MHC. Default: 3."
        ),
    )

    parser.add_argument(
        "--extensions",
        nargs="+",
        default=list(DEFAULT_EXTENSIONS),
        help=(
            "Accepted protein FASTA extensions. "
            "Default: .faa .faa.gz .fasta .fasta.gz"
        ),
    )

    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search recursively within the input directory.",
    )

    parser.add_argument(
        "--expected-files",
        type=int,
        default=None,
        help="Optional expected number of protein FASTA files.",
    )

    return parser.parse_args()


###############################################################################
# GENERAL FUNCTIONS
###############################################################################

def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def normalize_extensions(extensions: list[str]) -> tuple[str, ...]:
    normalized = []

    for extension in extensions:
        extension = extension.strip()

        if not extension:
            continue

        if not extension.startswith("."):
            extension = f".{extension}"

        normalized.append(extension)

    if not normalized:
        die("no valid protein FASTA extensions were provided")

    # Longest first so that .faa.gz is removed before .gz.
    return tuple(sorted(set(normalized), key=len, reverse=True))


def discover_protein_fastas(
    input_directory: Path,
    extensions: tuple[str, ...],
    recursive: bool,
) -> list[Path]:

    iterator = (
        input_directory.rglob("*")
        if recursive
        else input_directory.glob("*")
    )

    fasta_files = sorted(
        path
        for path in iterator
        if path.is_file()
        and any(path.name.endswith(extension) for extension in extensions)
    )

    return fasta_files


def remove_known_extension(
    filename: str,
    extensions: tuple[str, ...],
) -> str:

    for extension in extensions:
        if filename.endswith(extension):
            return filename[: -len(extension)]

    return filename


def make_source_id(
    fasta_path: Path,
    input_directory: Path,
    extensions: tuple[str, ...],
) -> str:

    relative_path = fasta_path.relative_to(input_directory)
    relative_text = str(relative_path)

    source_id = remove_known_extension(relative_text, extensions)

    # Preserve nested-directory information while generating safe identifiers.
    source_id = source_id.replace("/", "__")
    source_id = source_id.replace("\\", "__")
    source_id = source_id.replace(" ", "_")

    if not source_id:
        die(f"could not generate source identifier for {fasta_path}")

    return source_id


def open_text_file(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")

    return path.open("r")


###############################################################################
# FASTA FUNCTIONS
###############################################################################

def normalize_sequence(
    raw_sequence: str,
    protein_id: str,
    fasta_path: Path,
) -> str:

    sequence = raw_sequence.upper()

    invalid_character = re.search(
        r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO*]",
        sequence,
    )

    if invalid_character:
        die(
            f"invalid character {invalid_character.group(0)!r} in "
            f"{protein_id} from {fasta_path}"
        )

    if sequence.endswith("*"):
        sequence = sequence[:-1]

    if "*" in sequence:
        die(
            f"internal stop codon in {protein_id} from {fasta_path}"
        )

    if not sequence:
        die(
            f"empty protein sequence for {protein_id} in {fasta_path}"
        )

    return sequence


def read_fasta(
    fasta_path: Path,
) -> Iterator[tuple[str, str]]:

    header = None
    sequence_parts: list[str] = []
    records_seen = 0

    with open_text_file(fasta_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    protein_id = header.split()[0]

                    yield (
                        header,
                        normalize_sequence(
                            "".join(sequence_parts),
                            protein_id,
                            fasta_path,
                        ),
                    )

                    records_seen += 1

                header = line[1:].strip()
                sequence_parts = []

                if not header:
                    die(
                        f"empty FASTA header at line {line_number} "
                        f"in {fasta_path}"
                    )

            else:
                if header is None:
                    die(
                        f"sequence before the first FASTA header at "
                        f"line {line_number} in {fasta_path}"
                    )

                sequence_parts.append(line)

    if header is not None:
        protein_id = header.split()[0]

        yield (
            header,
            normalize_sequence(
                "".join(sequence_parts),
                protein_id,
                fasta_path,
            ),
        )

        records_seen += 1

    if records_seen == 0:
        die(f"no FASTA records found in {fasta_path}")


def infer_contig_id(protein_id: str) -> str:
    """
    Infer the originating contig from a standard Prodigal protein identifier.

    For example:
        contig_15_4 -> contig_15

    If the identifier does not end in an underscore followed by an integer,
    the complete protein identifier is retained.
    """

    match = re.fullmatch(r"(.+)_([0-9]+)", protein_id)

    if match:
        return match.group(1)

    return protein_id


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(
        sequence[start : start + width]
        for start in range(0, len(sequence), width)
    )


###############################################################################
# OUTPUT FUNCTIONS
###############################################################################

def write_table(
    output_path: Path,
    rows: list[dict[str, object]],
) -> None:

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


def write_fasta(
    output_path: Path,
    records: list[tuple[str, str, dict[str, object]]],
) -> None:

    with output_path.open("w") as handle:
        for combined_id, sequence, row in records:
            handle.write(
                f">{combined_id} "
                f"contig_id={row['contig_id']} "
                f"num_CxxCH={row['num_CxxCH']} "
                f"num_CxxCnH={row['num_CxxCnH']} "
                f"num_CXnCH={row['num_CXnCH']}\n"
            )

            handle.write(wrap_sequence(sequence))
            handle.write("\n")


###############################################################################
# MAIN
###############################################################################

def main() -> None:
    args = parse_arguments()

    input_directory = args.input_directory.resolve()
    output_directory = args.output_directory.resolve()

    if not input_directory.is_dir():
        die(f"input directory does not exist: {input_directory}")

    if args.min_motifs < 1:
        die("--min-motifs must be at least 1")

    if args.expected_files is not None and args.expected_files < 1:
        die("--expected-files must be a positive integer")

    extensions = normalize_extensions(args.extensions)

    protein_fastas = discover_protein_fastas(
        input_directory,
        extensions,
        args.recursive,
    )

    if not protein_fastas:
        die(
            f"no protein FASTA files with extensions "
            f"{', '.join(extensions)} were found in {input_directory}"
        )

    if (
        args.expected_files is not None
        and len(protein_fastas) != args.expected_files
    ):
        die(
            f"expected {args.expected_files} protein FASTA files, "
            f"but found {len(protein_fastas)}"
        )

    if output_directory.exists() and any(output_directory.iterdir()):
        die(
            f"output directory already exists and is not empty: "
            f"{output_directory}"
        )

    output_directory.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []

    all_records: list[tuple[str, str, dict[str, object]]] = []
    selected_records: list[tuple[str, str, dict[str, object]]] = []

    summary_rows = []
    seen_combined_ids = set()

    total_proteins = 0

    for file_number, fasta_path in enumerate(
        protein_fastas,
        start=1,
    ):
        source_id = make_source_id(
            fasta_path,
            input_directory,
            extensions,
        )

        proteins_scanned = 0
        proteins_with_motif = 0
        candidate_mhcs = 0
        seen_protein_ids = set()

        for header, sequence in read_fasta(fasta_path):
            protein_id = header.split()[0]

            if protein_id in seen_protein_ids:
                die(
                    f"duplicated protein ID {protein_id!r} "
                    f"in {fasta_path}"
                )

            seen_protein_ids.add(protein_id)

            combined_id = f"{source_id}|{protein_id}"

            if combined_id in seen_combined_ids:
                die(
                    f"duplicated combined protein ID: {combined_id}"
                )

            seen_combined_ids.add(combined_id)
            proteins_scanned += 1

            motif_hits = {
                motif_type: [
                    match.group(0)
                    for match in pattern.finditer(sequence)
                ]
                for motif_type, pattern in MOTIF_PATTERNS.items()
            }

            motif_counts = {
                motif_type: len(hits)
                for motif_type, hits in motif_hits.items()
            }

            present_types = [
                motif_type
                for motif_type, count in motif_counts.items()
                if count > 0
            ]

            if not present_types:
                continue

            qualifying_types = [
                motif_type
                for motif_type, count in motif_counts.items()
                if count >= args.min_motifs
            ]

            contig_id = infer_contig_id(protein_id)

            row = {
                "source_id": source_id,
                "source_file": str(
                    fasta_path.relative_to(input_directory)
                ),
                "protein_id": protein_id,
                "contig_id": contig_id,
                "combined_id": combined_id,
                "sequence_length": len(sequence),
                "num_CxxCH": motif_counts["CxxCH"],
                "num_CxxCnH": motif_counts["CxxCnH"],
                "num_CXnCH": motif_counts["CXnCH"],
                "motif_types_present": ";".join(present_types),
                "qualifying_motif_types": ";".join(
                    qualifying_types
                ),
                "CxxCH_motifs": ";".join(
                    motif_hits["CxxCH"]
                ),
                "CxxCnH_motifs": ";".join(
                    motif_hits["CxxCnH"]
                ),
                "CXnCH_motifs": ";".join(
                    motif_hits["CXnCH"]
                ),
            }

            all_rows.append(row)
            all_records.append(
                (combined_id, sequence, row)
            )

            proteins_with_motif += 1

            if qualifying_types:
                selected_rows.append(row)
                selected_records.append(
                    (combined_id, sequence, row)
                )

                candidate_mhcs += 1

        total_proteins += proteins_scanned

        summary_rows.append(
            {
                "source_id": source_id,
                "source_file": str(
                    fasta_path.relative_to(input_directory)
                ),
                "proteins_scanned": proteins_scanned,
                "proteins_with_any_target_motif": (
                    proteins_with_motif
                ),
                "candidate_MHCs": candidate_mhcs,
            }
        )

        print(
            f"[{file_number}/{len(protein_fastas)}] "
            f"{source_id}: "
            f"{proteins_scanned} proteins, "
            f"{candidate_mhcs} candidate MHCs"
        )

    all_tsv = output_directory / "MHC_all_motif_candidates.tsv"
    all_faa = output_directory / "MHC_all_motif_candidates.faa"

    selected_tsv = (
        output_directory
        / f"MHC_ge{args.min_motifs}_same_motif_type.tsv"
    )

    selected_faa = (
        output_directory
        / f"MHC_ge{args.min_motifs}_same_motif_type.faa"
    )

    summary_tsv = output_directory / "MHC_summary_per_source.tsv"
    completion_file = output_directory / "MHC_search.done"

    write_table(all_tsv, all_rows)
    write_fasta(all_faa, all_records)

    write_table(selected_tsv, selected_rows)
    write_fasta(selected_faa, selected_records)

    with summary_tsv.open("w", newline="") as handle:
        fieldnames = [
            "source_id",
            "source_file",
            "proteins_scanned",
            "proteins_with_any_target_motif",
            "candidate_MHCs",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(summary_rows)

    with completion_file.open("w") as handle:
        handle.write(
            f"protein_FASTA_files\t{len(protein_fastas)}\n"
        )
        handle.write(
            f"proteins_scanned\t{total_proteins}\n"
        )
        handle.write(
            f"proteins_with_any_target_motif\t{len(all_rows)}\n"
        )
        handle.write(
            f"candidate_MHCs\t{len(selected_rows)}\n"
        )
        handle.write(
            f"minimum_motifs_same_type\t{args.min_motifs}\n"
        )

    print()
    print("MHC motif search completed successfully.")
    print(f"Protein FASTA files: {len(protein_fastas)}")
    print(f"Proteins scanned: {total_proteins}")
    print(
        f"Proteins with at least one target motif: "
        f"{len(all_rows)}"
    )
    print(
        f"Proteins with >= {args.min_motifs} motifs "
        f"of at least one same type: {len(selected_rows)}"
    )
    print(f"Candidate MHC table: {selected_tsv}")
    print(f"Candidate MHC FASTA: {selected_faa}")
    print(f"Per-source summary: {summary_tsv}")


if __name__ == "__main__":
    main()
