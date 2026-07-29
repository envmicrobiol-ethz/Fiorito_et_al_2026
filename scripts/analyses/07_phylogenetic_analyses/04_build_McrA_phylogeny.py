#!/usr/bin/env python3

# DESCRIPTION
# Reconstructs the environmental McrA phylogeny used in this study.
#
# The workflow:
#   1. builds provenance-preserving protein databases from 1,081 dereplicated
#      MAGs and from proteins predicted on >=1-kb contigs from 32 metagenomes;
#   2. searches both databases with the McrA HMM from the GraftM package;
#   3. extracts candidate proteins and applies the reference-calibrated score
#      and HMM-coverage thresholds used for the primary phylogeny;
#   4. collapses exact amino-acid duplicates while retaining all MAG and
#      contig provenance;
#   5. adds the nonredundant environmental sequences to a published aligned
#      McrA/MCR-alpha reference set with MAFFT --addfragments;
#   6. trims the combined alignment with trimAl -gappyout;
#   7. infers a maximum-likelihood tree with IQ-TREE.
#
# INPUT
# 1. A text file containing one MAG protein FASTA path per line.
# 2. A TSV describing the metagenome-contig protein FASTAs, with columns:
#      protein_fasta, raw_sample_name, canonical_sample_name
# 3. The McrA search HMM from the GraftM McrA/K00399 package.
# 4. A pre-aligned McrA/MCR-alpha reference protein FASTA.
#
# OUTPUT
# Search databases, HMMER results, candidate and nonredundant sequence tables,
# combined and trimmed alignments, and the IQ-TREE phylogeny.
#
# SOFTWARE
# Python >=3.11, HMMER 3.4, MAFFT 7.505, trimAl 1.5.rev0 and IQ-TREE 2.2.2.7.
#
# USAGE
# python3 04_build_McrA_phylogeny.py \
#   --mag-protein-list MAG_protein_files.txt \
#   --contig-protein-mapping contig_protein_file_mapping.tsv \
#   --mcrA-hmm graftm_search.hmm \
#   --reference-alignment published_McrA_reference_alignment.faa \
#   --output-dir McrA_phylogeny \
#   --threads 8

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


###############################################################################
# CONSTANTS USED IN THIS STUDY
###############################################################################

DEFAULT_HMM_EVALUE = 1e-5
DEFAULT_MIN_FULL_SCORE = 280.8
DEFAULT_MIN_HMM_COVERAGE = 0.37

DEFAULT_EXPECTED_MAG_FILES = 1081
DEFAULT_EXPECTED_CONTIG_FILES = 32
DEFAULT_EXPECTED_REFERENCE_SEQUENCES = 177
DEFAULT_EXPECTED_PRIMARY_SEQUENCES = 128

SOURCE_NAME_CORRECTIONS = {
    "NorraRomyren_8_325-340cm": "NorraRomyren_8_325-350cm",
    "Lungsmossen_9_525-540cm": "Lungsmossen_9_525-550cm",
}


###############################################################################
# DATA STRUCTURES
###############################################################################

@dataclass(frozen=True)
class Source:
    source_type: str
    protein_fasta: Path
    raw_source_name: str
    canonical_source_name: str


@dataclass(frozen=True)
class ProteinMetadata:
    sequence_id: str
    source_type: str
    protein_fasta: str
    raw_source_name: str
    canonical_source_name: str
    original_protein_id: str


@dataclass(frozen=True)
class HMMHit:
    sequence_id: str
    target_length: int
    query_length: int
    full_evalue: float
    full_score: float
    domain_i_evalue: float
    domain_score: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int

    @property
    def hmm_coverage(self) -> float:
        return (self.hmm_to - self.hmm_from + 1) / self.query_length

    @property
    def aligned_length(self) -> int:
        return self.ali_to - self.ali_from + 1


###############################################################################
# ARGUMENTS
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Identify McrA proteins in MAG and metagenomic-contig protein "
            "databases, deduplicate exact sequences, add them to a reference "
            "alignment and infer the primary McrA phylogeny."
        )
    )

    parser.add_argument(
        "--mag-protein-list",
        required=True,
        type=Path,
        help="Text file containing one MAG protein FASTA path per line.",
    )
    parser.add_argument(
        "--contig-protein-mapping",
        required=True,
        type=Path,
        help=(
            "TSV with protein_fasta, raw_sample_name and "
            "canonical_sample_name columns."
        ),
    )
    parser.add_argument(
        "--mcrA-hmm",
        required=True,
        type=Path,
        help="McrA search HMM from the GraftM McrA/K00399 package.",
    )
    parser.add_argument(
        "--reference-alignment",
        required=True,
        type=Path,
        help="Published aligned McrA/MCR-alpha reference protein FASTA.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads used by HMMER, MAFFT and IQ-TREE. Default: 8.",
    )

    parser.add_argument(
        "--hmm-evalue",
        type=float,
        default=DEFAULT_HMM_EVALUE,
        help="HMMER sequence and domain E-value threshold. Default: 1e-5.",
    )
    parser.add_argument(
        "--min-full-score",
        type=float,
        default=DEFAULT_MIN_FULL_SCORE,
        help=(
            "Reference-calibrated minimum HMMER full-sequence score for the "
            "primary phylogeny. Default: 280.8."
        ),
    )
    parser.add_argument(
        "--min-hmm-coverage",
        type=float,
        default=DEFAULT_MIN_HMM_COVERAGE,
        help=(
            "Reference-calibrated minimum HMM coverage for the primary "
            "phylogeny. Default: 0.37."
        ),
    )

    parser.add_argument(
        "--expected-mag-files",
        type=int,
        default=DEFAULT_EXPECTED_MAG_FILES,
        help="Expected MAG protein FASTAs. Set to 0 to disable. Default: 1081.",
    )
    parser.add_argument(
        "--expected-contig-files",
        type=int,
        default=DEFAULT_EXPECTED_CONTIG_FILES,
        help=(
            "Expected metagenome-contig protein FASTAs. Set to 0 to disable. "
            "Default: 32."
        ),
    )
    parser.add_argument(
        "--expected-reference-sequences",
        type=int,
        default=DEFAULT_EXPECTED_REFERENCE_SEQUENCES,
        help=(
            "Expected reference alignment sequences. Set to 0 to disable. "
            "Default: 177."
        ),
    )
    parser.add_argument(
        "--expected-primary-sequences",
        type=int,
        default=DEFAULT_EXPECTED_PRIMARY_SEQUENCES,
        help=(
            "Expected nonredundant environmental sequences in the primary "
            "phylogeny. Set to 0 to disable. Default: 128."
        ),
    )

    parser.add_argument(
        "--hmmsearch-executable",
        default="hmmsearch",
        help="hmmsearch executable. Default: hmmsearch.",
    )
    parser.add_argument(
        "--mafft-executable",
        default="mafft",
        help="MAFFT executable. Default: mafft.",
    )
    parser.add_argument(
        "--trimal-executable",
        default="trimal",
        help="trimAl executable. Default: trimal.",
    )
    parser.add_argument(
        "--iqtree-executable",
        default=None,
        help=(
            "IQ-TREE executable. Default: automatically select iqtree2 or "
            "iqtree."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing non-empty output directory.",
    )

    return parser.parse_args()


###############################################################################
# GENERAL HELPERS
###############################################################################

def fail(message: str) -> None:
    raise RuntimeError(message)


def require_nonempty_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is missing or empty: {path}")


def canonicalize_name(name: str) -> str:
    corrected = name
    for old, new in SOURCE_NAME_CORRECTIONS.items():
        corrected = corrected.replace(old, new)
    return corrected


def fasta_stem(path: Path) -> str:
    name = path.name
    for suffix in (
        ".faa.gz",
        ".fasta.gz",
        ".fa.gz",
        ".faa",
        ".fasta",
        ".fa",
    ):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    require_nonempty_file(path, "Protein FASTA")

    header: str | None = None
    sequence_parts: list[str] = []

    with open_text(path) as handle:
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


def write_fasta_record(
    handle: TextIO,
    sequence_id: str,
    sequence: str,
    width: int = 80,
) -> None:
    handle.write(f">{sequence_id}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def count_fasta_records(path: Path) -> int:
    return sum(1 for _header, _sequence in read_fasta(path))


def normalize_protein_sequence(
    raw_sequence: str,
    source: Path,
    protein_id: str,
) -> tuple[str, bool]:
    sequence = re.sub(r"\s+", "", raw_sequence).upper()
    terminal_stop_removed = False

    if sequence.endswith("*"):
        sequence = sequence[:-1]
        terminal_stop_removed = True

    if not sequence:
        raise ValueError(f"Empty protein sequence: {source}::{protein_id}")

    if "*" in sequence:
        raise ValueError(
            f"Internal stop symbol in protein: {source}::{protein_id}"
        )

    invalid = re.search(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]", sequence)
    if invalid:
        raise ValueError(
            f"Invalid amino-acid character {invalid.group(0)!r} in "
            f"{source}::{protein_id}"
        )

    return sequence, terminal_stop_removed


def resolve_path(value: str, table_path: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = table_path.parent / path
    return path.resolve()


def run_command(
    command: list[str],
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
) -> None:
    stdout_handle = (
        stdout_path.open("w", encoding="utf-8")
        if stdout_path is not None
        else None
    )
    stderr_handle = (
        stderr_path.open("w", encoding="utf-8")
        if stderr_path is not None
        else None
    )

    try:
        completed = subprocess.run(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            check=False,
        )
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()

    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit status "
            f"{completed.returncode}: {' '.join(command)}"
        )


def require_executable(executable: str) -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise FileNotFoundError(
            f"Required executable was not found in PATH: {executable}"
        )
    return resolved


###############################################################################
# INPUT SOURCES
###############################################################################

def load_mag_sources(
    path: Path,
    expected_count: int,
) -> list[Source]:
    require_nonempty_file(path, "MAG protein file list")

    sources: list[Source] = []
    seen_paths: set[Path] = set()

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            value = raw_line.strip()

            if not value or value.startswith("#"):
                continue

            protein_fasta = resolve_path(value, path)
            require_nonempty_file(
                protein_fasta,
                f"MAG protein FASTA at line {line_number}",
            )

            if protein_fasta in seen_paths:
                raise ValueError(
                    f"Duplicated MAG protein FASTA path: {protein_fasta}"
                )
            seen_paths.add(protein_fasta)

            raw_name = fasta_stem(protein_fasta)
            sources.append(
                Source(
                    source_type="MAG",
                    protein_fasta=protein_fasta,
                    raw_source_name=raw_name,
                    canonical_source_name=canonicalize_name(raw_name),
                )
            )

    if expected_count and len(sources) != expected_count:
        raise ValueError(
            f"Expected {expected_count} MAG protein FASTAs, "
            f"found {len(sources)}."
        )

    return sources


def load_contig_sources(
    path: Path,
    expected_count: int,
) -> list[Source]:
    require_nonempty_file(path, "Contig protein mapping table")

    sources: list[Source] = []
    seen_paths: set[Path] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []

        required = {
            "protein_fasta",
            "raw_sample_name",
            "canonical_sample_name",
        }
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(
                "Contig protein mapping table is missing columns: "
                + ", ".join(sorted(missing))
            )

        for line_number, row in enumerate(reader, start=2):
            protein_fasta = resolve_path(
                row["protein_fasta"].strip(),
                path,
            )
            require_nonempty_file(
                protein_fasta,
                f"Contig protein FASTA at line {line_number}",
            )

            if protein_fasta in seen_paths:
                raise ValueError(
                    f"Duplicated contig protein FASTA path: {protein_fasta}"
                )
            seen_paths.add(protein_fasta)

            raw_name = row["raw_sample_name"].strip()
            canonical_name = row["canonical_sample_name"].strip()

            if not raw_name or not canonical_name:
                raise ValueError(
                    f"Empty sample name at line {line_number} in {path}"
                )

            if canonicalize_name(raw_name) != canonical_name:
                raise ValueError(
                    "Canonical sample-name disagreement at line "
                    f"{line_number}: raw={raw_name!r}, "
                    f"provided={canonical_name!r}, "
                    f"expected={canonicalize_name(raw_name)!r}"
                )

            sources.append(
                Source(
                    source_type="CONTIG",
                    protein_fasta=protein_fasta,
                    raw_source_name=raw_name,
                    canonical_source_name=canonical_name,
                )
            )

    if expected_count and len(sources) != expected_count:
        raise ValueError(
            f"Expected {expected_count} contig protein FASTAs, "
            f"found {len(sources)}."
        )

    return sources


###############################################################################
# DATABASE CONSTRUCTION
###############################################################################

def build_protein_database(
    sources: list[Source],
    output_fasta: Path,
    metadata_output: Path,
) -> dict[str, ProteinMetadata]:
    metadata: dict[str, ProteinMetadata] = {}
    source_rows: list[dict[str, object]] = []

    total_records = 0
    terminal_stops_removed = 0

    with output_fasta.open("w", encoding="utf-8") as output_handle:
        for source in sources:
            source_records = 0
            source_terminal_stops = 0
            seen_original_ids: set[str] = set()

            for original_header, raw_sequence in read_fasta(
                source.protein_fasta
            ):
                original_protein_id = original_header.split()[0]

                if original_protein_id in seen_original_ids:
                    raise ValueError(
                        "Duplicate protein identifier within source: "
                        f"{source.protein_fasta}::{original_protein_id}"
                    )
                seen_original_ids.add(original_protein_id)

                sequence, removed_stop = normalize_protein_sequence(
                    raw_sequence,
                    source.protein_fasta,
                    original_protein_id,
                )

                sequence_id = (
                    f"{source.source_type}__"
                    f"{source.raw_source_name}__"
                    f"{original_protein_id}"
                )

                if sequence_id in metadata:
                    raise ValueError(
                        f"Duplicate combined sequence identifier: {sequence_id}"
                    )

                write_fasta_record(
                    output_handle,
                    sequence_id,
                    sequence,
                )

                metadata[sequence_id] = ProteinMetadata(
                    sequence_id=sequence_id,
                    source_type=source.source_type,
                    protein_fasta=str(source.protein_fasta),
                    raw_source_name=source.raw_source_name,
                    canonical_source_name=source.canonical_source_name,
                    original_protein_id=original_protein_id,
                )

                total_records += 1
                source_records += 1

                if removed_stop:
                    terminal_stops_removed += 1
                    source_terminal_stops += 1

            source_rows.append(
                {
                    "source_type": source.source_type,
                    "protein_fasta": str(source.protein_fasta),
                    "raw_source_name": source.raw_source_name,
                    "canonical_source_name": source.canonical_source_name,
                    "protein_records_written": source_records,
                    "terminal_stops_removed": source_terminal_stops,
                }
            )

    if total_records == 0:
        raise ValueError(f"No protein records were written to {output_fasta}")

    with metadata_output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sequence_id",
                "source_type",
                "protein_fasta",
                "raw_source_name",
                "canonical_source_name",
                "original_protein_id",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for sequence_id in sorted(metadata):
            writer.writerow(metadata[sequence_id].__dict__)

    source_stats_path = metadata_output.with_name(
        metadata_output.stem + "_source_counts.tsv"
    )
    with source_stats_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "source_type",
                "protein_fasta",
                "raw_source_name",
                "canonical_source_name",
                "protein_records_written",
                "terminal_stops_removed",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(source_rows)

    print(
        f"Built {output_fasta.name}: "
        f"{len(sources)} sources, {total_records} proteins, "
        f"{terminal_stops_removed} terminal stops removed."
    )

    return metadata


###############################################################################
# HMMER SEARCH AND PARSING
###############################################################################

def run_hmmsearch(
    executable: str,
    hmm_path: Path,
    database_path: Path,
    output_prefix: Path,
    threads: int,
    evalue: float,
) -> tuple[Path, Path]:
    tblout = output_prefix.with_suffix(".tblout")
    domtblout = output_prefix.with_suffix(".domtblout")
    log = output_prefix.with_suffix(".hmmsearch.log")

    command = [
        executable,
        "--cpu",
        str(threads),
        "--noali",
        "-E",
        str(evalue),
        "--domE",
        str(evalue),
        "--tblout",
        str(tblout),
        "--domtblout",
        str(domtblout),
        str(hmm_path),
        str(database_path),
    ]

    run_command(command, stdout_path=log)

    require_nonempty_file(tblout, "HMMER tblout")
    require_nonempty_file(domtblout, "HMMER domtblout")

    return tblout, domtblout


def parse_domtblout(path: Path) -> dict[str, HMMHit]:
    require_nonempty_file(path, "HMMER domtblout")

    best_hits: dict[str, HMMHit] = {}

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue

            fields = raw_line.rstrip("\n").split(maxsplit=22)

            if len(fields) < 22:
                raise ValueError(
                    f"Malformed domtblout line {line_number} in {path}"
                )

            try:
                hit = HMMHit(
                    sequence_id=fields[0],
                    target_length=int(fields[2]),
                    query_length=int(fields[5]),
                    full_evalue=float(fields[6]),
                    full_score=float(fields[7]),
                    domain_i_evalue=float(fields[12]),
                    domain_score=float(fields[13]),
                    hmm_from=int(fields[15]),
                    hmm_to=int(fields[16]),
                    ali_from=int(fields[17]),
                    ali_to=int(fields[18]),
                )
            except ValueError as exc:
                raise ValueError(
                    f"Invalid numeric value at domtblout line "
                    f"{line_number} in {path}"
                ) from exc

            current = best_hits.get(hit.sequence_id)

            if current is None or (
                hit.domain_score,
                -hit.domain_i_evalue,
                hit.hmm_coverage,
            ) > (
                current.domain_score,
                -current.domain_i_evalue,
                current.hmm_coverage,
            ):
                best_hits[hit.sequence_id] = hit

    return best_hits


def write_candidate_table(
    hits: dict[str, HMMHit],
    metadata: dict[str, ProteinMetadata],
    output_path: Path,
) -> list[dict[str, object]]:
    missing_metadata = sorted(set(hits) - set(metadata))
    if missing_metadata:
        raise ValueError(
            "HMMER hits missing from provenance metadata. Examples: "
            + ", ".join(missing_metadata[:5])
        )

    rows: list[dict[str, object]] = []

    for sequence_id, hit in hits.items():
        meta = metadata[sequence_id]

        rows.append(
            {
                "sequence_id": sequence_id,
                "source_type": meta.source_type,
                "protein_fasta": meta.protein_fasta,
                "raw_source_name": meta.raw_source_name,
                "canonical_source_name": meta.canonical_source_name,
                "original_protein_id": meta.original_protein_id,
                "target_length": hit.target_length,
                "query_length": hit.query_length,
                "full_evalue": hit.full_evalue,
                "full_score": hit.full_score,
                "domain_i_evalue": hit.domain_i_evalue,
                "domain_score": hit.domain_score,
                "hmm_from": hit.hmm_from,
                "hmm_to": hit.hmm_to,
                "hmm_coverage": hit.hmm_coverage,
                "ali_from": hit.ali_from,
                "ali_to": hit.ali_to,
                "aligned_length": hit.aligned_length,
            }
        )

    rows.sort(
        key=lambda row: (
            -float(row["full_score"]),
            -float(row["domain_score"]),
            str(row["sequence_id"]),
        )
    )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        fieldnames = list(rows[0].keys()) if rows else [
            "sequence_id",
            "source_type",
            "protein_fasta",
            "raw_source_name",
            "canonical_source_name",
            "original_protein_id",
            "target_length",
            "query_length",
            "full_evalue",
            "full_score",
            "domain_i_evalue",
            "domain_score",
            "hmm_from",
            "hmm_to",
            "hmm_coverage",
            "ali_from",
            "ali_to",
            "aligned_length",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    return rows


###############################################################################
# EXTRACTION AND EXACT-SEQUENCE DEDUPLICATION
###############################################################################

def extract_selected_sequences(
    fasta_path: Path,
    selected_ids: set[str],
) -> dict[str, str]:
    extracted: dict[str, str] = {}

    for header, raw_sequence in read_fasta(fasta_path):
        sequence_id = header.split()[0]

        if sequence_id not in selected_ids:
            continue

        sequence, _removed_stop = normalize_protein_sequence(
            raw_sequence,
            fasta_path,
            sequence_id,
        )

        if sequence_id in extracted:
            raise ValueError(
                f"Duplicate FASTA identifier in {fasta_path}: {sequence_id}"
            )

        extracted[sequence_id] = sequence

    missing = selected_ids - set(extracted)
    if missing:
        raise ValueError(
            f"Failed to extract {len(missing)} selected sequences from "
            f"{fasta_path}. Examples: {sorted(missing)[:5]}"
        )

    return extracted


def deduplicate_candidates(
    candidate_rows: list[dict[str, object]],
    sequences: dict[str, str],
    output_dir: Path,
    min_full_score: float,
    min_hmm_coverage: float,
) -> tuple[Path, int]:
    all_occurrences_fasta = output_dir / "McrA_all_candidate_occurrences.faa"
    all_occurrences_tsv = output_dir / "McrA_all_candidate_occurrences.tsv"
    unique_fasta = output_dir / "McrA_all_unique_candidate_sequences.faa"
    unique_metadata = output_dir / "McrA_all_unique_candidate_sequences.tsv"
    occurrence_mapping = output_dir / "McrA_occurrence_to_unique_sequence.tsv"
    primary_fasta = output_dir / "McrA_primary_phylogeny_unique_sequences.faa"
    primary_metadata = output_dir / "McrA_primary_phylogeny_unique_sequences.tsv"
    excluded_metadata = (
        output_dir / "McrA_sequences_excluded_from_primary_phylogeny.tsv"
    )

    occurrence_rows: list[dict[str, object]] = []

    for row in candidate_rows:
        sequence_id = str(row["sequence_id"])
        sequence = sequences[sequence_id]

        eligible = (
            float(row["full_score"]) >= min_full_score
            and float(row["hmm_coverage"]) >= min_hmm_coverage
        )

        reasons: list[str] = []
        if float(row["full_score"]) < min_full_score:
            reasons.append("score_below_reference_minimum")
        if float(row["hmm_coverage"]) < min_hmm_coverage:
            reasons.append("hmm_coverage_below_reference_minimum")

        occurrence_rows.append(
            {
                **row,
                "sequence_length": len(sequence),
                "sequence_sha256": hashlib.sha256(
                    sequence.encode("utf-8")
                ).hexdigest(),
                "primary_phylogeny_eligible": "yes" if eligible else "no",
                "exclusion_reason": ";".join(reasons) if reasons else "none",
            }
        )

    occurrence_rows.sort(
        key=lambda row: (
            -float(row["full_score"]),
            str(row["sequence_id"]),
        )
    )

    with all_occurrences_fasta.open("w", encoding="utf-8") as handle:
        for row in occurrence_rows:
            sequence_id = str(row["sequence_id"])
            write_fasta_record(handle, sequence_id, sequences[sequence_id])

    with all_occurrences_tsv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(occurrence_rows[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(occurrence_rows)

    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in occurrence_rows:
        groups[sequences[str(row["sequence_id"])]].append(row)

    sorted_groups = sorted(
        groups.items(),
        key=lambda item: (
            -max(float(row["full_score"]) for row in item[1]),
            -len(item[0]),
            hashlib.sha256(item[0].encode("utf-8")).hexdigest(),
        ),
    )

    query_id_by_sequence = {
        sequence: f"QRY_{index:04d}"
        for index, (sequence, _rows) in enumerate(
            sorted_groups,
            start=1,
        )
    }

    mapping_rows: list[dict[str, object]] = []
    unique_rows: list[dict[str, object]] = []

    with unique_fasta.open("w", encoding="utf-8") as handle:
        for sequence, rows in sorted_groups:
            query_id = query_id_by_sequence[sequence]
            write_fasta_record(handle, query_id, sequence)

            eligible = any(
                row["primary_phylogeny_eligible"] == "yes"
                for row in rows
            )
            best_row = max(
                rows,
                key=lambda row: (
                    float(row["full_score"]),
                    float(row["hmm_coverage"]),
                    int(row["aligned_length"]),
                ),
            )

            mag_names = sorted(
                {
                    str(row["canonical_source_name"])
                    for row in rows
                    if row["source_type"] == "MAG"
                }
            )
            contig_samples = sorted(
                {
                    str(row["canonical_source_name"])
                    for row in rows
                    if row["source_type"] == "CONTIG"
                }
            )
            source_types = sorted(
                {str(row["source_type"]) for row in rows}
            )

            unique_rows.append(
                {
                    "query_id": query_id,
                    "sequence_sha256": hashlib.sha256(
                        sequence.encode("utf-8")
                    ).hexdigest(),
                    "sequence_length": len(sequence),
                    "occurrence_count": len(rows),
                    "source_types": ";".join(source_types),
                    "MAG_occurrences": sum(
                        row["source_type"] == "MAG"
                        for row in rows
                    ),
                    "CONTIG_occurrences": sum(
                        row["source_type"] == "CONTIG"
                        for row in rows
                    ),
                    "distinct_MAGs": len(mag_names),
                    "MAG_names": ";".join(mag_names),
                    "distinct_contig_samples": len(contig_samples),
                    "contig_samples": ";".join(contig_samples),
                    "maximum_full_score": max(
                        float(row["full_score"])
                        for row in rows
                    ),
                    "minimum_full_score": min(
                        float(row["full_score"])
                        for row in rows
                    ),
                    "maximum_hmm_coverage": max(
                        float(row["hmm_coverage"])
                        for row in rows
                    ),
                    "minimum_hmm_coverage": min(
                        float(row["hmm_coverage"])
                        for row in rows
                    ),
                    "best_occurrence_id": best_row["sequence_id"],
                    "best_source_type": best_row["source_type"],
                    "primary_phylogeny_eligible": (
                        "yes" if eligible else "no"
                    ),
                }
            )

            for row in rows:
                mapping_rows.append(
                    {
                        "query_id": query_id,
                        "sequence_id": row["sequence_id"],
                        "source_type": row["source_type"],
                        "raw_source_name": row["raw_source_name"],
                        "canonical_source_name": row[
                            "canonical_source_name"
                        ],
                        "original_protein_id": row[
                            "original_protein_id"
                        ],
                        "full_score": row["full_score"],
                        "hmm_coverage": row["hmm_coverage"],
                        "primary_phylogeny_eligible": row[
                            "primary_phylogeny_eligible"
                        ],
                        "exclusion_reason": row["exclusion_reason"],
                    }
                )

    unique_rows.sort(key=lambda row: str(row["query_id"]))
    mapping_rows.sort(
        key=lambda row: (
            str(row["query_id"]),
            str(row["source_type"]),
            str(row["canonical_source_name"]),
            str(row["sequence_id"]),
        )
    )

    with unique_metadata.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(unique_rows[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(unique_rows)

    with occurrence_mapping.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(mapping_rows[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(mapping_rows)

    primary_rows = [
        row
        for row in unique_rows
        if row["primary_phylogeny_eligible"] == "yes"
    ]
    excluded_rows = [
        row
        for row in unique_rows
        if row["primary_phylogeny_eligible"] == "no"
    ]

    sequence_by_query = {
        query_id_by_sequence[sequence]: sequence
        for sequence in groups
    }

    with primary_fasta.open("w", encoding="utf-8") as handle:
        for row in primary_rows:
            query_id = str(row["query_id"])
            write_fasta_record(
                handle,
                query_id,
                sequence_by_query[query_id],
            )

    for table_path, rows in (
        (primary_metadata, primary_rows),
        (excluded_metadata, excluded_rows),
    ):
        with table_path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(unique_rows[0].keys()),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    return primary_fasta, len(primary_rows)


###############################################################################
# ALIGNMENT AND TREE
###############################################################################

def select_iqtree_executable(user_choice: str | None) -> str:
    if user_choice:
        return require_executable(user_choice)

    for candidate in ("iqtree2", "iqtree"):
        resolved = shutil.which(candidate)
        if resolved is not None:
            return resolved

    raise FileNotFoundError(
        "Neither iqtree2 nor iqtree was found in PATH."
    )


def build_alignment_and_tree(
    primary_fasta: Path,
    reference_alignment: Path,
    output_dir: Path,
    threads: int,
    mafft: str,
    trimal: str,
    iqtree: str,
    expected_reference_count: int,
    expected_primary_count: int,
) -> dict[str, object]:
    alignment_dir = output_dir / "alignment"
    tree_dir = output_dir / "trees"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)

    reference_count = count_fasta_records(reference_alignment)
    primary_count = count_fasta_records(primary_fasta)

    if expected_reference_count and reference_count != expected_reference_count:
        raise ValueError(
            f"Expected {expected_reference_count} reference sequences, "
            f"found {reference_count}."
        )

    if expected_primary_count and primary_count != expected_primary_count:
        raise ValueError(
            f"Expected {expected_primary_count} primary environmental "
            f"sequences, found {primary_count}."
        )

    untrimmed = (
        alignment_dir / "McrA_reference_plus_queries_untrimmed.faa"
    )
    trimmed = (
        alignment_dir / "McrA_reference_plus_queries_trimal_gappyout.faa"
    )
    mafft_log = alignment_dir / "McrA_mafft_addfragments.log"
    trimal_log = alignment_dir / "McrA_trimal_gappyout.log"

    run_command(
        [
            mafft,
            "--auto",
            "--thread",
            str(threads),
            "--addfragments",
            str(primary_fasta),
            str(reference_alignment),
        ],
        stdout_path=untrimmed,
        stderr_path=mafft_log,
    )

    expected_combined = reference_count + primary_count
    combined_count = count_fasta_records(untrimmed)

    if combined_count != expected_combined:
        raise RuntimeError(
            f"Expected {expected_combined} aligned sequences, "
            f"found {combined_count}."
        )

    run_command(
        [
            trimal,
            "-in",
            str(untrimmed),
            "-out",
            str(trimmed),
            "-fasta",
            "-gappyout",
        ],
        stdout_path=trimal_log,
        stderr_path=trimal_log.with_suffix(".stderr.log"),
    )

    trimmed_count = count_fasta_records(trimmed)
    if trimmed_count != expected_combined:
        raise RuntimeError(
            f"Expected {expected_combined} trimmed sequences, "
            f"found {trimmed_count}."
        )

    tree_prefix = tree_dir / "McrA_gappyout"

    run_command(
        [
            iqtree,
            "-s",
            str(trimmed),
            "-st",
            "AA",
            "-m",
            "MFP",
            "-B",
            "1000",
            "--alrt",
            "1000",
            "-bnni",
            "-T",
            str(threads),
            "-seed",
            "12345",
            "--prefix",
            str(tree_prefix),
            "-redo",
        ]
    )

    treefile = Path(str(tree_prefix) + ".treefile")
    iqtree_report = Path(str(tree_prefix) + ".iqtree")

    require_nonempty_file(treefile, "IQ-TREE treefile")
    require_nonempty_file(iqtree_report, "IQ-TREE report")

    return {
        "reference_sequences": reference_count,
        "environmental_sequences": primary_count,
        "combined_sequences": expected_combined,
        "untrimmed_alignment": str(untrimmed),
        "gappyout_alignment": str(trimmed),
        "treefile": str(treefile),
        "iqtree_report": str(iqtree_report),
    }


###############################################################################
# MAIN
###############################################################################

def main() -> None:
    args = parse_args()

    if args.threads < 1:
        raise ValueError("--threads must be at least 1.")
    if not 0 < args.hmm_evalue:
        raise ValueError("--hmm-evalue must be greater than zero.")
    if args.min_full_score < 0:
        raise ValueError("--min-full-score cannot be negative.")
    if not 0 <= args.min_hmm_coverage <= 1:
        raise ValueError("--min-hmm-coverage must be between 0 and 1.")

    mag_list = args.mag_protein_list.resolve()
    contig_mapping = args.contig_protein_mapping.resolve()
    hmm_path = args.mcrA_hmm.resolve()
    reference_alignment = args.reference_alignment.resolve()
    output_dir = args.output_dir.resolve()

    require_nonempty_file(hmm_path, "McrA search HMM")
    require_nonempty_file(reference_alignment, "Reference alignment")

    if output_dir.exists() and any(output_dir.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Output directory exists and is not empty: {output_dir}"
            )
        shutil.rmtree(output_dir)

    database_dir = output_dir / "search_databases"
    hmmsearch_dir = output_dir / "hmmsearch"
    extracted_dir = output_dir / "extracted_sequences"

    database_dir.mkdir(parents=True, exist_ok=True)
    hmmsearch_dir.mkdir(parents=True, exist_ok=True)
    extracted_dir.mkdir(parents=True, exist_ok=True)

    hmmsearch = require_executable(args.hmmsearch_executable)
    mafft = require_executable(args.mafft_executable)
    trimal = require_executable(args.trimal_executable)
    iqtree = select_iqtree_executable(args.iqtree_executable)

    mag_sources = load_mag_sources(
        mag_list,
        args.expected_mag_files,
    )
    contig_sources = load_contig_sources(
        contig_mapping,
        args.expected_contig_files,
    )

    mag_database = database_dir / "MAG_proteins_for_hmmsearch.faa"
    contig_database = (
        database_dir / "contig_proteins_for_hmmsearch.faa"
    )

    mag_metadata = build_protein_database(
        mag_sources,
        mag_database,
        database_dir / "MAG_protein_provenance.tsv",
    )
    contig_metadata = build_protein_database(
        contig_sources,
        contig_database,
        database_dir / "CONTIG_protein_provenance.tsv",
    )

    mag_tblout, mag_domtblout = run_hmmsearch(
        hmmsearch,
        hmm_path,
        mag_database,
        hmmsearch_dir / "MAG_McrA",
        args.threads,
        args.hmm_evalue,
    )
    contig_tblout, contig_domtblout = run_hmmsearch(
        hmmsearch,
        hmm_path,
        contig_database,
        hmmsearch_dir / "CONTIG_McrA",
        args.threads,
        args.hmm_evalue,
    )

    mag_hits = parse_domtblout(mag_domtblout)
    contig_hits = parse_domtblout(contig_domtblout)

    mag_rows = write_candidate_table(
        mag_hits,
        mag_metadata,
        hmmsearch_dir / "MAG_McrA_candidates_detailed.tsv",
    )
    contig_rows = write_candidate_table(
        contig_hits,
        contig_metadata,
        hmmsearch_dir / "CONTIG_McrA_candidates_detailed.tsv",
    )

    candidate_rows = [*mag_rows, *contig_rows]

    if not candidate_rows:
        raise RuntimeError("No McrA candidate proteins were recovered.")

    mag_sequences = extract_selected_sequences(
        mag_database,
        set(mag_hits),
    )
    contig_sequences = extract_selected_sequences(
        contig_database,
        set(contig_hits),
    )
    candidate_sequences = {
        **mag_sequences,
        **contig_sequences,
    }

    primary_fasta, primary_count = deduplicate_candidates(
        candidate_rows,
        candidate_sequences,
        extracted_dir,
        args.min_full_score,
        args.min_hmm_coverage,
    )

    tree_summary = build_alignment_and_tree(
        primary_fasta=primary_fasta,
        reference_alignment=reference_alignment,
        output_dir=output_dir,
        threads=args.threads,
        mafft=mafft,
        trimal=trimal,
        iqtree=iqtree,
        expected_reference_count=args.expected_reference_sequences,
        expected_primary_count=args.expected_primary_sequences,
    )

    summary = {
        "description": "Environmental McrA identification and phylogeny.",
        "inputs": {
            "MAG_protein_list": str(mag_list),
            "contig_protein_mapping": str(contig_mapping),
            "McrA_HMM": str(hmm_path),
            "reference_alignment": str(reference_alignment),
        },
        "software_executables": {
            "hmmsearch": hmmsearch,
            "mafft": mafft,
            "trimal": trimal,
            "iqtree": iqtree,
        },
        "parameters": {
            "threads": args.threads,
            "hmm_sequence_Evalue": args.hmm_evalue,
            "hmm_domain_Evalue": args.hmm_evalue,
            "reference_calibrated_minimum_full_score": (
                args.min_full_score
            ),
            "reference_calibrated_minimum_HMM_coverage": (
                args.min_hmm_coverage
            ),
            "trimAl_mode": "gappyout",
            "IQ_TREE_model": "MFP",
            "ultrafast_bootstrap_replicates": 1000,
            "SH_aLRT_replicates": 1000,
            "IQ_TREE_BNNI": True,
            "IQ_TREE_seed": 12345,
        },
        "counts": {
            "MAG_protein_files": len(mag_sources),
            "contig_protein_files": len(contig_sources),
            "MAG_HMM_candidates": len(mag_rows),
            "CONTIG_HMM_candidates": len(contig_rows),
            "all_candidate_occurrences": len(candidate_rows),
            "primary_nonredundant_environmental_sequences": primary_count,
            **tree_summary,
        },
        "outputs": {
            "MAG_database": str(mag_database),
            "CONTIG_database": str(contig_database),
            "MAG_tblout": str(mag_tblout),
            "CONTIG_tblout": str(contig_tblout),
            "primary_environmental_FASTA": str(primary_fasta),
            "treefile": tree_summary["treefile"],
        },
    }

    summary_path = output_dir / "McrA_phylogeny_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print()
    print("McrA phylogeny completed successfully.")
    print(f"MAG HMM candidates: {len(mag_rows)}")
    print(f"Contig HMM candidates: {len(contig_rows)}")
    print(
        "Primary nonredundant environmental sequences: "
        f"{primary_count}"
    )
    print(f"Tree: {tree_summary['treefile']}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\n[ERROR] {exc}", file=sys.stderr)
        raise
