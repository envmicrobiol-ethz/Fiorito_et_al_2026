#!/usr/bin/env python3

# DESCRIPTION
# Extracts genomic regions spanning a target MHC gene and up to five
# neighboring CDS features on each side from MAGs or reference genomes.
#
# INPUT
# 1. A genome manifest TSV with columns: genome, fasta, gff
# 2. A target TSV with columns: genome, protein_id
#    An optional gff_id column may be supplied when the identifier used in
#    the GFF differs from the protein identifier.
#
# OUTPUT
# One nucleotide FASTA per target, a combined context FASTA, a context
# manifest, a gene-membership table and a JSON summary.
#
# USAGE
# python3 01_extract_MHC_genomic_contexts.py \
#   --genome-manifest genome_manifest.tsv \
#   --targets selected_MHCs.tsv \
#   --output-dir MHC_genomic_contexts \
#   --flank-genes 5

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


@dataclass(frozen=True)
class CDS:
    contig: str
    start: int
    end: int
    strand: str
    phase: str
    attributes_raw: str
    attributes: dict[str, str]

    @property
    def left(self) -> int:
        return min(self.start, self.end)

    @property
    def right(self) -> int:
        return max(self.start, self.end)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract nucleotide regions containing a target MHC and a defined "
            "number of neighboring genes on the same contig."
        )
    )
    parser.add_argument(
        "--genome-manifest",
        required=True,
        type=Path,
        help="TSV with columns: genome, fasta, gff.",
    )
    parser.add_argument(
        "--targets",
        required=True,
        type=Path,
        help="TSV with columns: genome, protein_id; optional: gff_id.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--flank-genes",
        type=int,
        default=5,
        help="Maximum number of CDS features selected on each side. Default: 5.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing non-empty output directory.",
    )
    return parser.parse_args()


def require_nonempty(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is missing or empty: {path}")


def resolve_from_table(value: str, table_path: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = table_path.parent / path
    return path.resolve()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require_nonempty(path, "TSV input")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def read_genome_manifest(path: Path) -> dict[str, tuple[Path, Path]]:
    fieldnames, rows = read_tsv(path)
    required = {"genome", "fasta", "gff"}
    missing = required - set(fieldnames)
    if missing:
        raise ValueError(
            f"Genome manifest is missing columns: {', '.join(sorted(missing))}"
        )

    genomes: dict[str, tuple[Path, Path]] = {}
    for row_number, row in enumerate(rows, start=2):
        genome = row["genome"].strip()
        if not genome:
            raise ValueError(f"Empty genome identifier at line {row_number}")
        if genome in genomes:
            raise ValueError(f"Duplicated genome identifier: {genome}")

        fasta = resolve_from_table(row["fasta"].strip(), path)
        gff = resolve_from_table(row["gff"].strip(), path)
        require_nonempty(fasta, f"Genome FASTA for {genome}")
        require_nonempty(gff, f"GFF for {genome}")
        genomes[genome] = (fasta, gff)

    if not genomes:
        raise ValueError("Genome manifest contains no genomes.")

    return genomes


def read_targets(path: Path) -> list[dict[str, str]]:
    fieldnames, rows = read_tsv(path)
    required = {"genome", "protein_id"}
    missing = required - set(fieldnames)
    if missing:
        raise ValueError(
            f"Target table is missing columns: {', '.join(sorted(missing))}"
        )

    targets: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for row_number, row in enumerate(rows, start=2):
        genome = row["genome"].strip()
        protein_id = row["protein_id"].strip()
        gff_id = row.get("gff_id", "").strip()

        if not genome or not protein_id:
            raise ValueError(
                f"Empty genome or protein_id in target table at line {row_number}"
            )

        key = (genome, protein_id)
        if key in seen:
            raise ValueError(
                f"Duplicated target in target table: {genome}, {protein_id}"
            )
        seen.add(key)

        targets.append(
            {
                "genome": genome,
                "protein_id": protein_id,
                "gff_id": gff_id,
            }
        )

    if not targets:
        raise ValueError("Target table contains no targets.")

    return targets


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    header: str | None = None
    parts: list[str] = []

    def finish() -> None:
        nonlocal header, parts
        if header is None:
            return
        identifier = header.split()[0]
        if identifier in sequences:
            raise ValueError(f"Duplicated FASTA identifier in {path}: {identifier}")
        sequence = re.sub(r"\s+", "", "".join(parts)).upper()
        if not sequence:
            raise ValueError(f"Empty FASTA sequence in {path}: {identifier}")
        invalid = re.search(r"[^ACGTURYSWKMBDHVN]", sequence)
        if invalid:
            raise ValueError(
                f"Invalid nucleotide {invalid.group(0)!r} in {path}: {identifier}"
            )
        sequences[identifier] = sequence.replace("U", "T")

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                parts = []
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
                parts.append(line)
        finish()

    if not sequences:
        raise ValueError(f"No FASTA records found in {path}")

    return sequences


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[unquote(key.strip())] = unquote(value.strip())
    return attributes


def read_gff(path: Path) -> list[CDS]:
    features: list[CDS] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Expected 9 GFF columns at line {line_number} in {path}"
                )

            contig, _, feature_type, start, end, _, strand, phase, attrs = fields
            if feature_type != "CDS":
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid coordinates at line {line_number} in {path}"
                ) from exc

            if start_i < 1 or end_i < 1:
                raise ValueError(
                    f"Non-positive coordinates at line {line_number} in {path}"
                )

            if strand not in {"+", "-", ".", "?"}:
                raise ValueError(
                    f"Invalid strand at line {line_number} in {path}: {strand}"
                )

            features.append(
                CDS(
                    contig=contig,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    phase=phase,
                    attributes_raw=attrs,
                    attributes=parse_attributes(attrs),
                )
            )

    if not features:
        raise ValueError(f"No CDS features found in GFF: {path}")

    return features


def normalize_identifier(value: str) -> str:
    value = value.strip().lstrip(">")
    value = value.split()[0]
    for prefix in ("cds-", "gene-", "rna-"):
        if value.startswith(prefix):
            value = value[len(prefix):]
    return value


def candidate_identifiers(feature: CDS) -> set[str]:
    values: set[str] = set()
    for key in ("ID", "protein_id", "Name", "locus_tag", "Parent"):
        raw = feature.attributes.get(key, "")
        if not raw:
            continue
        for value in raw.split(","):
            value = value.strip()
            if value:
                values.add(value)
                values.add(normalize_identifier(value))
    return {value for value in values if value}


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return value.strip("._") or "context"


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(
        sequence[start : start + width]
        for start in range(0, len(sequence), width)
    )


def main() -> None:
    args = parse_args()

    if args.flank_genes < 0:
        raise ValueError("--flank-genes must be zero or greater.")

    genome_manifest = args.genome_manifest.resolve()
    targets_path = args.targets.resolve()
    output_dir = args.output_dir.resolve()

    genomes = read_genome_manifest(genome_manifest)
    targets = read_targets(targets_path)

    unknown_genomes = sorted(
        {target["genome"] for target in targets} - set(genomes)
    )
    if unknown_genomes:
        raise ValueError(
            "Targets refer to genomes absent from the manifest: "
            + ", ".join(unknown_genomes)
        )

    if output_dir.exists() and any(output_dir.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Output directory exists and is not empty: {output_dir}"
            )
        shutil.rmtree(output_dir)

    fasta_dir = output_dir / "regions_fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    targets_by_genome: dict[str, list[dict[str, str]]] = defaultdict(list)
    for target in targets:
        targets_by_genome[target["genome"]].append(target)

    manifest_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    combined_records: list[tuple[str, str]] = []
    missing_targets: list[str] = []
    full_contexts = 0

    for genome in sorted(targets_by_genome):
        fasta_path, gff_path = genomes[genome]
        contigs = read_fasta(fasta_path)
        features = read_gff(gff_path)

        by_contig: dict[str, list[CDS]] = defaultdict(list)
        for feature in features:
            if feature.contig not in contigs:
                raise ValueError(
                    f"GFF contig {feature.contig!r} is absent from {fasta_path}"
                )
            by_contig[feature.contig].append(feature)

        for contig_features in by_contig.values():
            contig_features.sort(
                key=lambda feature: (feature.left, feature.right)
            )

        identifier_index: dict[str, list[CDS]] = defaultdict(list)
        for feature in features:
            for identifier in candidate_identifiers(feature):
                identifier_index[identifier].append(feature)

        for target in targets_by_genome[genome]:
            protein_id = target["protein_id"]
            lookup_id = target["gff_id"] or protein_id
            lookup_candidates = {
                lookup_id,
                normalize_identifier(lookup_id),
            }

            matches: list[CDS] = []
            seen_features: set[CDS] = set()
            for identifier in lookup_candidates:
                for feature in identifier_index.get(identifier, []):
                    if feature not in seen_features:
                        seen_features.add(feature)
                        matches.append(feature)

            if len(matches) == 0:
                missing_targets.append(f"{genome}\t{protein_id}\t{lookup_id}")
                continue
            if len(matches) > 1:
                raise ValueError(
                    f"Target {genome}|{protein_id} matches multiple CDS features. "
                    "Supply an unambiguous gff_id."
                )

            target_feature = matches[0]
            contig_features = by_contig[target_feature.contig]
            target_index = contig_features.index(target_feature)

            first_index = max(0, target_index - args.flank_genes)
            last_index = min(
                len(contig_features),
                target_index + args.flank_genes + 1,
            )
            context_features = contig_features[first_index:last_index]

            region_start = min(feature.left for feature in context_features)
            region_end = max(feature.right for feature in context_features)
            region_sequence = contigs[target_feature.contig][
                region_start - 1 : region_end
            ]

            expected_genes = 2 * args.flank_genes + 1
            is_full = len(context_features) == expected_genes
            if is_full:
                full_contexts += 1

            region_id = safe_name(f"{genome}__{protein_id}")
            output_fasta = fasta_dir / f"{region_id}.fna"

            header = (
                f"{region_id} genome={genome} target={protein_id} "
                f"contig={target_feature.contig} region={region_start}-{region_end} "
                f"target_strand={target_feature.strand}"
            )

            with output_fasta.open("w", encoding="utf-8") as handle:
                handle.write(f">{header}\n")
                handle.write(wrap(region_sequence))
                handle.write("\n")

            combined_records.append((header, region_sequence))

            manifest_rows.append(
                {
                    "region_id": region_id,
                    "genome": genome,
                    "protein_id": protein_id,
                    "gff_lookup_id": lookup_id,
                    "contig": target_feature.contig,
                    "region_start": region_start,
                    "region_end": region_end,
                    "region_length": len(region_sequence),
                    "target_start": target_feature.left,
                    "target_end": target_feature.right,
                    "target_strand": target_feature.strand,
                    "genes_in_context": len(context_features),
                    "expected_full_context_genes": expected_genes,
                    "full_context": is_full,
                    "output_fasta": str(output_fasta),
                }
            )

            for local_index, feature in enumerate(context_features):
                genomic_relative = first_index + local_index - target_index
                target_relative = (
                    genomic_relative
                    if target_feature.strand != "-"
                    else -genomic_relative
                )

                identifiers = sorted(candidate_identifiers(feature))
                gene_rows.append(
                    {
                        "region_id": region_id,
                        "genome": genome,
                        "target_protein_id": protein_id,
                        "contig": feature.contig,
                        "gene_start": feature.left,
                        "gene_end": feature.right,
                        "gene_strand": feature.strand,
                        "relative_position_genomic": genomic_relative,
                        "relative_position_target_orientation": target_relative,
                        "is_target": feature == target_feature,
                        "feature_identifiers": ";".join(identifiers),
                        "gff_attributes": feature.attributes_raw,
                    }
                )

    if missing_targets:
        missing_path = output_dir / "missing_targets.tsv"
        with missing_path.open("w", encoding="utf-8") as handle:
            handle.write("genome\tprotein_id\tgff_lookup_id\n")
            for row in missing_targets:
                handle.write(row + "\n")
        raise ValueError(
            f"{len(missing_targets)} targets were not found in the GFF files. "
            f"See: {missing_path}"
        )

    combined_fasta = output_dir / "all_MHC_contexts.fna"
    with combined_fasta.open("w", encoding="utf-8") as handle:
        for header, sequence in combined_records:
            handle.write(f">{header}\n")
            handle.write(wrap(sequence))
            handle.write("\n")

    manifest_path = output_dir / "MHC_context_manifest.tsv"
    manifest_fields = [
        "region_id",
        "genome",
        "protein_id",
        "gff_lookup_id",
        "contig",
        "region_start",
        "region_end",
        "region_length",
        "target_start",
        "target_end",
        "target_strand",
        "genes_in_context",
        "expected_full_context_genes",
        "full_context",
        "output_fasta",
    ]
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=manifest_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    genes_path = output_dir / "MHC_context_genes.tsv"
    gene_fields = [
        "region_id",
        "genome",
        "target_protein_id",
        "contig",
        "gene_start",
        "gene_end",
        "gene_strand",
        "relative_position_genomic",
        "relative_position_target_orientation",
        "is_target",
        "feature_identifiers",
        "gff_attributes",
    ]
    with genes_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=gene_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(gene_rows)

    summary = {
        "targets": len(targets),
        "contexts_extracted": len(manifest_rows),
        "flank_genes_each_side": args.flank_genes,
        "full_contexts": full_contexts,
        "contig_boundary_contexts": len(manifest_rows) - full_contexts,
        "combined_context_fasta": str(combined_fasta),
        "context_manifest": str(manifest_path),
        "context_gene_table": str(genes_path),
    }
    summary_path = output_dir / "MHC_context_extraction_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print("MHC genomic-context extraction completed successfully.")
    print(f"Targets processed: {len(targets)}")
    print(f"Full contexts: {full_contexts}")
    print(
        "Contig-boundary contexts: "
        f"{len(manifest_rows) - full_contexts}"
    )
    print(f"Per-target FASTA directory: {fasta_dir}")
    print(f"Combined FASTA: {combined_fasta}")
    print(f"Context manifest: {manifest_path}")


if __name__ == "__main__":
    main()
