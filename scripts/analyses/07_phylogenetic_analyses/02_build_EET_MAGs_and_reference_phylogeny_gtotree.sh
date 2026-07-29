#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=EET_MAGs_references_GToTree
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=40G
#SBATCH --time=120:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# DESCRIPTION
# Builds a phylogeny containing the EET-bearing MAGs from this study and
# electroactive reference genomes using GToTree v1.8.6 and the combined
# bacterial and archaeal single-copy marker set.
#
# INPUT
# 1. A tab-separated EET MAG table whose first column contains MAG names or
#    paths. Additional columns, such as EET1-EET4 assignments, are retained
#    only in the user-supplied source table and are not passed to GToTree.
# 2. A directory containing all species-level dereplicated MAG FASTA files.
# 3. A directory containing the electroactive reference genome FASTA files.
# 4. A run directory for validation files and GToTree output.
#
# OUTPUT
# Resolved MAG and reference path lists, a combined GToTree input list,
# a genome manifest and the inferred genome tree.
#
# USAGE
# sbatch 02_build_EET_MAGs_and_reference_phylogeny_gtotree.sh \
#   EET_MAGs.tsv \
#   dereplicated_MAG_directory \
#   reference_genome_directory \
#   output_run_directory
#
# OPTIONAL ENVIRONMENT VARIABLES
# GTOTREE_ACTIVATE=/path/to/gtotree/environment/bin/activate
# EXPECTED_MAG_COUNT=120
# EXPECTED_REFERENCE_COUNT=60
#
# Set either expected count to 0 to disable that count check.

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch $0 <EET_MAGs.tsv> <MAG_dir> <reference_dir> <run_dir>" >&2
    exit 1
fi

EET_MAG_TABLE="$1"
DEREP_MAG_DIR="$2"
REFERENCE_DIR="$3"
RUN_DIR="$4"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
EXPECTED_MAG_COUNT="${EXPECTED_MAG_COUNT:-120}"
EXPECTED_REFERENCE_COUNT="${EXPECTED_REFERENCE_COUNT:-60}"

MAG_PATHS="${RUN_DIR}/EET_MAG_paths.txt"
REFERENCE_PATHS="${RUN_DIR}/reference_genome_paths.txt"
INPUT_LIST="${RUN_DIR}/GToTree_input_genomes.txt"
MANIFEST="${RUN_DIR}/genome_manifest.tsv"
MISSING_MAGS="${RUN_DIR}/missing_EET_MAGs.txt"
DUPLICATE_BASENAMES="${RUN_DIR}/duplicate_basenames.txt"
OUTPUT_DIR="${RUN_DIR}/output"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -s "$EET_MAG_TABLE" ]]; then
    echo "ERROR: EET MAG table was not found or is empty: $EET_MAG_TABLE" >&2
    exit 1
fi

if [[ ! -d "$DEREP_MAG_DIR" ]]; then
    echo "ERROR: Dereplicated MAG directory does not exist: $DEREP_MAG_DIR" >&2
    exit 1
fi

if [[ ! -d "$REFERENCE_DIR" ]]; then
    echo "ERROR: Reference-genome directory does not exist: $REFERENCE_DIR" >&2
    exit 1
fi

if [[ -e "$RUN_DIR" ]] && [[ -n "$(find "$RUN_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Run directory already exists and is not empty: $RUN_DIR" >&2
    exit 1
fi

mkdir -p "$RUN_DIR"

: > "$MAG_PATHS"
: > "$REFERENCE_PATHS"
: > "$MISSING_MAGS"
: > "$DUPLICATE_BASENAMES"

###############################################################################
# RESOLVE EET MAG PATHS
###############################################################################

N_MAG_ROWS=0

while IFS= read -r LINE || [[ -n "$LINE" ]]; do
    LINE="${LINE%$'\r'}"

    [[ -z "$LINE" ]] && continue
    [[ "$LINE" =~ ^[[:space:]]*# ]] && continue

    MAG_NAME="${LINE%%$'\t'*}"
    MAG_NAME="$(
        printf '%s' "$MAG_NAME" \
            | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
    )"

    [[ -z "$MAG_NAME" ]] && continue

    N_MAG_ROWS=$((N_MAG_ROWS + 1))
    MAG_PATH=""

    if [[ "$MAG_NAME" = /* ]] && [[ -s "$MAG_NAME" ]]; then
        MAG_PATH="$MAG_NAME"
    elif [[ -s "${DEREP_MAG_DIR}/${MAG_NAME}" ]]; then
        MAG_PATH="${DEREP_MAG_DIR}/${MAG_NAME}"
    else
        for EXTENSION in .fa .fna .fasta; do
            if [[ -s "${DEREP_MAG_DIR}/${MAG_NAME}${EXTENSION}" ]]; then
                MAG_PATH="${DEREP_MAG_DIR}/${MAG_NAME}${EXTENSION}"
                break
            fi
        done
    fi

    if [[ -n "$MAG_PATH" ]]; then
        printf '%s\n' "$MAG_PATH" >> "$MAG_PATHS"
    else
        printf '%s\n' "$MAG_NAME" >> "$MISSING_MAGS"
    fi
done < "$EET_MAG_TABLE"

###############################################################################
# FIND REFERENCE GENOMES
###############################################################################

find "$REFERENCE_DIR" \
    -maxdepth 1 \
    -type f \
    \( \
        -name '*.fa' \
        -o -name '*.fna' \
        -o -name '*.fasta' \
    \) \
    -print \
    | sort \
    > "$REFERENCE_PATHS"

###############################################################################
# VALIDATE COUNTS
###############################################################################

N_MISSING="$(wc -l < "$MISSING_MAGS" | tr -d ' ')"
N_MAG_PATHS="$(wc -l < "$MAG_PATHS" | tr -d ' ')"
N_UNIQUE_MAGS="$(sort -u "$MAG_PATHS" | wc -l | tr -d ' ')"
N_REFERENCES="$(wc -l < "$REFERENCE_PATHS" | tr -d ' ')"

if [[ "$N_MISSING" -gt 0 ]]; then
    echo "ERROR: $N_MISSING EET MAGs could not be resolved." >&2
    cat "$MISSING_MAGS" >&2
    exit 1
fi

if [[ "$N_MAG_PATHS" -ne "$N_UNIQUE_MAGS" ]]; then
    echo "ERROR: Duplicated EET MAG paths were detected." >&2
    exit 1
fi

if [[ "$EXPECTED_MAG_COUNT" -ne 0 ]] \
    && [[ "$N_MAG_PATHS" -ne "$EXPECTED_MAG_COUNT" ]]; then
    echo "ERROR: Expected $EXPECTED_MAG_COUNT EET MAGs, found $N_MAG_PATHS." >&2
    exit 1
fi

if [[ "$EXPECTED_REFERENCE_COUNT" -ne 0 ]] \
    && [[ "$N_REFERENCES" -ne "$EXPECTED_REFERENCE_COUNT" ]]; then
    echo "ERROR: Expected $EXPECTED_REFERENCE_COUNT references, found $N_REFERENCES." >&2
    exit 1
fi

###############################################################################
# BUILD COMBINED INPUT AND MANIFEST
###############################################################################

cat "$MAG_PATHS" "$REFERENCE_PATHS" > "$INPUT_LIST"

N_TOTAL="$(wc -l < "$INPUT_LIST" | tr -d ' ')"
N_UNIQUE_PATHS="$(sort -u "$INPUT_LIST" | wc -l | tr -d ' ')"

if [[ "$N_TOTAL" -ne "$N_UNIQUE_PATHS" ]]; then
    echo "ERROR: Duplicated full genome paths were detected." >&2
    exit 1
fi

if grep -n '[[:space:]]' "$INPUT_LIST"; then
    echo "ERROR: One or more genome paths contain spaces or tabs." >&2
    exit 1
fi

sed 's|.*/||' "$INPUT_LIST" \
    | sort \
    | uniq -d \
    > "$DUPLICATE_BASENAMES"

if [[ -s "$DUPLICATE_BASENAMES" ]]; then
    echo "ERROR: Duplicated genome filenames were detected:" >&2
    cat "$DUPLICATE_BASENAMES" >&2
    exit 1
fi

printf 'genome_group\tgenome_name\tgenome_path\n' > "$MANIFEST"

while IFS= read -r GENOME_PATH; do
    printf 'EET_MAG\t%s\t%s\n' \
        "$(basename "$GENOME_PATH")" \
        "$GENOME_PATH"
done < "$MAG_PATHS" >> "$MANIFEST"

while IFS= read -r GENOME_PATH; do
    printf 'reference_genome\t%s\t%s\n' \
        "$(basename "$GENOME_PATH")" \
        "$GENOME_PATH"
done < "$REFERENCE_PATHS" >> "$MANIFEST"

###############################################################################
# SOFTWARE
###############################################################################

if [[ -n "${GTOTREE_ACTIVATE:-}" ]]; then
    if [[ ! -f "$GTOTREE_ACTIVATE" ]]; then
        echo "ERROR: GTOTREE_ACTIVATE does not exist: $GTOTREE_ACTIVATE" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "$GTOTREE_ACTIVATE"
fi

if ! command -v GToTree >/dev/null 2>&1; then
    echo "ERROR: GToTree was not found in PATH." >&2
    echo "Activate the GToTree v1.8.6 environment or set GTOTREE_ACTIVATE." >&2
    exit 1
fi

###############################################################################
# RUN GTOTREE
###############################################################################

echo "============================================================"
echo "Building EET MAG and reference-genome phylogeny"
echo "EET MAGs: $N_MAG_PATHS"
echo "Reference genomes: $N_REFERENCES"
echo "Total genomes: $N_TOTAL"
echo "Marker set: Bacteria_and_Archaea"
echo "Output directory: $OUTPUT_DIR"
echo "============================================================"

GToTree \
    -f "$INPUT_LIST" \
    -H Bacteria_and_Archaea \
    -o "$OUTPUT_DIR" \
    -j "$THREADS"

echo "============================================================"
echo "GToTree completed successfully."
echo "Input list: $INPUT_LIST"
echo "Manifest: $MANIFEST"
echo "Output directory: $OUTPUT_DIR"
echo "============================================================"
