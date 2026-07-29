#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=all_MAGs_GToTree
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=40G
#SBATCH --time=120:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# DESCRIPTION
# Builds a genome-wide phylogeny for all species-level dereplicated MAGs
# recovered in this study using GToTree v1.8.6 and the combined bacterial
# and archaeal single-copy marker set.
#
# INPUT
# A text file containing one MAG FASTA path per line.
#
# OUTPUT
# A GToTree output directory containing the concatenated marker alignment,
# intermediate files and the inferred genome tree.
#
# USAGE
# sbatch 01_build_all_MAGs_phylogeny_gtotree.sh \
#   input_MAG_paths.txt \
#   GToTree_all_MAGs_output
#
# OPTIONAL ENVIRONMENT VARIABLES
# GTOTREE_ACTIVATE=/path/to/gtotree/environment/bin/activate
# EXPECTED_GENOME_COUNT=1081
#
# Set EXPECTED_GENOME_COUNT=0 to disable the count check.

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 <MAG_paths.txt> <output_directory>" >&2
    exit 1
fi

GENOME_LIST="$1"
OUTPUT_DIR="$2"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
EXPECTED_GENOME_COUNT="${EXPECTED_GENOME_COUNT:-1081}"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -s "$GENOME_LIST" ]]; then
    echo "ERROR: Genome list was not found or is empty: $GENOME_LIST" >&2
    exit 1
fi

if [[ -e "$OUTPUT_DIR" ]]; then
    echo "ERROR: Output directory already exists: $OUTPUT_DIR" >&2
    exit 1
fi

N_GENOMES=0
N_MISSING=0
N_DUPLICATES=0

while IFS= read -r GENOME || [[ -n "$GENOME" ]]; do
    GENOME="${GENOME%$'\r'}"

    [[ -z "$GENOME" ]] && continue
    [[ "$GENOME" =~ ^[[:space:]]*# ]] && continue

    N_GENOMES=$((N_GENOMES + 1))

    if [[ ! -s "$GENOME" ]]; then
        echo "MISSING: $GENOME" >&2
        N_MISSING=$((N_MISSING + 1))
    fi
done < "$GENOME_LIST"

N_UNIQUE="$(
    grep -v '^[[:space:]]*$' "$GENOME_LIST" \
        | grep -v '^[[:space:]]*#' \
        | sed 's/\r$//' \
        | sort -u \
        | wc -l \
        | tr -d ' '
)"

N_DUPLICATES=$((N_GENOMES - N_UNIQUE))

if [[ "$N_MISSING" -ne 0 ]]; then
    echo "ERROR: $N_MISSING genome files are missing or empty." >&2
    exit 1
fi

if [[ "$N_DUPLICATES" -ne 0 ]]; then
    echo "ERROR: $N_DUPLICATES duplicated genome paths were detected." >&2
    exit 1
fi

if [[ "$EXPECTED_GENOME_COUNT" -ne 0 ]] \
    && [[ "$N_GENOMES" -ne "$EXPECTED_GENOME_COUNT" ]]; then
    echo "ERROR: Expected $EXPECTED_GENOME_COUNT genomes, found $N_GENOMES." >&2
    exit 1
fi

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
echo "Building phylogeny for all dereplicated MAGs"
echo "Genome count: $N_GENOMES"
echo "Marker set: Bacteria_and_Archaea"
echo "Threads: $THREADS"
echo "Output directory: $OUTPUT_DIR"
echo "============================================================"

GToTree \
    -f "$GENOME_LIST" \
    -H Bacteria_and_Archaea \
    -o "$OUTPUT_DIR" \
    -j "$THREADS"

echo "============================================================"
echo "GToTree completed successfully."
echo "Genomes analysed: $N_GENOMES"
echo "Output directory: $OUTPUT_DIR"
echo "============================================================"
