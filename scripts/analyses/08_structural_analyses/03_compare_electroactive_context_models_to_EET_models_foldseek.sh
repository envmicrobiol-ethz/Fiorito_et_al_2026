#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=Foldseek_contexts_vs_EET
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=40G
#SBATCH --time=300:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# DESCRIPTION
# Compares selected best AlphaFold3 models from electroactive reference
# microorganism genomic contexts against structural models of the peatland
# EET proteins using Foldseek v8-ef4e960.
#
# For this study, the query directory contained one selected best model for
# each of 55 proteins from five 11-gene reference contexts:
#   - Acidobacterium capsulatum
#   - Paludibaculum fermentans
#   - Pontiella desulfatans
#   - Pontiella sulfatireligans
#   - Geothrix fermentans
#
# The target directory contained 12 structural models representing the
# selected peatland EET proteins. Both directories must contain the final
# selected best models, not all alternative AlphaFold3 models.
#
# INPUT
# 1. Directory containing the 55 selected best reference-context models.
# 2. Directory containing the 12 selected peatland EET target models.
# 3. Output TSV path.
#
# OUTPUT
# A Foldseek all-hit table with structural similarity and alignment metrics.
#
# USAGE
# sbatch 03_compare_electroactive_context_models_to_EET_models_foldseek.sh \
#   best_reference_context_models \
#   peatland_EET_models \
#   foldseek_reference_contexts_vs_EET.tsv
#
# OPTIONAL ENVIRONMENT VARIABLES
# EXPECTED_QUERY_COUNT=55
# EXPECTED_TARGET_COUNT=12
#
# Set either expected count to 0 to disable that count check.

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch $0 <query_structure_dir> <target_structure_dir> <output.tsv>" >&2
    exit 1
fi

QUERY_DIR="$1"
TARGET_DIR="$2"
OUTPUT_TSV="$3"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
TMP_DIR="${OUTPUT_TSV}.tmp"

EXPECTED_QUERY_COUNT="${EXPECTED_QUERY_COUNT:-55}"
EXPECTED_TARGET_COUNT="${EXPECTED_TARGET_COUNT:-12}"

COVERAGE_MODE=2
MIN_COVERAGE=0.90

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -d "$QUERY_DIR" ]]; then
    echo "ERROR: Query structure directory does not exist: $QUERY_DIR" >&2
    exit 1
fi

if [[ ! -d "$TARGET_DIR" ]]; then
    echo "ERROR: Target structure directory does not exist: $TARGET_DIR" >&2
    exit 1
fi

N_QUERY="$(
    find "$QUERY_DIR" \
        -maxdepth 1 \
        -type f \
        \( \
            -name '*.pdb' \
            -o -name '*.pdb.gz' \
            -o -name '*.cif' \
            -o -name '*.mmcif' \
        \) \
        | wc -l \
        | tr -d ' '
)"

N_TARGET="$(
    find "$TARGET_DIR" \
        -maxdepth 1 \
        -type f \
        \( \
            -name '*.pdb' \
            -o -name '*.pdb.gz' \
            -o -name '*.cif' \
            -o -name '*.mmcif' \
        \) \
        | wc -l \
        | tr -d ' '
)"

if [[ "$N_QUERY" -eq 0 ]]; then
    echo "ERROR: No supported query structures were found in: $QUERY_DIR" >&2
    exit 1
fi

if [[ "$N_TARGET" -eq 0 ]]; then
    echo "ERROR: No supported target structures were found in: $TARGET_DIR" >&2
    exit 1
fi

if [[ "$EXPECTED_QUERY_COUNT" -ne 0 ]] && [[ "$N_QUERY" -ne "$EXPECTED_QUERY_COUNT" ]]; then
    echo "ERROR: Expected $EXPECTED_QUERY_COUNT query structures, found $N_QUERY." >&2
    exit 1
fi

if [[ "$EXPECTED_TARGET_COUNT" -ne 0 ]] && [[ "$N_TARGET" -ne "$EXPECTED_TARGET_COUNT" ]]; then
    echo "ERROR: Expected $EXPECTED_TARGET_COUNT target structures, found $N_TARGET." >&2
    exit 1
fi

if [[ -e "$OUTPUT_TSV" ]]; then
    echo "ERROR: Output already exists: $OUTPUT_TSV" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_TSV")"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

###############################################################################
# SOFTWARE
###############################################################################

if ! command -v foldseek >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load stack/2024-06
        module load gcc/12.2.0
        module load openmpi/4.1.6
        module load foldseek/8-ef4e960
    fi
fi

if ! command -v foldseek >/dev/null 2>&1; then
    echo "ERROR: Foldseek was not found in PATH." >&2
    exit 1
fi

###############################################################################
# RUN FOLDSEEK
###############################################################################

echo "============================================================"
echo "Running Foldseek reference-context versus EET comparison"
echo "Selected best query models: $N_QUERY"
echo "Selected EET target models: $N_TARGET"
echo "Coverage mode: $COVERAGE_MODE"
echo "Minimum query coverage: $MIN_COVERAGE"
echo "Query directory: $QUERY_DIR"
echo "Target directory: $TARGET_DIR"
echo "Output table: $OUTPUT_TSV"
echo "============================================================"

foldseek easy-search \
    "$QUERY_DIR" \
    "$TARGET_DIR" \
    "$OUTPUT_TSV" \
    "$TMP_DIR" \
    --cov-mode "$COVERAGE_MODE" \
    -c "$MIN_COVERAGE" \
    --max-seqs 1000 \
    --threads "$THREADS" \
    --format-output \
    "query,target,evalue,bits,prob,alntmscore,qtmscore,ttmscore,lddt,rmsd,alnlen,qlen,tlen,qcov,tcov,fident"

if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo "ERROR: Foldseek output was not created: $OUTPUT_TSV" >&2
    exit 1
fi

N_HITS="$(wc -l < "$OUTPUT_TSV" | tr -d ' ')"

echo "============================================================"
echo "Foldseek comparison completed successfully."
echo "Selected best query models: $N_QUERY"
echo "Selected EET target models: $N_TARGET"
echo "Reported query-target hits: $N_HITS"
echo "Output table: $OUTPUT_TSV"
echo "============================================================"
