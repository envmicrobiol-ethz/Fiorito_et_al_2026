#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=MHC_contexts_clinker
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=120:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
shopt -s nullglob

# DESCRIPTION
# Compares all Bakta-annotated MHC genomic contexts together using
# clinker v0.0.32, without separating Acidobacteriota and
# Verrucomicrobiota contexts.
#
# INPUT
# A directory containing all genomic-context GenBank files.
#
# Generated from nucleotide context FASTAs produced by:
# scripts/analyses/06_EET_genomic_contexts/
# 02_export_MHC_context_nucleotide_FASTAs.py
#
# OUTPUT
# A single interactive clinker HTML alignment containing all genomic contexts.
#
# USAGE
# sbatch 03_compare_all_MHC_contexts_clinker.sh \
#   genomic_context_genbank_directory \
#   output_alignment.html

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 <genbank_directory> <output_alignment.html>" >&2
    exit 1
fi

GENBANK_DIR="$1"
OUTPUT_HTML="$2"

CLINKER_ACTIVATE="${CLINKER_ACTIVATE:-}"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -d "$GENBANK_DIR" ]]; then
    echo "ERROR: GenBank directory not found: $GENBANK_DIR" >&2
    exit 1
fi

if [[ -e "$OUTPUT_HTML" ]]; then
    echo "ERROR: Output file already exists: $OUTPUT_HTML" >&2
    exit 1
fi

###############################################################################
# SOFTWARE
###############################################################################

if [[ -n "$CLINKER_ACTIVATE" ]]; then

    if [[ ! -f "$CLINKER_ACTIVATE" ]]; then
        echo "ERROR: Clinker activation script not found:" >&2
        echo "$CLINKER_ACTIVATE" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "$CLINKER_ACTIVATE"
fi

if ! command -v clinker >/dev/null 2>&1; then
    echo "ERROR: clinker was not found in PATH." >&2
    exit 1
fi

###############################################################################
# COLLECT GENBANK FILES
###############################################################################

mapfile -d '' GENBANK_FILES < <(
    find "$GENBANK_DIR" \
        -maxdepth 1 \
        -type f \
        \( \
            -name "*.gb" \
            -o -name "*.gbk" \
            -o -name "*.gbff" \
        \) \
        -print0 |
        sort -z
)

if [[ ${#GENBANK_FILES[@]} -lt 2 ]]; then
    echo "ERROR: At least two GenBank files are required." >&2
    echo "Files found: ${#GENBANK_FILES[@]}" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_HTML")"

###############################################################################
# RUN CLINKER
###############################################################################

echo "============================================================"
echo "Starting clinker comparison"
echo "GenBank files: ${#GENBANK_FILES[@]}"
echo "Input directory: $GENBANK_DIR"
echo "Output alignment: $OUTPUT_HTML"
echo "============================================================"

printf '  %s\n' "${GENBANK_FILES[@]}"

clinker \
    "${GENBANK_FILES[@]}" \
    --dont_set_origin \
    -p "$OUTPUT_HTML"

if [[ ! -s "$OUTPUT_HTML" ]]; then
    echo "ERROR: Clinker HTML output was not generated." >&2
    exit 1
fi

echo "============================================================"
echo "Clinker comparison completed successfully."
echo "Genomic contexts compared: ${#GENBANK_FILES[@]}"
echo "Output alignment: $OUTPUT_HTML"
echo "============================================================"
