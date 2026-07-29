#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=SignalP_MHC_representatives
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=200:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# DESCRIPTION
# Predicts signal peptides for the conserved MHC cluster representatives
# selected for structural analysis using SignalP 6.
#
# This script is intended for a FASTA containing one representative protein
# sequence per conserved MHC cluster, rather than all candidate MHC proteins.
#
# INPUT
# A protein FASTA containing the conserved cluster representatives.
#
# OUTPUT
# SignalP 6 results written to the requested output directory.
#
# USAGE
# sbatch 01_predict_signal_peptides_cluster_representatives.sh \
#   conserved_MHC_cluster_representatives.faa \
#   signalp_output_directory
#
# OPTIONAL ENVIRONMENT VARIABLE
# SIGNALP_ACTIVATE=/path/to/signalp/environment/bin/activate

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 <cluster_representatives.faa> <output_directory>" >&2
    exit 1
fi

INPUT_FASTA="$1"
OUTPUT_DIR="$2"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "ERROR: Input FASTA was not found or is empty: $INPUT_FASTA" >&2
    exit 1
fi

N_SEQUENCES="$(grep -c '^>' "$INPUT_FASTA" || true)"

if [[ "$N_SEQUENCES" -eq 0 ]]; then
    echo "ERROR: No protein records were found in: $INPUT_FASTA" >&2
    exit 1
fi

if [[ -e "$OUTPUT_DIR" ]] && [[ -n "$(find "$OUTPUT_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Output directory already exists and is not empty: $OUTPUT_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

###############################################################################
# SOFTWARE
###############################################################################

if [[ -n "${SIGNALP_ACTIVATE:-}" ]]; then
    if [[ ! -f "$SIGNALP_ACTIVATE" ]]; then
        echo "ERROR: SIGNALP_ACTIVATE does not exist: $SIGNALP_ACTIVATE" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "$SIGNALP_ACTIVATE"
fi

if ! command -v signalp6 >/dev/null 2>&1; then
    echo "ERROR: signalp6 was not found in PATH." >&2
    echo "Activate the SignalP 6 environment or set SIGNALP_ACTIVATE." >&2
    exit 1
fi

###############################################################################
# RUN SIGNALP 6
###############################################################################

echo "============================================================"
echo "Running SignalP 6"
echo "Conserved cluster representatives: $N_SEQUENCES"
echo "Input FASTA: $INPUT_FASTA"
echo "Output directory: $OUTPUT_DIR"
echo "============================================================"

signalp6 \
    -fmt none \
    -ff "$INPUT_FASTA" \
    -od "$OUTPUT_DIR"

echo "============================================================"
echo "SignalP 6 completed successfully."
echo "Sequences analysed: $N_SEQUENCES"
echo "Output directory: $OUTPUT_DIR"
echo "============================================================"
