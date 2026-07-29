#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=ED_Figure7_clinker
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=120:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# DESCRIPTION
# Reproduces the four clinker comparisons used for ED Figure 7 by comparing
# representative EET1-EET4 loci from this study with characterized EET loci
# from model organisms.
#
# INPUT
# Four representative loci from this study and four reference GenBank files:
# Shewanella oneidensis MR-1, Geobacter sulfurreducens PCA OmcE,
# Shewanella baltica OS678 and Thermincola potens JR.
#
# OUTPUT
# Four clinker HTML alignments corresponding to EET1-EET4.
#
# USAGE
# sbatch 04_ED_Figure7_reference_context_comparison_clinker.sh \
#   output_directory \
#   EET1_representative.gbff \
#   EET2_representative.gbff \
#   EET3_representative.gbff \
#   EET4_representative.gbff \
#   Shewanella_oneidensis_MR1.gb \
#   Geobacter_sulfurreducens_PCA_OmcE.gb \
#   Shewanella_baltica_OS678.gb \
#   Thermincola_potens_JR.gb

if [[ $# -ne 9 ]]; then
    echo "Usage:" >&2
    echo "  sbatch $0 <output_dir> <EET1_rep> <EET2_rep> <EET3_rep> <EET4_rep> <S_oneidensis> <G_sulfurreducens_OmcE> <S_baltica> <T_potens>" >&2
    exit 1
fi

OUTPUT_DIR="$1"
REP1="$2"
REP2="$3"
REP3="$4"
REP4="$5"
SHEW_ONEIDENSIS="$6"
GEO_OMCE="$7"
SHEW_BALTICA="$8"
THERMINCOLA="$9"

if [[ -e "$OUTPUT_DIR" ]] && [[ -n "$(find "$OUTPUT_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Output directory already exists and is not empty: $OUTPUT_DIR" >&2
    exit 1
fi

for FILE in \
    "$REP1" \
    "$REP2" \
    "$REP3" \
    "$REP4" \
    "$SHEW_ONEIDENSIS" \
    "$GEO_OMCE" \
    "$SHEW_BALTICA" \
    "$THERMINCOLA"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: GenBank file not found or empty: $FILE" >&2
        exit 1
    fi
done

if [[ -n "${CLINKER_ACTIVATE:-}" ]]; then
    if [[ ! -f "$CLINKER_ACTIVATE" ]]; then
        echo "ERROR: CLINKER_ACTIVATE does not exist: $CLINKER_ACTIVATE" >&2
        exit 1
    fi
    # shellcheck disable=SC1090
    source "$CLINKER_ACTIVATE"
fi

if ! command -v clinker >/dev/null 2>&1; then
    echo "ERROR: clinker was not found in PATH." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Running ED Figure 7 comparison: EET1"
clinker \
    "$SHEW_BALTICA" \
    "$REP1" \
    "$GEO_OMCE" \
    --dont_set_origin \
    -ufo \
    -p "${OUTPUT_DIR}/ED_Figure7_EET1.html"

echo "Running ED Figure 7 comparison: EET2"
clinker \
    "$SHEW_ONEIDENSIS" \
    "$REP2" \
    "$SHEW_BALTICA" \
    --dont_set_origin \
    -ufo \
    -p "${OUTPUT_DIR}/ED_Figure7_EET2.html"

echo "Running ED Figure 7 comparison: EET3"
clinker \
    "$REP3" \
    "$SHEW_BALTICA" \
    --dont_set_origin \
    -p "${OUTPUT_DIR}/ED_Figure7_EET3.html"

echo "Running ED Figure 7 comparison: EET4"
clinker \
    "$THERMINCOLA" \
    "$REP4" \
    "$SHEW_BALTICA" \
    --dont_set_origin \
    -ufo \
    -p "${OUTPUT_DIR}/ED_Figure7_EET4.html"

echo "All ED Figure 7 clinker comparisons completed successfully."
echo "Output directory: $OUTPUT_DIR"
