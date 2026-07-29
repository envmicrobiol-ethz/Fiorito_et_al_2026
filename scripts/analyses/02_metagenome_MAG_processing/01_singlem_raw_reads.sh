#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=iMG_SingleM_raw
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G
#SBATCH --time=140:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

###############################################################################
# DESCRIPTION
# Estimates taxonomic profiles and microbial read fractions from raw
# paired-end metagenomic reads using SingleM.
#
# INPUT
# Sample directories containing:
#   <sample>_PE1.fq.gz
#   <sample>_PE2.fq.gz
#
# SingleM metapackage.
#
# OUTPUT
# Per-sample SingleM profiles, Krona reports and microbial-fraction tables.
#
# USAGE
# sbatch 01_singlem_raw_reads.sh \
#   input_directory \
#   output_directory \
#   SingleM_metapackage
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch $0 <input_directory> <output_directory> <SingleM_metapackage>" >&2
    exit 1
fi

BASE_DIR="$1"
OUT_DIR="$2"
METAPACKAGE="$3"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

if [ ! -d "$BASE_DIR" ]; then
    echo "ERROR: Input directory not found: $BASE_DIR" >&2
    exit 1
fi

if [ ! -f "$METAPACKAGE" ]; then
    echo "ERROR: SingleM metapackage not found: $METAPACKAGE" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate singlem_0.18.3

export SINGLEM_METAPACKAGE_PATH="$METAPACKAGE"

for SAMPLE_DIR in "$BASE_DIR"/*; do

    [ -d "$SAMPLE_DIR" ] || continue

    SAMPLE_NAME=$(basename "$SAMPLE_DIR")

    PE1="$SAMPLE_DIR/${SAMPLE_NAME}_PE1.fq.gz"
    PE2="$SAMPLE_DIR/${SAMPLE_NAME}_PE2.fq.gz"

    if [ ! -f "$PE1" ] || [ ! -f "$PE2" ]; then
        echo "WARNING: Paired-end reads not found for $SAMPLE_NAME; skipping."
        continue
    fi

    PROFILE="$OUT_DIR/sample_${SAMPLE_NAME}_profile.tsv"
    KRONA="$OUT_DIR/sample_${SAMPLE_NAME}_profile.html"
    MICROBIAL_FRACTION="$OUT_DIR/sample_${SAMPLE_NAME}_microbial_fraction.tsv"

    echo "Processing: $SAMPLE_NAME"

    singlem pipe \
        -1 "$PE1" \
        -2 "$PE2" \
        -p "$PROFILE" \
        --threads "$THREADS"

    if [ ! -s "$PROFILE" ]; then
        echo "ERROR: SingleM profile not generated for $SAMPLE_NAME." >&2
        continue
    fi

    singlem summarise \
        --input-taxonomic-profile "$PROFILE" \
        --output-taxonomic-profile-krona "$KRONA"

    singlem microbial_fraction \
        --forward "$PE1" \
        --reverse "$PE2" \
        -p "$PROFILE" \
        > "$MICROBIAL_FRACTION"

    echo "SingleM completed: $SAMPLE_NAME"

done

echo "All SingleM analyses completed."
