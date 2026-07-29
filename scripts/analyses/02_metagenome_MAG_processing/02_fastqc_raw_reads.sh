#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=iMG_FastQC_raw
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G
#SBATCH --time=140:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

###############################################################################
# DESCRIPTION
# Runs FastQC on raw paired-end metagenomic reads before trimming.
#
# INPUT
# Sample directories containing:
#   <sample>_PE1.fq.gz
#   <sample>_PE2.fq.gz
#
# OUTPUT
# Per-sample FastQC reports in:
#   <sample>/qc_before_trimming/
#
# USAGE
# sbatch 02_fastqc_raw_reads.sh input_directory
###############################################################################

if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch $0 <input_directory>" >&2
    exit 1
fi

BASE_DIR="$1"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

if [ ! -d "$BASE_DIR" ]; then
    echo "ERROR: Input directory not found: $BASE_DIR" >&2
    exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate fastqc_0.12.1

for SAMPLE_DIR in "$BASE_DIR"/*; do

    [ -d "$SAMPLE_DIR" ] || continue

    SAMPLE_NAME=$(basename "$SAMPLE_DIR")

    PE1="$SAMPLE_DIR/${SAMPLE_NAME}_PE1.fq.gz"
    PE2="$SAMPLE_DIR/${SAMPLE_NAME}_PE2.fq.gz"

    if [ ! -f "$PE1" ] || [ ! -f "$PE2" ]; then
        echo "WARNING: Paired-end reads not found for $SAMPLE_NAME; skipping."
        continue
    fi

    QC_DIR="$SAMPLE_DIR/qc_before_trimming"

    mkdir -p "$QC_DIR"

    echo "Running FastQC: $SAMPLE_NAME"

    fastqc \
        "$PE1" \
        "$PE2" \
        --outdir "$QC_DIR" \
        --threads "$THREADS"

    echo "FastQC completed: $SAMPLE_NAME"

done

echo "All FastQC analyses completed."
