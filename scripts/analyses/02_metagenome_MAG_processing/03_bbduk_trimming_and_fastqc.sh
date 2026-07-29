#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=iMG_bbduk_qc
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=120:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

###############################################################################
# DESCRIPTION
# Trims paired-end metagenomic reads with BBDuk and runs FastQC on the
# resulting cleaned reads.
#
# INPUT
# Sample directories containing:
#   <sample>_PE1.fq.gz
#   <sample>_PE2.fq.gz
#
# Adapter FASTA file.
#
# OUTPUT
# Cleaned paired-end reads and per-sample post-trimming FastQC reports.
#
# USAGE
# sbatch 03_bbduk_trimming_and_fastqc.sh \
#   input_directory \
#   adapters.fa
###############################################################################

set -euo pipefail

module load bbmap/39.01

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate fastqc_0.12.1

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 <input_directory> <adapters.fa>" >&2
    exit 1
fi

BASE_DIR="$1"
ADAPTERS="$2"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

if [ ! -d "$BASE_DIR" ]; then
    echo "ERROR: Input directory not found: $BASE_DIR" >&2
    exit 1
fi

if [ ! -f "$ADAPTERS" ]; then
    echo "ERROR: Adapter file not found: $ADAPTERS" >&2
    exit 1
fi

for SAMPLE_DIR in "$BASE_DIR"/*; do

    [ -d "$SAMPLE_DIR" ] || continue

    SAMPLE_NAME=$(basename "$SAMPLE_DIR")

    PE1="$SAMPLE_DIR/${SAMPLE_NAME}_PE1.fq.gz"
    PE2="$SAMPLE_DIR/${SAMPLE_NAME}_PE2.fq.gz"

    OUT1="$SAMPLE_DIR/${SAMPLE_NAME}_PE1_trim_clean.fq.gz"
    OUT2="$SAMPLE_DIR/${SAMPLE_NAME}_PE2_trim_clean.fq.gz"

    QC_DIR="$SAMPLE_DIR/qc_after_trimming"

    if [ ! -f "$PE1" ] || [ ! -f "$PE2" ]; then
        echo "WARNING: Paired-end reads not found for $SAMPLE_NAME; skipping."
        continue
    fi

    mkdir -p "$QC_DIR"

    echo "Processing: $SAMPLE_NAME"

    bbduk.sh -Xmx20g \
        threads="$THREADS" \
        in1="$PE1" \
        in2="$PE2" \
        out1="$OUT1" \
        out2="$OUT2" \
        ref="$ADAPTERS" \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        trimpolyg=20 \
        trimq=20 \
        qtrim=rl \
        minlen=105 \
        tpe \
        tbo

    fastqc \
        "$OUT1" \
        "$OUT2" \
        --outdir "$QC_DIR" \
        --threads "$THREADS"

    echo "BBDuk and FastQC completed: $SAMPLE_NAME"

done

echo "All samples completed."
