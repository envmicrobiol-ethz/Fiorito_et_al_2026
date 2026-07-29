#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=protein_phylogeny
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# DESCRIPTION
# Builds a maximum-likelihood protein phylogeny for selected MHC or porin
# homologs using MAFFT L-INS-i, trimAl and IQ-TREE 3.
#
# The input FASTA is expected to contain the study proteins together with the
# publicly available homologs already retrieved and filtered using the
# sequence-identity and alignment-coverage criteria described in the Methods.
#
# INPUT
# A protein FASTA containing unique sequence identifiers.
#
# OUTPUT
# MAFFT alignment, trimAl-filtered alignment and IQ-TREE 3 result files.
#
# USAGE
# sbatch 03_build_MHC_or_porin_protein_phylogeny.sh \
#   selected_MHC_or_porin_homologs.faa \
#   output_directory \
#   tree_prefix

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch $0 <protein_fasta> <output_directory> <tree_prefix>" >&2
    exit 1
fi

INPUT_FASTA="$1"
OUTPUT_DIR="$2"
TREE_PREFIX="$3"

THREADS="${SLURM_CPUS_PER_TASK:-40}"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "ERROR: Protein FASTA was not found or is empty: $INPUT_FASTA" >&2
    exit 1
fi

if [[ -e "$OUTPUT_DIR" ]] && [[ -n "$(find "$OUTPUT_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Output directory already exists and is not empty: $OUTPUT_DIR" >&2
    exit 1
fi

N_SEQUENCES="$(grep -c '^>' "$INPUT_FASTA" || true)"

if [[ "$N_SEQUENCES" -lt 3 ]]; then
    echo "ERROR: At least three protein sequences are required." >&2
    exit 1
fi

N_UNIQUE_IDS="$(
    grep '^>' "$INPUT_FASTA" \
        | sed 's/^>//' \
        | awk '{print $1}' \
        | sort -u \
        | wc -l \
        | tr -d ' '
)"

if [[ "$N_UNIQUE_IDS" -ne "$N_SEQUENCES" ]]; then
    echo "ERROR: Duplicate FASTA identifiers were detected." >&2
    echo "Provide unique and traceable sequence identifiers before tree building." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

ALIGNED="${OUTPUT_DIR}/${TREE_PREFIX}_aligned.fasta"
TRIMMED="${OUTPUT_DIR}/${TREE_PREFIX}_aligned_trimmed.fasta"
IQTREE_PREFIX="${OUTPUT_DIR}/${TREE_PREFIX}"

###############################################################################
# SOFTWARE
###############################################################################

for PROGRAM in mafft trimal iqtree3; do
    if ! command -v "$PROGRAM" >/dev/null 2>&1; then
        echo "ERROR: Required program was not found in PATH: $PROGRAM" >&2
        exit 1
    fi
done

###############################################################################
# STEP 1: MAFFT L-INS-I
###############################################################################

echo "============================================================"
echo "Building protein phylogeny"
echo "Input sequences: $N_SEQUENCES"
echo "Input FASTA: $INPUT_FASTA"
echo "Output directory: $OUTPUT_DIR"
echo "Tree prefix: $TREE_PREFIX"
echo "============================================================"

echo "Step 1/3: MAFFT L-INS-i alignment"

mafft \
    --localpair \
    --maxiterate 1000 \
    --thread "$THREADS" \
    "$INPUT_FASTA" \
    > "$ALIGNED"

if [[ ! -s "$ALIGNED" ]]; then
    echo "ERROR: MAFFT alignment was not generated." >&2
    exit 1
fi

###############################################################################
# STEP 2: TRIMAL
###############################################################################

echo "Step 2/3: trimAl automated trimming"

trimal \
    -in "$ALIGNED" \
    -out "$TRIMMED" \
    -automated1

if [[ ! -s "$TRIMMED" ]]; then
    echo "ERROR: trimAl output was not generated." >&2
    exit 1
fi

###############################################################################
# STEP 3: IQ-TREE 3
###############################################################################

echo "Step 3/3: IQ-TREE 3 maximum-likelihood inference"

iqtree3 \
    -s "$TRIMMED" \
    -m MFP \
    -bb 1000 \
    -nt "$THREADS" \
    -pre "$IQTREE_PREFIX"

if [[ ! -s "${IQTREE_PREFIX}.treefile" ]]; then
    echo "ERROR: IQ-TREE treefile was not generated." >&2
    exit 1
fi

echo "============================================================"
echo "Protein phylogeny completed successfully."
echo "Alignment: $ALIGNED"
echo "Trimmed alignment: $TRIMMED"
echo "Maximum-likelihood tree: ${IQTREE_PREFIX}.treefile"
echo "Consensus tree: ${IQTREE_PREFIX}.contree"
echo "IQ-TREE report: ${IQTREE_PREFIX}.iqtree"
echo "============================================================"
