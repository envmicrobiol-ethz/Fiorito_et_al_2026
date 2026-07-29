#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=external_MAGs_BLASTn
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=350:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
shopt -s nullglob

# DESCRIPTION
# Searches nucleotide sequences from the conserved EET gene cassettes against
# an external MAG dataset using BLASTn-plus v2.14.1.
#
# The script concatenates all MAG nucleotide FASTA files in the supplied
# directory, checks that contig identifiers are unique, builds a nucleotide
# BLAST database and runs BLASTn. Identity, query-coverage and gap filters are
# applied by 03_call_credible_EET_cassettes_from_BLASTn.py.
#
# INPUT
# 1. Nucleotide FASTA containing the EET cassette query genes.
# 2. Directory containing one external MAG nucleotide FASTA per genome.
# 3. Dataset label, for example SPRUCE or STORDALEN MIRE.
# 4. Output directory.
#
# OUTPUT
# Concatenated target FASTA, BLAST database and raw BLASTn hit table.
#
# USAGE
# sbatch 01_run_BLASTn_external_MAGs.sh \
#   EET_query_genes.fna \
#   external_MAG_directory \
#   SPRUCE \
#   SPRUCE_BLASTn

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch $0 <query.fna> <MAG_directory> <dataset_label> <output_dir>" >&2
    exit 1
fi

QUERY_FASTA="$1"
MAG_DIR="$2"
DATASET="$3"
OUTPUT_DIR="$4"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -s "$QUERY_FASTA" ]]; then
    echo "ERROR: Query FASTA was not found or is empty: $QUERY_FASTA" >&2
    exit 1
fi

if [[ ! -d "$MAG_DIR" ]]; then
    echo "ERROR: MAG directory does not exist: $MAG_DIR" >&2
    exit 1
fi

if [[ -z "$DATASET" ]]; then
    echo "ERROR: Dataset label cannot be empty." >&2
    exit 1
fi

if [[ -e "$OUTPUT_DIR" ]] && [[ -n "$(find "$OUTPUT_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Output directory already exists and is not empty: $OUTPUT_DIR" >&2
    exit 1
fi

MAG_FASTAS=(
    "$MAG_DIR"/*.fna
    "$MAG_DIR"/*.fa
    "$MAG_DIR"/*.fasta
)

if [[ ${#MAG_FASTAS[@]} -eq 0 ]]; then
    echo "ERROR: No .fna, .fa or .fasta MAG files were found in: $MAG_DIR" >&2
    exit 1
fi

###############################################################################
# SOFTWARE
###############################################################################

if ! command -v makeblastdb >/dev/null 2>&1 \
    || ! command -v blastn >/dev/null 2>&1; then

    if command -v module >/dev/null 2>&1; then
        module load stack/2024-06
        module load gcc/12.2.0
        module load blast-plus/2.14.1
    fi
fi

if ! command -v makeblastdb >/dev/null 2>&1; then
    echo "ERROR: makeblastdb was not found in PATH." >&2
    exit 1
fi

if ! command -v blastn >/dev/null 2>&1; then
    echo "ERROR: blastn was not found in PATH." >&2
    exit 1
fi

###############################################################################
# OUTPUT PATHS
###############################################################################

FASTA_DIR="${OUTPUT_DIR}/database_fasta"
DB_DIR="${OUTPUT_DIR}/blast_database"
RAW_DIR="${OUTPUT_DIR}/raw_hits"

mkdir -p "$FASTA_DIR" "$DB_DIR" "$RAW_DIR"

COMBINED_FASTA="${FASTA_DIR}/${DATASET}_all_MAGs.fna"
CONTIG_IDS="${FASTA_DIR}/${DATASET}_contig_ids.txt"
DUPLICATE_IDS="${FASTA_DIR}/${DATASET}_duplicate_contig_ids.txt"
DB_PREFIX="${DB_DIR}/${DATASET}_MAGs"
BLAST_OUTPUT="${RAW_DIR}/${DATASET}_EET_genes_BLASTn.tsv"

###############################################################################
# CONCATENATE TARGET MAG FASTAS
###############################################################################

printf '%s\n' "${MAG_FASTAS[@]}" | sort | while IFS= read -r FASTA; do
    cat "$FASTA"
done > "$COMBINED_FASTA"

if [[ ! -s "$COMBINED_FASTA" ]]; then
    echo "ERROR: Combined target FASTA was not generated." >&2
    exit 1
fi

grep '^>' "$COMBINED_FASTA" \
    | sed 's/^>//' \
    | awk '{print $1}' \
    > "$CONTIG_IDS"

N_CONTIGS="$(wc -l < "$CONTIG_IDS" | tr -d ' ')"

if [[ "$N_CONTIGS" -eq 0 ]]; then
    echo "ERROR: No contig identifiers were found in the target FASTA." >&2
    exit 1
fi

sort "$CONTIG_IDS" | uniq -d > "$DUPLICATE_IDS"

if [[ -s "$DUPLICATE_IDS" ]]; then
    echo "ERROR: Duplicate contig identifiers were found across the MAG FASTAs." >&2
    echo "Examples:" >&2
    head -20 "$DUPLICATE_IDS" >&2
    exit 1
fi

###############################################################################
# BUILD DATABASE AND RUN BLASTN
###############################################################################

echo "============================================================"
echo "Running external-dataset BLASTn search"
echo "Dataset: $DATASET"
echo "MAG files: ${#MAG_FASTAS[@]}"
echo "Target contigs: $N_CONTIGS"
echo "Query FASTA: $QUERY_FASTA"
echo "Output: $BLAST_OUTPUT"
echo "============================================================"

makeblastdb \
    -in "$COMBINED_FASTA" \
    -dbtype nucl \
    -parse_seqids \
    -out "$DB_PREFIX"

blastn \
    -query "$QUERY_FASTA" \
    -db "$DB_PREFIX" \
    -out "$BLAST_OUTPUT" \
    -evalue 1e-10 \
    -max_target_seqs 5000 \
    -num_threads "$THREADS" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

if [[ ! -f "$BLAST_OUTPUT" ]]; then
    echo "ERROR: BLASTn output was not created: $BLAST_OUTPUT" >&2
    exit 1
fi

N_HITS="$(wc -l < "$BLAST_OUTPUT" | tr -d ' ')"

echo "============================================================"
echo "BLASTn completed successfully."
echo "Raw hits: $N_HITS"
echo "Raw output: $BLAST_OUTPUT"
echo "============================================================"
