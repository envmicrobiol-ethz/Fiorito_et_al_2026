#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=Foldseek_MHC_vs_PDB100
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=40G
#SBATCH --time=300:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# DESCRIPTION
# Searches the selected structural models of peat-metagenome MHC cluster
# representatives against a local PDB100 structural database using
# Foldseek v8-ef4e960.
#
# The query directory should contain the selected best structural model for
# each conserved MHC cluster representative. PDB, PDB.GZ, CIF and mmCIF query
# structures are accepted by the input check.
#
# INPUT
# 1. Directory containing the selected MHC representative structures.
# 2. Prefix of a local Foldseek PDB100 database.
# 3. Output TSV path.
#
# OUTPUT
# A Foldseek all-hit table reporting structural coverage, probability,
# E-value, aligned TM-score, LDDT and sequence-identity statistics.
#
# USAGE
# sbatch 02_search_MHC_structures_against_PDB100_foldseek.sh \
#   best_MHC_representative_models \
#   /path/to/foldseek/PDB100_database \
#   foldseek_MHC_representatives_vs_PDB100.tsv

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch $0 <query_structure_dir> <PDB100_db_prefix> <output.tsv>" >&2
    exit 1
fi

QUERY_DIR="$1"
PDB100_DB="$2"
OUTPUT_TSV="$3"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
TMP_DIR="${OUTPUT_TSV}.tmp"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -d "$QUERY_DIR" ]]; then
    echo "ERROR: Query structure directory does not exist: $QUERY_DIR" >&2
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

if [[ "$N_QUERY" -eq 0 ]]; then
    echo "ERROR: No supported query structures were found in: $QUERY_DIR" >&2
    exit 1
fi

if [[ ! -e "$PDB100_DB" ]] && [[ ! -e "${PDB100_DB}.dbtype" ]]; then
    echo "ERROR: Foldseek PDB100 database was not found: $PDB100_DB" >&2
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
# RUN FOLDSEEK AGAINST PDB100
###############################################################################

echo "============================================================"
echo "Running Foldseek against PDB100"
echo "Query structures: $N_QUERY"
echo "Query directory: $QUERY_DIR"
echo "PDB100 database: $PDB100_DB"
echo "Output table: $OUTPUT_TSV"
echo "Threads: $THREADS"
echo "============================================================"

foldseek easy-search \
    "$QUERY_DIR" \
    "$PDB100_DB" \
    "$OUTPUT_TSV" \
    "$TMP_DIR" \
    --threads "$THREADS" \
    --format-output \
    "query,target,qcov,tcov,alnlen,prob,evalue,alntmscore,lddt,pident,nident"

if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo "ERROR: Foldseek output was not created: $OUTPUT_TSV" >&2
    exit 1
fi

N_HITS="$(wc -l < "$OUTPUT_TSV" | tr -d ' ')"

echo "============================================================"
echo "Foldseek PDB100 search completed successfully."
echo "Query structures: $N_QUERY"
echo "Reported hits: $N_HITS"
echo "Output table: $OUTPUT_TSV"
echo "============================================================"
