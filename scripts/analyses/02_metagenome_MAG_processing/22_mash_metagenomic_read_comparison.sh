#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=iMG_Mash_reads
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=40G
#SBATCH --time=350:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
shopt -s nullglob

# DESCRIPTION
# Compares metagenomic read datasets from Värmland, Stordalen Mire and SPRUCE
# using Mash sketches and pairwise Mash distances.
#
# INPUT
# Directories containing metagenomic FASTQ files from the three datasets.
# Värmland reads may occur inside per-sample subdirectories.
#
# OUTPUT
# Individual read-file sketches, one combined sketch per dataset, a combined
# sketch containing all datasets and a pairwise Mash distance table.
#
# USAGE
# sbatch 22_mash_metagenomic_read_comparison.sh \
#   Varmland_reads_directory \
#   Stordalen_reads_directory \
#   SPRUCE_reads_directory \
#   output_directory

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch $0 <Varmland_reads_dir> <Stordalen_reads_dir> <SPRUCE_reads_dir> <output_dir>" >&2
    exit 1
fi

VARMLAND_READS_DIR="$1"
STORDALEN_READS_DIR="$2"
SPRUCE_READS_DIR="$3"
OUTPUT_DIR="$4"

MASH_ENV="${MASH_ENV:-mash_2.1}"

SKETCH_SIZE=10000
MIN_COPIES=2

###############################################################################
# INPUT VALIDATION
###############################################################################

for INPUT_DIR in \
    "$VARMLAND_READS_DIR" \
    "$STORDALEN_READS_DIR" \
    "$SPRUCE_READS_DIR"
do
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "ERROR: Input directory not found: $INPUT_DIR" >&2
        exit 1
    fi
done

if [[ -e "$OUTPUT_DIR" ]]; then
    echo "ERROR: Output directory already exists: $OUTPUT_DIR" >&2
    echo "Remove it or provide a different output directory." >&2
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda was not found." >&2
    exit 1
fi

###############################################################################
# SOFTWARE
###############################################################################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$MASH_ENV"

###############################################################################
# OUTPUT DIRECTORIES
###############################################################################

VARMLAND_SKETCH_DIR="${OUTPUT_DIR}/Varmland_sketches"
STORDALEN_SKETCH_DIR="${OUTPUT_DIR}/Stordalen_sketches"
SPRUCE_SKETCH_DIR="${OUTPUT_DIR}/SPRUCE_sketches"

mkdir -p \
    "$VARMLAND_SKETCH_DIR" \
    "$STORDALEN_SKETCH_DIR" \
    "$SPRUCE_SKETCH_DIR"

###############################################################################
# FUNCTION: FIND FASTQ FILES
###############################################################################

find_fastq_files() {

    local input_directory="$1"

    find "$input_directory" \
        -type f \
        \( \
            -name "*.fq" \
            -o -name "*.fq.gz" \
            -o -name "*.fastq" \
            -o -name "*.fastq.gz" \
        \) \
        -print0 |
        sort -z
}

###############################################################################
# FUNCTION: GENERATE INDIVIDUAL SKETCHES
###############################################################################

generate_sketches() {

    local dataset_name="$1"
    local input_directory="$2"
    local sketch_directory="$3"

    local read_files=()
    local read_file
    local relative_path
    local sketch_name
    local sketch_output
    local count=0

    mapfile -d '' read_files < <(
        find_fastq_files "$input_directory"
    )

    if [[ ${#read_files[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found for $dataset_name in: $input_directory" >&2
        exit 1
    fi

    echo "============================================================"
    echo "Generating Mash sketches: $dataset_name"
    echo "FASTQ files: ${#read_files[@]}"
    echo "Sketch size: $SKETCH_SIZE"
    echo "Minimum multiplicity: $MIN_COPIES"
    echo "============================================================"

    for read_file in "${read_files[@]}"; do

        relative_path="${read_file#"$input_directory"/}"

        sketch_name="$relative_path"
        sketch_name="${sketch_name%.fastq.gz}"
        sketch_name="${sketch_name%.fq.gz}"
        sketch_name="${sketch_name%.fastq}"
        sketch_name="${sketch_name%.fq}"

        # Preserve sample-directory information while producing safe filenames.
        sketch_name="${sketch_name//\//__}"
        sketch_name="${sketch_name// /_}"

        sketch_output="${sketch_directory}/${sketch_name}"

        if [[ -e "${sketch_output}.msh" ]]; then
            echo "ERROR: Duplicate sketch output: ${sketch_output}.msh" >&2
            exit 1
        fi

        echo "Sketching: $read_file"

        mash sketch \
            -r \
            -s "$SKETCH_SIZE" \
            -m "$MIN_COPIES" \
            "$read_file" \
            -o "$sketch_output"

        if [[ ! -s "${sketch_output}.msh" ]]; then
            echo "ERROR: Mash sketch was not generated: ${sketch_output}.msh" >&2
            exit 1
        fi

        count=$((count + 1))
    done

    echo "$dataset_name sketches generated: $count"
}

###############################################################################
# GENERATE SKETCHES
###############################################################################

generate_sketches \
    "Värmland" \
    "$VARMLAND_READS_DIR" \
    "$VARMLAND_SKETCH_DIR"

generate_sketches \
    "Stordalen Mire" \
    "$STORDALEN_READS_DIR" \
    "$STORDALEN_SKETCH_DIR"

generate_sketches \
    "SPRUCE" \
    "$SPRUCE_READS_DIR" \
    "$SPRUCE_SKETCH_DIR"

###############################################################################
# COMBINE SKETCHES WITHIN EACH DATASET
###############################################################################

VARMLAND_SKETCHES=("$VARMLAND_SKETCH_DIR"/*.msh)
STORDALEN_SKETCHES=("$STORDALEN_SKETCH_DIR"/*.msh)
SPRUCE_SKETCHES=("$SPRUCE_SKETCH_DIR"/*.msh)

VARMLAND_COMBINED="${OUTPUT_DIR}/Varmland_metagenomic_reads"
STORDALEN_COMBINED="${OUTPUT_DIR}/Stordalen_metagenomic_reads"
SPRUCE_COMBINED="${OUTPUT_DIR}/SPRUCE_metagenomic_reads"

mash paste \
    "$VARMLAND_COMBINED" \
    "${VARMLAND_SKETCHES[@]}"

mash paste \
    "$STORDALEN_COMBINED" \
    "${STORDALEN_SKETCHES[@]}"

mash paste \
    "$SPRUCE_COMBINED" \
    "${SPRUCE_SKETCHES[@]}"

###############################################################################
# COMBINE ALL DATASETS
###############################################################################

ALL_DATASETS="${OUTPUT_DIR}/all_metagenomic_reads"

mash paste \
    "$ALL_DATASETS" \
    "${VARMLAND_COMBINED}.msh" \
    "${STORDALEN_COMBINED}.msh" \
    "${SPRUCE_COMBINED}.msh"

if [[ ! -s "${ALL_DATASETS}.msh" ]]; then
    echo "ERROR: Combined Mash sketch was not generated." >&2
    exit 1
fi

###############################################################################
# CALCULATE PAIRWISE MASH DISTANCES
###############################################################################

DISTANCE_TABLE="${OUTPUT_DIR}/mash_read_distances.tsv"

mash dist \
    "${ALL_DATASETS}.msh" \
    "${ALL_DATASETS}.msh" \
    > "$DISTANCE_TABLE"

if [[ ! -s "$DISTANCE_TABLE" ]]; then
    echo "ERROR: Mash distance table was not generated." >&2
    exit 1
fi

echo "============================================================"
echo "Mash read comparison completed successfully."
echo "Värmland sketches: ${#VARMLAND_SKETCHES[@]}"
echo "Stordalen Mire sketches: ${#STORDALEN_SKETCHES[@]}"
echo "SPRUCE sketches: ${#SPRUCE_SKETCHES[@]}"
echo "Combined sketch: ${ALL_DATASETS}.msh"
echo "Pairwise distances: $DISTANCE_TABLE"
echo "============================================================"
