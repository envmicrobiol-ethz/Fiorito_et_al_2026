#!/bin/bash

set -euo pipefail

###############################################################################
# 16S rRNA AMPLICON PROCESSING
#
# This script processes raw paired-end 16S rRNA amplicon reads, including:
#   1. Addition of sample identifiers to FASTQ headers
#   2. Raw-read quality control
#   3. PhiX removal
#   4. Primer trimming
#   5. Paired-end read merging
#   6. Quality filtering and dereplication
#   7. ASV inference and chimera removal
#   8. Target verification with Metaxa2
#   9. OTU clustering at 97% sequence identity
#  10. ASV and OTU count-table generation
#  11. Taxonomic classification against SILVA v138
#
# Input files must be named:
#   <sample>_1.fastq.gz
#   <sample>_2.fastq.gz
#
# Usage:
#   bash 01_process_16S_amplicon_reads.sh \
#       <raw_read_directory> \
#       <output_directory> \
#       <PhiX_Bowtie2_index> \
#       <SILVA_database_fasta> \
#       [threads]
#
# Example:
#   bash 01_process_16S_amplicon_reads.sh \
#       1_raw \
#       2_processed \
#       /path/to/phix \
#       /path/to/SILVA138_RESCRIPt.fasta \
#       8
###############################################################################

###############################################################################
# ARGUMENTS
###############################################################################

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage:"
    echo "  bash $0 <raw_dir> <output_dir> <PhiX_index> <SILVA_fasta> [threads]"
    exit 1
fi

RAW_DIR="$1"
OUTPUT_DIR="$2"
PHIX_INDEX="$3"
TAXONOMY_DB="$4"
THREADS="${5:-8}"

###############################################################################
# PRIMERS
###############################################################################

FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"

###############################################################################
# SOFTWARE ENVIRONMENTS
###############################################################################

CONDA_SH="${CONDA_SH:-/cluster/project/umbiol/software/miniconda3/etc/profile.d/conda.sh}"

FASTQC_ENV="fastqc_0.12.1"
CUTADAPT_ENV="cutadapt_4.9"
METAXA_ENV="metaxa_2.2.3"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -d "$RAW_DIR" ]]; then
    echo "[ERROR] Raw-read directory not found: $RAW_DIR" >&2
    exit 1
fi

if [[ ! -f "$TAXONOMY_DB" ]]; then
    echo "[ERROR] SILVA taxonomy database not found: $TAXONOMY_DB" >&2
    exit 1
fi

if [[ ! -f "$CONDA_SH" ]]; then
    echo "[ERROR] Conda initialization script not found: $CONDA_SH" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] Threads must be a positive integer." >&2
    exit 1
fi

# Convert paths to absolute paths before changing directory.
RAW_DIR="$(cd "$RAW_DIR" && pwd)"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
TAXONOMY_DB="$(readlink -f "$TAXONOMY_DB")"

PHIX_INDEX_DIR="$(cd "$(dirname "$PHIX_INDEX")" && pwd)"
PHIX_INDEX="${PHIX_INDEX_DIR}/$(basename "$PHIX_INDEX")"

###############################################################################
# HELPER FUNCTION
###############################################################################

append_labeled_fastq() {
    local input_file="$1"
    local sample_name="$2"
    local output_file="$3"

    if [[ "$input_file" == *.gz ]]; then
        gzip -cd -- "$input_file"
    else
        cat -- "$input_file"
    fi |
        awk -v sample="$sample_name" '
            NR % 4 == 1 {
                print $1 ";sample=" sample ";"
                next
            }
            {
                print
            }
        ' >> "$output_file"
}

###############################################################################
# PREPARE COMBINED FASTQ FILES
###############################################################################

echo "[INFO] Preparing combined FASTQ files."

mapfile -t R1_FILES < <(
    find "$RAW_DIR" -maxdepth 1 -type f \
        \( -name "*_1.fastq.gz" -o -name "*_1.fastq" \) |
        sort
)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No files matching *_1.fastq.gz or *_1.fastq were found." >&2
    exit 1
fi

RAW_R1="$OUTPUT_DIR/1.all.raw.R1.fastq"
RAW_R2="$OUTPUT_DIR/1.all.raw.R2.fastq"

# Empty existing combined files before appending.
: > "$RAW_R1"
: > "$RAW_R2"

for R1_FILE in "${R1_FILES[@]}"; do

    R1_BASENAME="$(basename "$R1_FILE")"

    if [[ "$R1_BASENAME" == *_1.fastq.gz ]]; then
        SAMPLE_NAME="${R1_BASENAME%_1.fastq.gz}"
        R2_FILE="$RAW_DIR/${SAMPLE_NAME}_2.fastq.gz"
    else
        SAMPLE_NAME="${R1_BASENAME%_1.fastq}"
        R2_FILE="$RAW_DIR/${SAMPLE_NAME}_2.fastq"
    fi

    if [[ ! -f "$R2_FILE" ]]; then
        echo "[ERROR] R2 file missing for sample: $SAMPLE_NAME" >&2
        exit 1
    fi

    echo "[INFO] Adding sample: $SAMPLE_NAME"

    append_labeled_fastq "$R1_FILE" "$SAMPLE_NAME" "$RAW_R1"
    append_labeled_fastq "$R2_FILE" "$SAMPLE_NAME" "$RAW_R2"
done

cd "$OUTPUT_DIR"

###############################################################################
# RAW-READ QUALITY CONTROL
###############################################################################

echo "[INFO] Running FastQC."

source "$CONDA_SH"
conda activate "$FASTQC_ENV"

fastqc \
    --threads "$THREADS" \
    1.all.raw.R1.fastq \
    1.all.raw.R2.fastq

conda deactivate

module load stack/.2024-06-silent gcc/12.2.0 vsearch/2.22.1

vsearch \
    --fastq_eestats 1.all.raw.R1.fastq \
    --output 1.all.raw.R1.eestats.txt

vsearch \
    --fastq_eestats 1.all.raw.R2.fastq \
    --output 1.all.raw.R2.eestats.txt

###############################################################################
# PHIX REMOVAL
###############################################################################

echo "[INFO] Removing PhiX reads."

module load stack/.2024-06-silent gcc/12.2.0 bowtie2/2.5.1-u2j3omo

bowtie2 \
    -x "$PHIX_INDEX" \
    -1 1.all.raw.R1.fastq \
    -2 1.all.raw.R2.fastq \
    --un-conc all.nophix \
    --fast \
    --threads "$THREADS" \
    > /dev/null

mv all.1.nophix 2.all.nophix.R1.fastq
mv all.2.nophix 2.all.nophix.R2.fastq

###############################################################################
# PRIMER TRIMMING
###############################################################################

echo "[INFO] Trimming primers."

source "$CONDA_SH"
conda activate "$CUTADAPT_ENV"

cutadapt \
    --cores "$THREADS" \
    -g "^${FORWARD_PRIMER}" \
    -G "^${REVERSE_PRIMER}" \
    --error-rate 0.15 \
    --minimum-length 1 \
    --output 3.all.trim.R1.fastq \
    --paired-output 3.all.trim.R2.fastq \
    2.all.nophix.R1.fastq \
    2.all.nophix.R2.fastq

conda deactivate

###############################################################################
# PAIRED-END READ MERGING
###############################################################################

echo "[INFO] Merging paired-end reads."

module load stack/.2024-06-silent gcc/12.2.0 vsearch/2.22.1

vsearch \
    --fastq_mergepairs 3.all.trim.R1.fastq \
    --reverse 3.all.trim.R2.fastq \
    --fastqout 4.all.merge.fastq \
    --fastaout 4.all.merge.fasta \
    --fastq_truncqual 5 \
    --fastq_allowmergestagger \
    --fastq_minovlen 20 \
    --fastq_minmergelen 250 \
    --threads "$THREADS"

vsearch \
    --fastq_eestats 4.all.merge.fastq \
    --output 4.all.merge.eestats.txt

###############################################################################
# QUALITY FILTERING AND DEREPLICATION
###############################################################################

echo "[INFO] Filtering and dereplicating merged reads."

vsearch \
    --fastq_filter 4.all.merge.fastq \
    --fastaout 5.all.eefilter.fasta \
    --fastq_maxee 2 \
    --threads "$THREADS"

vsearch \
    --derep_fulllength 5.all.eefilter.fasta \
    --sizeout \
    --relabel Uniq \
    --output 6.all.uniq.fasta

###############################################################################
# ASV INFERENCE AND CHIMERA REMOVAL
###############################################################################

echo "[INFO] Inferring ASVs."

vsearch \
    --cluster_unoise 6.all.uniq.fasta \
    --centroids 6.all.ASV.fasta \
    --uc 6.all.ASV.uc \
    --relabel ASV \
    --sizeorder \
    --sizein \
    --sizeout \
    --minsize 4 \
    --threads "$THREADS"

echo "[INFO] Removing chimeric ASVs."

vsearch \
    --uchime3_denovo 6.all.ASV.fasta \
    --nonchimeras 7.all.ASV_nochim.fasta \
    --uchimeout 7.all.ASV_nochim.txt \
    --uchimealns 7.all.ASV_chimaln.txt \
    --abskew 16 \
    --sizein \
    --sizeout

###############################################################################
# TARGET VERIFICATION WITH METAXA2
###############################################################################

echo "[INFO] Verifying 16S rRNA targets with Metaxa2."

source "$CONDA_SH"
conda activate "$METAXA_ENV"

metaxa2_x \
    -i 7.all.ASV_nochim.fasta \
    -o 8.all.ASV_metaxa \
    --allow_reorder F \
    --graphical F \
    --truncate F \
    --complement F \
    --cpu "$THREADS"

conda deactivate

###############################################################################
# REFORMAT AND SORT METAXA2 OUTPUT
###############################################################################

echo "[INFO] Reformatting Metaxa2 output."

awk -F "|" '{print $1}' \
    8.all.ASV_metaxa.extraction.fasta \
    > 8.all.ASV_metaxa.extraction.tmp.fasta

mv \
    8.all.ASV_metaxa.extraction.tmp.fasta \
    8.all.ASV_metaxa.extraction.fasta

module load stack/.2024-06-silent gcc/12.2.0 seqkit/0.10.1

seqkit sort \
    --by-name \
    --line-width 80 \
    8.all.ASV_metaxa.extraction.fasta \
    --out-file 8.all.ASV_metaxa.fasta

###############################################################################
# OTU CLUSTERING
###############################################################################

echo "[INFO] Clustering ASVs into OTUs at 97% sequence identity."

module load stack/.2024-06-silent gcc/12.2.0 vsearch/2.22.1

vsearch \
    --cluster_size 8.all.ASV_metaxa.fasta \
    --centroids 8.all.OTU_metaxa.fasta \
    --uc 8.all.OTU_metaxa.uc \
    --id 0.97 \
    --relabel OTU \
    --sizeorder \
    --sizein \
    --sizeout \
    --threads "$THREADS"

###############################################################################
# REMOVE SIZE ANNOTATIONS FROM FASTA HEADERS
###############################################################################

echo "[INFO] Removing size annotations from FASTA headers."

awk -F ";size" '{print $1}' \
    8.all.ASV_metaxa.fasta \
    > 8.all.ASV_metaxa.tmp.fasta

mv \
    8.all.ASV_metaxa.tmp.fasta \
    8.all.ASV_metaxa.fasta

awk -F ";size" '{print $1}' \
    8.all.OTU_metaxa.fasta \
    > 8.all.OTU_metaxa.tmp.fasta

mv \
    8.all.OTU_metaxa.tmp.fasta \
    8.all.OTU_metaxa.fasta

###############################################################################
# GENERATE ASV AND OTU TABLES
###############################################################################

echo "[INFO] Mapping merged reads to ASVs."

vsearch \
    --usearch_global 4.all.merge.fasta \
    --db 8.all.ASV_metaxa.fasta \
    --id 0.8 \
    --maxhits 10 \
    --maxaccepts 20 \
    --maxrejects 0 \
    --uc 9.all.ASV_map.uc \
    --matched 9.all.ASV_map.fasta \
    --otutabout 9.all.ASV_map.txt \
    --threads "$THREADS"

echo "[INFO] Mapping merged reads to OTUs."

vsearch \
    --usearch_global 4.all.merge.fasta \
    --db 8.all.OTU_metaxa.fasta \
    --id 0.8 \
    --maxhits 10 \
    --maxaccepts 20 \
    --maxrejects 0 \
    --uc 9.all.OTU_map.uc \
    --matched 9.all.OTU_map.fasta \
    --otutabout 9.all.OTU_map.txt \
    --threads "$THREADS"

sort -n -k1.4 \
    --output 9.all.ASV_map.txt \
    9.all.ASV_map.txt

sort -n -k1.4 \
    --output 9.all.OTU_map.txt \
    9.all.OTU_map.txt

###############################################################################
# TAXONOMIC CLASSIFICATION WITH SINTAX
###############################################################################

echo "[INFO] Assigning taxonomy with SINTAX."

vsearch \
    --sintax 8.all.ASV_metaxa.fasta \
    --db "$TAXONOMY_DB" \
    --tabbedout 9.all.ASV_tax.sintax.txt \
    --sintax_cutoff 0.8 \
    --strand plus \
    --threads "$THREADS"

vsearch \
    --sintax 8.all.OTU_metaxa.fasta \
    --db "$TAXONOMY_DB" \
    --tabbedout 9.all.OTU_tax.sintax.txt \
    --sintax_cutoff 0.8 \
    --strand plus \
    --threads "$THREADS"

{
    printf "ASV\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"

    sort -n -k1.4 9.all.ASV_tax.sintax.txt |
        awk '{print $1 "\t" $4}' |
        sed 's/,/\t/g'
} > 9.all.ASV_tax.sintax.rf.txt

{
    printf "OTU\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"

    sort -n -k1.4 9.all.OTU_tax.sintax.txt |
        awk '{print $1 "\t" $4}' |
        sed 's/,/\t/g'
} > 9.all.OTU_tax.sintax.rf.txt

###############################################################################
# TAXONOMIC CLASSIFICATION WITH LCA
###############################################################################

echo "[INFO] Assigning taxonomy with LCA."

vsearch \
    --usearch_global 8.all.ASV_metaxa.fasta \
    --db "$TAXONOMY_DB" \
    --id 0.97 \
    --maxaccepts 20 \
    --uc_allhits \
    --uc 9.all.ASV_tax.lca.uc \
    --lca_cutoff 0.9 \
    --lcaout 9.all.ASV_tax.lca.txt \
    --threads "$THREADS"

vsearch \
    --usearch_global 8.all.OTU_metaxa.fasta \
    --db "$TAXONOMY_DB" \
    --id 0.97 \
    --maxaccepts 20 \
    --uc_allhits \
    --uc 9.all.OTU_tax.lca.uc \
    --lca_cutoff 0.9 \
    --lcaout 9.all.OTU_tax.lca.txt \
    --threads "$THREADS"

{
    printf "ASV\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"

    sort -n -k1.4 9.all.ASV_tax.lca.txt |
        sed 's/,/\t/g'
} > 9.all.ASV_tax.lca.rf.txt

{
    printf "OTU\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"

    sort -n -k1.4 9.all.OTU_tax.lca.txt |
        sed 's/,/\t/g'
} > 9.all.OTU_tax.lca.rf.txt

###############################################################################
# FINAL QUALITY CONTROL
###############################################################################

echo "[INFO] Calculating final ASV and OTU sequence statistics."

module load stack/.2024-06-silent gcc/12.2.0 seqkit/0.10.1

seqkit stats \
    8.all.ASV_metaxa.fasta \
    8.all.OTU_metaxa.fasta

echo "[INFO] 16S rRNA amplicon processing completed successfully."
echo "[INFO] Output directory: $OUTPUT_DIR"
