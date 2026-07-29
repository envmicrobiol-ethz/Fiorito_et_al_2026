## Scripts

- `01_singlem_raw_reads.sh`: estimates taxonomic profiles and microbial read fractions from raw paired-end metagenomic reads using SingleM.
- `02_fastqc_raw_reads.sh`: performs FastQC quality assessment of raw paired-end metagenomic reads before trimming.
- `03_bbduk_trimming_and_fastqc.sh`: performs adapter removal and quality trimming using BBDuk, followed by FastQC assessment of the cleaned reads.
- `04_metaspades_assembly.sh`: assembles cleaned paired-end metagenomic reads for each sample using metaSPAdes with k-mer sizes 21, 33, 55, 77, 99 and 127.
- `05a_bbmap_assembly_stats.sh`: calculates per-sample metaSPAdes assembly statistics using BBMap `stats.sh`.
