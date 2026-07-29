## Scripts

- `01_singlem_raw_reads.sh`: estimates taxonomic profiles and microbial read fractions from raw paired-end metagenomic reads using SingleM.
- `02_fastqc_raw_reads.sh`: performs FastQC quality assessment of raw paired-end metagenomic reads before trimming.
- `03_bbduk_trimming_and_fastqc.sh`: performs adapter removal and quality trimming using BBDuk, followed by FastQC assessment of the cleaned reads.
- `04_metaspades_assembly.sh`: assembles cleaned paired-end metagenomic reads for each sample using metaSPAdes with k-mer sizes 21, 33, 55, 77, 99 and 127.
- `05a_bbmap_assembly_stats.sh`: calculates per-sample metaSPAdes assembly statistics using BBMap `stats.sh`.
- `05b_seqkit_assembly_stats.sh`: calculates per-sample metaSPAdes assembly statistics using SeqKit.
- `06_filter_scaffolds_min1000.sh`: removes scaffolds shorter than 1,000 nucleotides from each metaSPAdes assembly using BBMap `reformat.sh`.
- `07_predict_genes_with_prodigal.sh`: predicts protein-coding genes from scaffolds ≥1,000 nucleotides using Prodigal v2.6.3 in metagenomic mode.
- `08_generate_differential_coverage_profiles.sh`: maps selected cleaned metagenomic read sets to each filtered assembly using BBMap with random assignment of ambiguously mapped reads, followed by BAM sorting and indexing with SAMtools.
- `08_differential_coverage_mapping_design.tsv`: defines the cross-sample read-to-assembly combinations used to generate differential scaffold-coverage profiles for genome binning.
