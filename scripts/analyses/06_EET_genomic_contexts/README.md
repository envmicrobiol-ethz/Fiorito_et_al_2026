## Scripts

- 01_extract_MHC_genomic_contexts.py: extracts nucleotide regions spanning each selected MHC and up to five neighboring genes on each side from MAGs or reference genomes. Script used for Acidobacteriota and Verrucomicrobiota MAGs coming from peat metagenomes of this study and electroactive reference microbial genomes.
- 02_export_MHC_context_nucleotide_FASTAs.py: retrieves the nucleotide sequence of each selected MHC genomic context from the original MAG or reference genome. 
- 03_compare_MHC_contexts_clinker.sh: compares Acidobacteriota and Verrucomicrobiota MHC loci and across the complete context collection using clinker v0.0.32.
