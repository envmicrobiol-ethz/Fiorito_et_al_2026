# Gene catalogue and normalization

This directory contains the workflows used to construct the non-redundant
gene catalogue, calculate gene coverage, and normalize metagenomic and
metatranscriptomic measurements.

## Source and attribution

The gene catalogue, gene coverage, and normalization scripts were developed
by Taylor Priest and obtained from:

https://github.com/tpriest0/Profiling_metagenomes_and_mags

The version included here corresponds to commit:

`1e2aae10da1bf01ab51052cf48f221b92bdc17e7`

The original repository was released under the CC0 1.0 Universal license.
The scripts were copied and organized here for the analyses presented in
Fiorito et al. (2026).

## Contents

- `gene_catalogue_pipeline/`: constructs the non-redundant gene catalogue.
- `gene_coverage_pipeline/`: maps reads to the gene catalogue and calculates gene coverage.
- `calc_norm_cov_and_exp-functions.py`: shared normalization functions.
- `calc_norm_cov_and_exp-gene_catalog.py`: normalizes gene-catalogue coverage and expression.
- `calc_norm_cov_and_exp-mags.py`: normalizes MAG-level coverage and expression.
