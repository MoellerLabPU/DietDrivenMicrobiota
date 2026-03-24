# Fig S2, S3, S4 — Diversity Analysis

Supplementary Figures S2–S4 present community diversity analyses including relative abundance, alpha diversity, and differential abundance testing.

## Code

- **`diversity.Rmd`** — R Markdown that:
  1. Calculates genome-size-normalized relative abundances from metagenomic profiles
  2. Generates family-level stacked barplots across diet groups
  3. Computes Shannon diversity with pairwise Wilcoxon tests
  4. Runs ANCOM-BC2 differential abundance analysis (control and fat groups separately)
  5. Produces volcano plots of log fold changes

## Required Inputs

- Per-sample InStrain profiles (loaded externally)
- `fastq_stats.tsv` — read counts per sample for normalization
- `representative_genome_size.tsv` — genome lengths for size normalization
- `mouse_MAGs.contree` — phylogenetic tree
- Mouse metadata tables
