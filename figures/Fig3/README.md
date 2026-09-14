# Fig 3 — Gene-Level Score Enrichment and Hypergeometric Tests

Figure 3 presents gene-level enrichment analysis testing whether specific functional categories (COG annotations) are over-represented among genes with significant allele frequency changes, across the SGBs that showed no strain replacement under the high-fat diet.

## Code

- **`gene_scores_and_hypergeometric_test.Rmd`** — R Markdown that:
  1. Loads per-gene AlleleFlux scores (divergence, fat-parallelism, control-parallelism) and joins them with reCOGnizer COG functional annotations
  2. Restricts to a fixed list of 29 SGBs with no strain replacement in HF by conANI (the `keep_mag_ids` vector in the first chunk, annotated with figure numbering and species)
  3. Flags genes carrying BH-significant sites, per test
  4. Runs hypergeometric enrichment tests per functional category — for each test separately and for the combined set, each with and without the mobilome category — with BH correction within each run
  5. Writes the enrichment tables and draws q-value dot plots (combined, mobilome excluded)
  6. Draws per-gene Manhattan plots with a chosen protein description highlighted
  7. Draws the MAG phylogeny with a presence-absence grid of significant genes

## Required Inputs

All read by bare filename from the working directory:

- `gene_scores_control/`, `gene_scores_fat/`, `gene_scores_divergence/` — per-MAG gene score TSVs from AlleleFlux (`*_single_sample_control_gene_scores_combined.tsv`, `*_single_sample_fat_gene_scores_combined.tsv`, `*_two_sample_paired_gene_scores_combined.tsv`)
- `reCOGnizer_results.tsv` — COG functional annotations
- `pre_end_summary_all_rows.tsv` — BH-corrected q-values per position
- `pre_end_summary_significant_rows.tsv` — significant rows after BH correction
- `mouse_MAGs.contree` — IQ-TREE phylogeny of the MAGs (phylogeny section only)

## Outputs

- `fig3_output/` — enrichment result TSVs and q-value plots
- `manhattan_plots/` — one PDF per gene
- `tree_presence_absence_29mags.pdf` — phylogeny with presence-absence grid

## R Dependencies

`dplyr`, `readr`, `stringr`, `tidyr`, `purrr`, `forcats`, `tibble`, `ggplot2`, `ggrepel`, `writexl`, `ape`, `phytools`, `ggtree`
