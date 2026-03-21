# Fig 3 — Gene-Level Score Enrichment and Hypergeometric Tests

Figure 3 presents gene-level enrichment analysis testing whether specific functional categories (COG annotations) are over-represented among genes with significant allele frequency changes.

## Code

- **`gene_scores_and_hypergeometric_test.Rmd`** — R Markdown that:
  1. Loads per-gene AlleleFlux scores (divergence, fat-parallelism, control-parallelism)
  2. Joins with reCOGnizer COG functional annotations
  3. Identifies genes with BH-significant sites
  4. Performs hypergeometric enrichment tests across functional categories
  5. Generates dot plots of enrichment results (with and without mobilome)

## Required Inputs

- Per-MAG gene score TSV files (control, fat, divergence)
- `reCOGnizer_results.tsv` — COG functional annotations
- `sweeps.tsv` — strain replacement results (for MAG filtering)
- `pre_end_summary_all_rows.tsv` — BH-corrected q-values per position
- `pre_end_summary_significant_rows.tsv` — significant rows after BH correction
