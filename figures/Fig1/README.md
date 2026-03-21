# Fig 1C-D — Phylogeny with Divergence/Parallelism Score Heatmaps

Panels C and D display a MAG phylogeny alongside heatmaps of AlleleFlux divergence and parallelism scores, with per-MAG minimum BH-corrected q-values indicating statistical significance.

- **`scores_figs.Rmd`** — R Markdown that builds the phylogeny figure with faceted heatmap panels. Generates variants for each statistical test (paired t-test, Wilcoxon, LMM, CMH) and includes a t-test vs. Wilcoxon score comparison plot.

## Workflow

1. Loads AlleleFlux score tables (divergence and parallelism) and a minimum BH q-value summary (`minBH_summary_real.tsv`).
2. Reads the MAG phylogenetic tree (`mouse_MAGs.contree`) and prunes to MAGs with paired test results.
3. Converts the tree to ultrametric form and colors branches by phylum.
4. Adds faceted heatmap panels for divergence/parallelism scores (viridis color scale) and per-test minimum q-value tiles with significance stars (`*` q < 0.05, `**` q < 0.01, `***` q < 0.001).

## Required Inputs

- `scores_single_sample-pre_end-fat_control_{group}-MAGs.tsv` — Single-sample parallelism scores per group
- `scores_two_sample_paired-pre_end-fat_control-MAGs.tsv` — Paired divergence scores
- `scores_lmm-pre_end-fat_control-MAGs.tsv` / `scores_lmm_across_time-*.tsv` — LMM scores
- `scores_cmh-pre_end-fat_control-MAGs-{pre,end}.tsv` / `scores_cmh_across_time-*.tsv` — CMH scores
- `minBH_summary_real.tsv` — Minimum BH-corrected q-values per MAG across tests
- `mouse_MAGs.contree` — IQ-TREE phylogeny of MAGs

## R Dependencies

`ggtree`, `treeio`, `ape`, `tidyverse`, `ggnewscale`, `viridis`
