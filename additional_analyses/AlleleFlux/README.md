# AlleleFlux revision run (breadth ≥ 0.5, coverage ≥ 1)

A second AlleleFlux run over the diet-manipulation data with stricter per-MAG QC
than the main-text run: coverage breadth ≥ 0.5 and mean depth ≥ 1 per MAG per
sample, MAPQ ≥ 20. It tests 22 MAGs in the paired `pre_end` divergence test
(the main-text run tests 62). Outputs live on scratch under
`revision_July_2026/AlleleFlux_revision/`.

**What the paper uses from this run:** the `pre_end` HF-vs-LF p-value heatmap,
shown as an Extended Data figure. Its self-contained code, the as-run config
and submission script, and the rendered figure are in
[`notebooks/AlleleFlux_breadth0.5_cov1/`](../../notebooks/AlleleFlux_breadth0.5_cov1/).
Everything else in this directory was produced alongside it and is not in the
paper.

## Files

- **`alleleflux_config.yaml`** — the run configuration (same file as the copy in
  `notebooks/AlleleFlux_breadth0.5_cov1/`). Enables LMM, CMH, dN/dS, regional
  contrast and gene scores in addition to the t-tests.
- **`alleleflux_config_perm1.yaml`** — a permutation null for the divergence
  score, in bring-your-own mode: the real run's profiles, QC and allele-frequency
  cache are reused via `input.reuse_from`, and group labels are relabelled from
  the Fig 1 `perm_group_swap_set1.tsv` sheet. One config = one null.
- **`slurm_scripts/run_alleleflux.sh`**, **`run_alleleflux_perm1.sh`** — the
  submission scripts for the real run and the null, as run.
- **`fix_bam_paths.py`** — rewrites the `bam_path` column of the original metadata
  sheet, whose BAM locations no longer exist; the config points at the corrected
  sheet it writes.
- **`notebooks/alleleflux_scores_and_pvalue_heatmaps.qmd`** — Quarto (R)
  notebook over the run: parallelism score plots (within-group), divergence
  score plots against the perm1 null, and the two-panel p-value heatmap for both
  periods (`pre_post` and `pre_end`). The `pre_end` heatmap is the one in the
  paper; the trimmed notebook under `notebooks/` draws only that panel.
