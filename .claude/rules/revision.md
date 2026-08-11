---
paths:
  - "Revision/**"
---

# Revision/ — July 2026 revision analyses (verified 2026-08-10)

Everything here targets the **revision run**:
`/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/AlleleFlux_revision/AlleleFlux`
— tighter QC than the published mapq20 run (22 paired pre_end MAGs vs 62; 39
MAGs in `allele_analysis_pre_end-fat_control/`), **parquet** allele_analysis
outputs (mapq20 is tsv.gz), and it ships ready-made
`p_value_summary/significant_sites_summary/` rollups that mapq20 lacks.

## AlleleFlux configs (`Revision/AlleleFlux/`)

- `alleleflux_config.yaml` — the real revision run (`run_name:
  "AlleleFlux_revision"`, root_dir above; `permutation.enabled: False`, so that
  block is inert).
- `alleleflux_config_perm1.yaml` — the permutation null in **BYO mode**:
  `permutation.enabled: True` + `permuted_metadata_dir` pointing at its own
  `permuted/perm1` root, `input.reuse_from` → the real run's `longitudinal/`
  dir, using the Fig-1 `perm_group_swap_set1.tsv` sheet. One BYO run = one
  sheet; for more nulls make more configs. Mechanics + footguns: memory note
  `alleleflux-permutation-runs`.
- `alleleflux_visualization_config_per_replicate.yaml` — copy of
  `figures/Fig2/alleleflux_visualization_config.yaml` with
  `combined_per_replicate: True`; its `output_dir` deliberately reuses the
  existing `plotting_SLG443_...` dir so cached `track_freqs` are reused, and
  its p_value inputs still point at the OLD mapq20 run (intentional).
- `slurm_scripts/` — sbatch wrappers for the above (`run_alleleflux.sh`,
  `run_alleleflux_perm1.sh`, `run_visualization_per_replicate.sh`).

## Notebooks (`Revision/AlleleFlux/notebooks/`)

- `alleleflux_scores_and_pvalue_heatmaps.qmd` — parallelism/divergence score
  plots + p-value heatmaps over the revision run.
- `per_litter_boxplots.ipynb` — per-litter anchor-allele boxplots from the
  5.2 GB SLG443_bin.96 long table; holds the repo's only HF/LF mapping
  (`GROUP_LABELS = {"fat": "HF", "control": "LF"}`).

## Sibling analyses

- `Revision/variable_sites/` — per-SGB variable-site counts and spacing
  (mapq20, Fig-1 SGBs). "Variable site" definitions and the three nested site
  sets: invoke the `alleleflux-internals` skill first.
- `Revision/relative_abundance/` — regression of AlleleFlux significance vs
  MAG relative abundance; its tables/figures live on scratch (see its README).
- `Revision/P.sartorii/` — isolate follow-up.
