---
paths:
  - "additional_analyses/**"
  - "notebooks/AlleleFlux_breadth0.5_cov1/**"
  - "notebooks/relative_abundance/**"
  - "notebooks/standing_variation_vs_de_novo/**"
---

# July 2026 revision analyses (restructured 2026-09-14)

`Revision/` was split: analyses that made it into the paper moved to
`notebooks/`, the rest were renamed to `additional_analyses/`. Fig 2D (per-litter
violins) and Fig 2E (five-mice-excluded trajectories) moved to `figures/Fig2/`.

## The revision AlleleFlux run

`/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/AlleleFlux_revision/AlleleFlux`
— tighter QC than the published mapq20 run (22 paired pre_end MAGs vs 62; 39
MAGs in `allele_analysis_pre_end-fat_control/`), **parquet** allele_analysis
outputs (mapq20 is tsv.gz), and it ships ready-made
`p_value_summary/significant_sites_summary/` rollups that mapq20 lacks.

- `notebooks/AlleleFlux_breadth0.5_cov1/` — the paper's extended-data heatmap:
  as-run `alleleflux_config.yaml` + `run_alleleflux.sh` (copies; the script's
  `CONFIG=` line points at the config beside it, everything else as run),
  `pvalue_heatmap_pre_end.qmd` (pre_end, fat_control; no permutation null needed),
  and the rendered PNG. Render with the `alleleflux-R` conda env + Positron's
  bundled quarto; `fig_dir` is relative.
- `additional_analyses/AlleleFlux/` — the same config again, the perm1 null in
  **BYO mode** (`permutation.enabled: True`, `permuted_metadata_dir` → its own
  `permuted/perm1` root, `input.reuse_from` → the real run's `longitudinal/`;
  one BYO run = one sheet), both slurm scripts, `fix_bam_paths.py`, and the full
  `alleleflux_scores_and_pvalue_heatmaps.qmd` (parallelism + divergence-vs-null
  score plots + both-period heatmaps). Memory note `alleleflux-permutation-runs`
  has the footguns.

## Paper-included analyses under `notebooks/`

- `relative_abundance/` — regression of AlleleFlux significance vs MAG relative
  abundance (mapq20, pre_end, 62 divergence-tested MAGs). One README now holds
  the how-to-run and the plain-language glossary; `DESIGN.md` is gitignored and
  gone. Tables/figures live on scratch.
- `standing_variation_vs_de_novo/` — three chained steps:
  `div_and_hf_sites_in_both.py` (sites significant in both divergence and HF
  parallelism, per non-replacing species) → `pre_allele_presence.ipynb` (raw
  profiles for 2 species: was the HF END allele present at PRE?) →
  `baseline_presence/` (same question for every significant site via
  `alleleflux-baseline-presence`; its sbatch logs to a relative `logs/`).

## Not in the paper, under `additional_analyses/`

- `variable_sites_spacing/` — variable-site counts and spacing at contig, SGB and
  summary level (mapq20, the 62 divergence-tested Fig-1 SGBs). A "variable site"
  means **a site that was TESTED** under four definitions (`div`/`hf`/`lf`/`union`).
  README (with the column glossary merged in) + `METHODS.md`, both hand-maintained;
  the generator that once wrote METHODS.md is gone. Settled gotchas: the uniform
  null is `(L+1)/(n+1)` not `L/n`; `p_value_summary` needs `test_type` **and**
  `group_analyzed` pinned (the code asserts uniqueness rather than deduping).
- `P.sartorii_ani/` — isolate follow-up: Snakemake CheckM2 → skani, then
  `annotate_ani.py` labels pairs by mouse/timepoint/diet for PRE→END persistence.

The repo's only HF/LF mapping (`GROUP_LABELS = {"fat": "HF", "control": "LF"}`)
is in `figures/Fig2/per_litter_violins.ipynb`.
