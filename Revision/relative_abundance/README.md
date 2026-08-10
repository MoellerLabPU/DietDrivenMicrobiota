# Relative abundance vs. AlleleFlux significance

Is the significance AlleleFlux reports for a MAG associated with **how abundant** that MAG is, or with
**how much its abundance changed** between PRE and END?

Scope: the `pre_end-fat_control` comparison of the `AlleleFlux_mapq20` run, `single_sample_tTest`
(parallelism, within one arm) and `two_sample_paired_tTest` (divergence, arms paired by cage).

`DESIGN.md` is the full write-up — method, validation of the inputs, and a step-by-step comparison
against the DRIDO `heatmap_support_summary.py` pipeline. Start there for the *why*; this file is the
*how to run*.

## Layout

```
DESIGN.md                        method, validation, DRIDO comparison
build_rel_abundance.py           Table S3 + metadata + cell stats -> tables/rel_abundance_by_cell.tsv
rel_abundance_regression.ipynb   the regression, tables and figures
tables/                          on scratch (see Regenerating)
  rel_abundance_by_cell.tsv      202 rows: one per heatmap cell
  regression_results.tsv         63 rows: one per fit
figures/                         PDFs, on scratch
```

## Regenerating

**Step 0 — the significance table** (only if it is missing; it is an input, not an output here):

```bash
python /home/su2806/alleleflux_benchmark/drido/significant_sites_heatmap/significant_sites_summary.py \
  --input-dir /scratch/gpfs/AMOELLER/diet_manip/AlleleFlux_mapq20/longitudinal/p_value_summary \
  --outdir    /scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/relative_abundance/significant_sites_summary
```

Reads ~3.5 GB and also writes a ~343 MB `significant_sites_sig_sites.tsv` that nothing here consumes.

**Step 1 — the abundance table:**

```bash
python build_rel_abundance.py          # all inputs have defaults; --help lists the overrides
```

**Step 2 — the regression:**

```bash
jupyter nbconvert --to notebook --execute --inplace rel_abundance_regression.ipynb
```

All paths in the notebook are absolute: inputs and outputs both live under
`/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/relative_abundance/` (`tables/` and
`figures/` there). Executing rewrites `tables/regression_results.tsv` in that directory.

## What the numbers mean

**The unit of observation is the cage, not the mouse.** AlleleFlux averages every mouse sharing a
`replicate` into one value before testing, so each test here has *n = 8*, not 15 or 30. Abundance is
aggregated the same way — per mouse, then into cages, then across cages. Averaging over mice instead
would weight the one-mouse cage differently and describe a different population than the p-value
does. This is the single most important thing to preserve if the code is adapted.

**Relative abundance is the published Table S3, unchanged** — a per-sample composition over all 160
genomes, summing to 100%. It is *not* renormalized to the subset AlleleFlux tested, which is where
this differs most from DRIDO.

**Every model is univariate.** Abundance and change are tested separately, never in one formula, so
each β is a total association. The absolute change/contrast predictors themselves track abundance
level (Spearman ρ 0.74–0.96), so a β on them largely restates the level models.

## Key columns

`rel_abundance_by_cell.tsv`, one row per (test_family, group, mag_id):

| Column | Meaning |
|---|---|
| `ra_pre`, `ra_end`, `ra_mean` | abundance level (%); `ra_mean` is the mean of the two |
| `delta_ra` | mean over cages of (END − PRE); signed |
| `abs_mean_delta_ra` | `abs(delta_ra)` — size of the net shift |
| `mean_abs_delta_ra` | mean of `abs(cage delta)` — size of the typical cage's shift |
| `delta_ra_contrast` | two-sample only: within-cage fat-minus-control difference-in-differences |
| `abs_mean_delta_ra_contrast` | two-sample only: `abs` of the above |
| `mean_abs_delta_ra_contrast` | two-sample only: mean of `abs(within-cage contrast)` |
| `n_replicates`, `n_mice` | 8, and 15 or 30 — assertion targets, constant by design |

`abs_mean_delta_ra` and `mean_abs_delta_ra` differ whenever cages move in opposite directions: the
first averages before taking the absolute value and so lets them cancel, the second does not.

## Assertions

`build_rel_abundance.py` raises rather than repairing. It checks that every Table S3 column sums to
100%, that every mouse has both timepoints and belongs to one cage and one arm, that every cage
appears in both arms, that the tested MAG set matches the pipeline's own eligibility flags, that
every eligible MAG is tested on all 8 replicates, that `delta_ra == ra_end − ra_pre`, that
`delta_ra_contrast == delta_ra_fat − delta_ra_control`, and that `ra_mean` equals a flat pool over the
(replicate, timepoint) values.

That last pair is the important one. `ra_mean` as a mean-of-means equals DRIDO's flat pool **only
while every cage appears at both timepoints**; the same condition is what makes the paired `delta_ra`
well defined. One assertion catches both.

The eligibility cross-check is scoped to `pre_end`, where the flags and the cell-stats MAG set agree
exactly (68 / 72 / 62). They do **not** agree in `pre_post` — see `DESIGN.md` section 3c(ii).

## If the run changes

The metadata-only shortcut here rests on `min_sample_num == 8 == the replicate count`, which forces
every tested position to carry all 8 cages. Lower the gate and positions with 5–7 cages become
testable, the tested cage set starts varying by MAG and by position, and abundance would have to be
derived per cell from the preprocessed files the way DRIDO does. `build_rel_abundance.py` asserts the
condition, so that change fails loudly instead of silently averaging over the wrong cages.
