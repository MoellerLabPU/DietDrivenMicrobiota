# Relative abundance vs. AlleleFlux significance

Is the significance AlleleFlux reports for a MAG associated with **how abundant** that MAG is, or with
**how much its abundance changed** between PRE and END? If yes, the evolutionary signal might be an
artifact of sequencing depth; if no, or if the association runs the "wrong" way, the signal stands
on its own.

**Scope.** The `pre_end-fat_control` comparison of the `AlleleFlux_mapq20` run, `single_sample_tTest`
(parallelism, within one arm) and `two_sample_paired_tTest` (divergence, arms paired by cage). Every
regression, including the parallelism ones, is restricted to the **62 MAGs the divergence test
actually tested**. `rel_abundance_by_cell.tsv` itself still carries all 202 tested cells; the subset
happens in the notebook, after the merge.

## Layout

```
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
alleleflux-significant-sites-summary \
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

## The three units of analysis

Each MAG gets up to three rows ("cells"), one per statistical test AlleleFlux ran on it:

| Unit | AlleleFlux test | What it asks |
|---|---|---|
| parallelism: fat | `single_sample_tTest` on the fat arm | did allele frequencies shift *consistently across the fat-arm cages*? |
| parallelism: control | `single_sample_tTest` on the control arm | same, within the control arm |
| divergence (fat vs control) | `two_sample_paired_tTest` | did the two arms' allele frequencies *move apart*, pairing arms within each cage? |

**The unit of observation is the cage (replicate), never the mouse.** AlleleFlux averages every mouse
sharing a `replicate` into one value before testing, so each test here has *n = 8*, not 15 or 30.
Abundance is aggregated the same way — per mouse, then into cages, then across cages. Averaging over
mice instead would weight the one-mouse cage differently and describe a different population than
the p-value does (cage 2 has one mouse per arm; the others have two). This is the single most
important thing to preserve if the code is adapted.

**Relative abundance is the published Table S3, unchanged** — a per-sample composition over all 160
genomes, summing to 100%. It is *not* renormalized to the subset AlleleFlux tested.

**Every model is univariate.** Abundance and change are tested separately, never in one formula, so
each β is a total association. The absolute change/contrast predictors themselves track abundance
level (Spearman ρ 0.74–0.96), so a β on them largely restates the level models.

## How the abundance numbers are built — the three-step ladder

Every derived quantity climbs the same ladder: **mouse → cage → MAG**.

> **Step 1 — per mouse:**
> `Δ(mouse) = RA(END) − RA(PRE)`. Each mouse is its own baseline. Signed.
>
> **Step 2 — per cage and arm:**
> `d(cage, arm) = mean of Δ(mouse) over that arm's mice in that cage` — a **plain, unweighted
> arithmetic mean**. **No absolute value is taken at this step — ever.** A cage where one mouse
> rose +2 and the other fell −2 gets a cage value of 0. Unweighted means a 1-mouse cage counts
> exactly as much as a 2-mouse cage, mirroring what AlleleFlux itself does before testing.
> (The same plain mean also produces per-cage PRE and END levels.)
>
> **Step 3 — per MAG:**
> the 8 cage values are combined into the reported columns. **This is the only step where an
> absolute value can appear**, and it appears in two deliberately different flavors (see the
> glossary).

### The worked example, end to end

The glossary below computes from one real MAG, `SLG191_DASTool_bins_89`. Step 1 first, shown for one
cage: in cage 4, control mouse 542 went 0.0634 → 0.0000 (Δ = −0.0634) and mouse 543 went
0.1122 → 4.0890 (Δ = +3.9768), so step 2 gives

d(cage 4, control) = (−0.0634 + 3.9768) / 2 = **+1.9567**.

Fat mice in the same cage: Δ = −0.4948 and −4.8235, so d(cage 4, fat) = (−0.4948 − 4.8235)/2 =
**−2.6592**. Plain means, signs kept. Doing this for every cage gives the step-2 table every
column below is computed from (RA in %; c is defined under Contrast):

| cage | pre (ctrl) | end (ctrl) | **d (ctrl)** | pre (fat) | end (fat) | **d (fat)** | **c = d(fat) − d(ctrl)** |
|---|---|---|---|---|---|---|---|
| 1 | 3.2925 | 2.6412 | −0.6513 | 1.1551 | 0.0000 | −1.1551 | −0.5038 |
| 2 | 3.2917 | 2.7260 | −0.5658 | 0.5250 | 0.0000 | −0.5250 | +0.0407 |
| 4 | 0.0878 | 2.0445 | +1.9567 | 2.6592 | 0.0000 | −2.6592 | −4.6159 |
| 5 | 1.0164 | 2.1792 | +1.1629 | 1.0032 | 0.0000 | −1.0032 | −2.1661 |
| 6 | 0.7194 | 0.1319 | −0.5875 | 1.1921 | 0.0000 | −1.1921 | −0.6046 |
| 7 | 2.5410 | 3.4055 | +0.8645 | 1.9019 | 0.0000 | −1.9019 | −2.7664 |
| 8 | 2.8877 | 0.8602 | −2.0275 | 0.4916 | 0.0000 | −0.4916 | +1.5359 |
| 9 | 3.7769 | 3.6973 | −0.0797 | 0.0000 | 0.0000 | 0.0000 | +0.0797 |
| **sum** | **17.6134** | **17.6858** | **+0.0723** | **8.9281** | **0.0000** | **−8.9281** | **−9.0004** |

(Cage values are shown to 4 decimals; the sums are computed from full precision, so adding the
rounded column can differ in the last digit. This MAG collapsed to 0% in every fat-arm mouse by
END, which is why its fat `end` column is all zeros.)

## Glossary of the abundance columns, with the calculations

`rel_abundance_by_cell.tsv` has one row per (test_family, group, mag_id). All values are in
percentage points of relative abundance. Each formula is followed by the worked number from the
table above; every result matches the shipped table exactly. `n_replicates` (8) and `n_mice`
(15 or 30) are assertion targets, constant by design.

### Level — how abundant is the MAG?

| Column | Formula (over the 8 cage values) | Reads as |
|---|---|---|
| `ra_pre` | mean of the cage-mean PRE levels | typical abundance before the diet switch |
| `ra_end` | same, at END | typical abundance at the end |
| `ra_mean` | (`ra_pre` + `ra_end`) / 2 | overall "how big is this organism" |

Worked, control arm: `ra_pre` = 17.6134 / 8 = **2.2017**; `ra_end` = 17.6858 / 8 = **2.2107**;
`ra_mean` = **2.2062**. For the divergence unit the same formulas run over all 16 (cage × arm)
values pooled: `ra_pre` = (17.6134 + 8.9281) / 16 = **1.6588**, `ra_end` = **1.1054**.

### Change — how much did it move?

| Column | Formula (d₁…d₈ are the signed cage deltas) | Reads as |
|---|---|---|
| `delta_ra` | (d₁ + … + d₈) / 8 | net shift, direction kept (equals `ra_end − ra_pre` exactly) |
| `abs_mean_delta_ra` | \| (d₁ + … + d₈) / 8 \| — **average first, then absolute** | size of the *net* shift; opposite-moving cages cancel before the absolute value |
| `mean_abs_delta_ra` | ( \|d₁\| + … + \|d₈\| ) / 8 — **absolute first, then average** | how far a *typical* cage moved, no cancelling |

Worked, control arm (deltas −0.6513, −0.5658, +1.9567, +1.1629, −0.5875, +0.8645, −2.0275, −0.0797):

- `delta_ra` = +0.0723 / 8 = **+0.0090** — the rises and falls almost perfectly cancel
- `abs_mean_delta_ra` = |+0.0090| = **0.0090**
- `mean_abs_delta_ra` = 7.8958 / 8 = **0.9870**

That 0.009-vs-0.987 gap is the whole reason both variants exist: net, this MAG barely moved in the
control arm; per cage, it typically moved by nearly a full percentage point, just in different
directions. In the fat arm every cage fell, so the three numbers coincide: `delta_ra` = −8.9281 / 8
= **−1.1160**, and both absolute variants = **1.1160**. Divergence (pooled over 16): `delta_ra` =
**−0.5535**; `mean_abs_delta_ra` = **1.0515**.

### Contrast — did fat move *differently* from control? (divergence unit only)

The divergence test never asks "did the MAG change"; it asks, within each cage, "did the fat side
change differently from the control side". The contrast is the abundance quantity with that exact
shape, a difference-in-differences, again computed per cage first: c(cage) = d(cage, fat) −
d(cage, control). E.g. cage 4: c = −2.6592 − (+1.9567) = **−4.6159**.

| Column | Formula (c₁…c₈ from the table) | Reads as |
|---|---|---|
| `delta_ra_fat`, `delta_ra_control` | mean of each arm's cage deltas | the two arms' net shifts, separately |
| `delta_ra_contrast` | (c₁ + … + c₈) / 8 | how much *more* fat moved than control, net |
| `abs_mean_delta_ra_contrast` | \| (c₁ + … + c₈) / 8 \| | size of that net difference |
| `mean_abs_delta_ra_contrast` | ( \|c₁\| + … + \|c₈\| ) / 8 | typical cage's fat-vs-control gap, no cancelling |

Worked: `delta_ra_fat` = **−1.1160**, `delta_ra_control` = **+0.0090**; `delta_ra_contrast` =
−9.0004 / 8 = **−1.1251** (equivalently −1.1160 − 0.0090); `abs_mean_delta_ra_contrast` =
**1.1251**; `mean_abs_delta_ra_contrast` = 12.3131 / 8 = **1.5391**.

Why the contrast matters: a MAG that doubled in every mouse *regardless of diet* has a large
`delta_ra` but a contrast near zero, exactly the MAG the divergence test should not flag. This MAG
is the opposite case: it crashed only under fat, so its contrast is almost entirely diet-driven.

### Where the absolute values are, in one sentence

The mean of each cage is always a plain signed mean of the mouse-level values within that cage and
arm. Absolute values enter only at the final across-cages step, and only in the columns whose names
say so: `abs_mean_*` = absolute value *of* the mean (cancellation allowed, then |·|), `mean_abs_*` =
mean of the absolute values (no cancellation). `delta_ra` and `delta_ra_contrast` never touch one.

## The significance side (responses)

Each cell also carries the AlleleFlux test results for that MAG, summarized three ways:

| Term | Formula | Reads as |
|---|---|---|
| `-log10(min p-value)` | −log₁₀( smallest per-site p in the MAG ) | strength of the single best site |
| `-log10(min q-value)` | −log₁₀( smallest BH/FDR-adjusted q ) | best site after multiple-testing correction (q is corrected genome-wide across all MAGs, not per MAG) |
| breadth: % of sites p<0.05 (`pct_sig_p`) | 100 × (sites with p < 0.05) / (sites tested in the MAG) | how *widespread* the signal is across the genome |

Worked, `SLG659_DASTool_bins_22`, divergence unit: 19 sites tested, 4 with p < 0.05, best site
p = 0.02467, best q = 0.11121: breadth = **21.05**, −log₁₀(p) = **1.61**, −log₁₀(q) = **0.95**.

## The regression terms

Every model is **univariate**: one abundance predictor x against one significance response y,
62 points (one per MAG), fit by ordinary least squares:

y = β₀ + β₁·z,  where z = (x − mean(x)) / SD(x)

- **z-scoring**: each predictor is centered and scaled before fitting, so **β₁ is in response units
  per 1 SD of the predictor** and βs are comparable across predictors. This rescales β and its
  standard error only; p and R² are identical to the unscaled fit.
- **log10 on level predictors**: abundance spans ~70-fold with a heavy right tail, so levels are
  log10-transformed before z-scoring. Change/contrast predictors stay linear (they can be ≤ 0).
- **p** tests β₁ = 0 and is reported raw, with no multiple-testing correction across the 63 fits
  (~3 of 63 are expected significant by chance).
- **R²** is the fraction of the response's variance the predictor explains.
- **Spearman ρ / p**: the same points correlated on ranks. Catches any monotone trend and is immune
  to the extreme-abundance MAGs. OLS and Spearman agreeing is the credibility check.

## Assertions

`build_rel_abundance.py` raises rather than repairing. It checks that every Table S3 column sums to
100%, that every mouse has both timepoints and belongs to one cage and one arm, that every cage
appears in both arms, that the tested MAG set matches the pipeline's own eligibility flags, that
every eligible MAG is tested on all 8 replicates, that `delta_ra == ra_end − ra_pre`, that
`delta_ra_contrast == delta_ra_fat − delta_ra_control`, and that `ra_mean` equals a flat pool over the
(replicate, timepoint) values.

That last pair is the important one. `ra_mean` as a mean-of-means equals a flat pool **only while
every cage appears at both timepoints**; the same condition is what makes the paired `delta_ra`
well defined. One assertion catches both.

The eligibility cross-check is scoped to `pre_end`, where the flags and the cell-stats MAG set agree
exactly (68 / 72 / 62). They do **not** agree in `pre_post`.

## If the run changes

The metadata-only shortcut here rests on `min_sample_num == 8 == the replicate count`, which forces
every tested position to carry all 8 cages. Lower the gate and positions with 5–7 cages become
testable, the tested cage set starts varying by MAG and by position, and abundance would have to be
derived per cell from the preprocessed files. `build_rel_abundance.py` asserts the condition, so
that change fails loudly instead of silently averaging over the wrong cages.
