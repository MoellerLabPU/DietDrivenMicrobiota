# Relative abundance vs. AlleleFlux significance — what every term means

This note explains every quantity in `rel_abundance_regression.ipynb` and how it is calculated,
using real numbers from the study throughout. The analysis asks: *is the significance AlleleFlux
reports for a MAG explained by how abundant that MAG is, or by how much its abundance changed?*
If yes, the evolutionary signal might be an artifact of sequencing depth; if no (or if the
association runs the "wrong" way), the signal stands on its own.

**Scope.** The `AlleleFlux_mapq20` run, PRE vs END, high-fat (`fat`) vs control diet. The analysis
is restricted to the **62 MAGs the divergence test actually tested** — every regression, including
the parallelism ones, runs on this same MAG set.

---

## 1. The three units of analysis

Each MAG gets up to three rows ("cells"), one per statistical test AlleleFlux ran on it:

| Unit | AlleleFlux test | What it asks |
|---|---|---|
| parallelism: fat | `single_sample_tTest` on the fat arm | did allele frequencies shift *consistently across the fat-arm cages*? |
| parallelism: control | `single_sample_tTest` on the control arm | same, within the control arm |
| divergence (fat vs control) | `two_sample_paired_tTest` | did the two arms' allele frequencies *move apart*, pairing arms within each cage? |

**The unit of observation is the cage (replicate), never the mouse.** AlleleFlux averages all mice
sharing a cage into one value before testing, so every test has n = 8 cages. The abundance numbers
below are aggregated the same way on purpose — otherwise the abundance would describe a different
population than the p-value does (cage 2 has one mouse per arm; the others have two).

## 2. How the abundance numbers are built — the three-step ladder

All abundance comes from the published Table S3: per-sample relative abundance in %, summing to
100% over all 160 genomes (not renormalized to the tested subset). Every derived quantity climbs
the same ladder: **mouse → cage → MAG**.

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
> absolute value can appear**, and it appears in two deliberately different flavors (§3).

### 2.1 The worked example, end to end

All of §3 computes from one real MAG, `SLG191_DASTool_bins_89`. Step 1 first, shown for one cage:
in cage 4, control mouse 542 went 0.0634 → 0.0000 (Δ = −0.0634) and mouse 543 went
0.1122 → 4.0890 (Δ = +3.9768), so step 2 gives

d(cage 4, control) = (−0.0634 + 3.9768) / 2 = **+1.9567**.

Fat mice in the same cage: Δ = −0.4948 and −4.8235, so d(cage 4, fat) = (−0.4948 − 4.8235)/2 =
**−2.6592**. Plain means, signs kept. Doing this for every cage gives the step-2 table every
column in §3 is computed from (RA in %; c is defined in §3):

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
rounded column can differ in the last digit.)

(This MAG collapsed to 0% in every fat-arm mouse by END — which is why its fat `end` column is
all zeros.)

## 3. Glossary of the abundance columns, with the calculations

All values in percentage points of relative abundance. Each formula is followed by the worked
number from the §2.1 table; every result matches the shipped `rel_abundance_by_cell.tsv` exactly.

### Level — how abundant is the MAG?

| Column | Formula (over the 8 cage values) | Reads as |
|---|---|---|
| `ra_pre` | mean of the cage-mean PRE levels | typical abundance before the diet switch |
| `ra_end` | same, at END | typical abundance at the end |
| `ra_mean` | (`ra_pre` + `ra_end`) / 2 | overall "how big is this organism" |

Worked, control arm:

- `ra_pre` = 17.6134 / 8 = **2.2017**
- `ra_end` = 17.6858 / 8 = **2.2107**
- `ra_mean` = (2.2017 + 2.2107) / 2 = **2.2062**

For the divergence unit the same formulas run over all 16 (cage × arm) values pooled:
`ra_pre` = (17.6134 + 8.9281) / 16 = **1.6588**, `ra_end` = (17.6858 + 0) / 16 = **1.1054**.

### Change — how much did it move?

| Column | Formula (d₁…d₈ are the signed cage deltas) | Reads as |
|---|---|---|
| `delta_ra` | (d₁ + … + d₈) / 8 | net shift, direction kept (equals `ra_end − ra_pre` exactly) |
| `abs_mean_delta_ra` | \| (d₁ + … + d₈) / 8 \| — **average first, then absolute** | size of the *net* shift; opposite-moving cages cancel before the absolute value |
| `mean_abs_delta_ra` | ( \|d₁\| + … + \|d₈\| ) / 8 — **absolute first, then average** | how far a *typical* cage moved, no cancelling |

Worked, control arm (deltas −0.6513, −0.5658, +1.9567, +1.1629, −0.5875, +0.8645, −2.0275, −0.0797):

- `delta_ra` = +0.0723 / 8 = **+0.0090** — the rises and falls almost perfectly cancel
- `abs_mean_delta_ra` = |+0.0090| = **0.0090**
- `mean_abs_delta_ra` = (0.6513 + 0.5658 + 1.9567 + 1.1629 + 0.5875 + 0.8645 + 2.0275 + 0.0797) / 8
  = 7.8958 / 8 = **0.9870**

That 0.009-vs-0.987 gap is the whole reason both variants exist: net, this MAG barely moved in the
control arm; per cage, it typically moved by nearly a full percentage point — just in different
directions. In the fat arm every cage fell, so the three numbers coincide:
`delta_ra` = −8.9281 / 8 = **−1.1160**, and both absolute variants = **1.1160**.
Divergence (pooled over 16): `delta_ra` = (+0.0723 − 8.9281) / 16 = **−0.5535**;
`mean_abs_delta_ra` = (7.8959 + 8.9281) / 16 = **1.0515**.

### Contrast — did fat move *differently* from control? (divergence unit only)

The divergence test never asks "did the MAG change"; it asks, within each cage, "did the fat side
change differently from the control side". The contrast is the abundance quantity with that exact
shape — a difference-in-differences, again computed per cage first:

c(cage) = d(cage, fat) − d(cage, control)

e.g. cage 4: c = −2.6592 − (+1.9567) = **−4.6159** (fat crashed while control rose — a large
diet-driven disagreement in that cage).

| Column | Formula (c₁…c₈ from the table) | Reads as |
|---|---|---|
| `delta_ra_fat`, `delta_ra_control` | mean of each arm's cage deltas | the two arms' net shifts, separately |
| `delta_ra_contrast` | (c₁ + … + c₈) / 8 | how much *more* fat moved than control, net |
| `abs_mean_delta_ra_contrast` | \| (c₁ + … + c₈) / 8 \| | size of that net difference |
| `mean_abs_delta_ra_contrast` | ( \|c₁\| + … + \|c₈\| ) / 8 | typical cage's fat-vs-control gap, no cancelling |

Worked:

- `delta_ra_fat` = **−1.1160**, `delta_ra_control` = **+0.0090** (from the change section)
- `delta_ra_contrast` = −9.0004 / 8 = **−1.1251**
  (equivalently −1.1160 − (+0.0090) — the two routes agree exactly)
- `abs_mean_delta_ra_contrast` = |−1.1251| = **1.1251**
- `mean_abs_delta_ra_contrast` = (0.5038 + 0.0407 + 4.6159 + 2.1661 + 0.6046 + 2.7664 + 1.5359
  + 0.0797) / 8 = 12.3131 / 8 = **1.5391**

Why the contrast matters: a MAG that doubled in every mouse *regardless of diet* has a large
`delta_ra` but a contrast near zero — exactly the MAG the divergence test should not flag. This
MAG is the opposite case: it crashed only under fat, so its contrast (−1.13) is almost entirely
diet-driven change.

## 4. The significance side (responses)

Each cell also carries the AlleleFlux test results for that MAG, summarized three ways:

| Term | Formula | Reads as |
|---|---|---|
| `-log10(min p-value)` | −log₁₀( smallest per-site p in the MAG ) | strength of the single best site (bigger = more significant) |
| `-log10(min q-value)` | −log₁₀( smallest BH/FDR-adjusted q ) | best site after multiple-testing correction (q is corrected genome-wide across all MAGs, not per MAG) |
| breadth: % of sites p<0.05 (`pct_sig_p`) | 100 × (sites with p < 0.05) / (sites tested in the MAG) | how *widespread* the signal is across the genome, not just its best site |

Worked, `SLG659_DASTool_bins_22`, divergence unit: 19 sites tested, 4 with p < 0.05, best site
p = 0.02467, best q = 0.11121:

- breadth = 100 × 4 / 19 = **21.05**
- −log₁₀(0.02467) = **1.61** (a p of 0.001 would give 3; each unit is a factor of 10)
- −log₁₀(0.11121) = **0.95**

The −log₁₀ transform is only a change of units so that "more significant" plots upward and tiny
p-values don't all squash against zero.

## 5. The regression terms

Every model is **univariate**: one abundance predictor x against one significance response y,
62 points (one per MAG), fit by ordinary least squares (OLS):

y = β₀ + β₁·z,  where z = (x − mean(x)) / SD(x)

- **z-scoring**: each predictor is centered and divided by its standard deviation before fitting,
  so **β₁ is in response units per 1 SD of the predictor** and βs are comparable across predictors.
  (This rescales β and its standard error only; p and R² are identical to the unscaled fit.)
  Example: β = −6.1 for "breadth vs PRE abundance" means one SD higher log-abundance predicts 6.1
  fewer percentage points of significant sites.
- **log10 on level predictors** — abundance spans ~70-fold with a heavy right tail, so levels are
  log10-transformed before z-scoring; equal fold-changes become equal distances and the few very
  abundant MAGs can't dominate the slope. Change/contrast predictors stay linear (they can be ≤ 0).
- **p** — tests β₁ = 0; reported raw, with no multiple-testing correction across the 63 fits
  (stated in the notebook, so ~3 of 63 are expected significant by chance).
- **R²** — fraction of the response's variance the predictor explains (0 = none, 1 = all).
- **Spearman ρ / Spearman p** — the same points, but correlated on **ranks**: both variables are
  replaced by their rank order (1st, 2nd, …) and the correlation is computed on those. Catches any
  monotone trend and is immune to the extreme-abundance MAGs (ranks don't care by *how much* the
  top MAG leads). OLS and Spearman agreeing is the credibility check.

## 6. Direct answer to the cage-mean question

**The mean of each cage is always a normal (signed, unweighted) mean — of the mouse-level values
within that cage and arm. No absolute value is involved at the cage step, for any column.**
Absolute values enter only at the final across-cages step, and only in the columns whose names say
so, in two flavors: `abs_mean_*` = absolute value *of* the mean (cancellation allowed, then |·|),
and `mean_abs_*` = mean of the absolute values (no cancellation). The plain `delta_ra` and
`delta_ra_contrast` never touch an absolute value at any step. §3's control-arm example is the
two flavors at their most different: 0.009 vs 0.987 from the same eight cage values.

---

*Code: `build_rel_abundance.py` (steps 1–3 and all integrity checks) and
`rel_abundance_regression.ipynb` (models and figures) in `Revision/relative_abundance/`;
method rationale in `DESIGN.md`. Every number in this note was recomputed from Table S3 and
matches the shipped `rel_abundance_by_cell.tsv`.*
