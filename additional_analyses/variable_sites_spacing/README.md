# Variable sites and their distribution across the Fig 1 SGBs

Per-SGB and **per-contig** counts of variable sites, how far apart they sit, and
how that compares with BH-significant divergence and parallelism sites — comparing
the END and PRE timepoints, over the 62 SGBs tested in the divergence test.

Not part of the paper.

## What a variable site is

**A site that was TESTED** — a row in `p_value_summary`, significant or not.
The between-group and within-group tests preprocess separately and so test
different sites (21.5% of HF-tested sites are not divergence-tested), which means
the definition is not single-valued. Four are reported side by side, as a column
prefix:

| key | tested by | sites | contigs |
|---|---|---|---|
| `div` | `two_sample_paired_tTest` | 393,475 | 3,008 |
| `hf` | `single_sample_tTest`, `group_analyzed == 'fat'` | 332,435 | 3,507 |
| `lf` | `single_sample_tTest`, `group_analyzed == 'control'` | 353,873 | 3,949 |
| `union` | any of the three | 575,812 | 4,402 |

Each significance question is judged against its **own** tested set — Q5/Q6 vs
`div`, Q7/Q8 vs `hf`, Q9/Q10 vs `lf` — so a hit can never fall outside the
denominator its fraction is expressed against. `union` has no `sig` columns at
all: it is a union of tested sets, so it has no test and no significance question
of its own.

This replaces an earlier definition based on AlleleFlux's zero-difference filter.
That filter removes almost nothing (the resulting set covered 1–60% of each
genome), so it measured coverage breadth rather than variability.

## Layout

```
variable_site_distribution.py   analysis -> four output tables
plots.py                        figures -> 6 PNGs + a 62-page per-MAG PDF
METHODS.md                      how each metric is derived, worked on real data
```

`METHODS.md` and the column glossary below are maintained by hand; keep them in
step with `variable_site_distribution.py` when the columns or the arithmetic change.

## Running it

Read-only w.r.t. the AlleleFlux run; ~7 seconds, well under 4 GB.

```bash
OUT=/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/variable_site_stats

python3 variable_site_distribution.py --outdir $OUT     # --limit N for a smoke test
python3 plots.py                                        # defaults to $OUT and $OUT/figs
```

### Figures

`plots.py` writes six PNGs and `per_mag_cards.pdf` into `--outdir`. The PDF is one page
per SGB, split by definition: the spacing-vs-null scatter and the gap ECDF, each as four
**contig-level** panels (`div` / `hf` / `lf` / `union`), beside that SGB's summary numbers.
Output goes to `$OUT/figs` on scratch rather than the repo; the per-MAG PDF alone is
~360 KB and regenerates in 25 s, so it is output, not source. It reads only the tables
and computes no new statistics, so a figure can never disagree with them. Colours are
the Okabe-Ito colourblind-safe set; `union` is deliberately neutral grey because it is
an aggregate of the other three, not a peer category.

## Outputs

| file | rows x cols | one row is | derived from |
|---|---|---|---|
| `contigs/<MAG>.tsv` | 62 files | one contig of that SGB | the run |
| `contig_level_all.tsv` | 4,402 x 38 | one contig, any SGB | the per-MAG files, concatenated |
| `sgb_level.tsv` | 62 x 43 | one SGB | that SGB's contig table |
| `summary.tsv` | 4 x 13 | one variable-site definition | the SGB table |

Each level is built from the level below it, never recomputed from the raw sites,
so the three cannot disagree. Every mean spacing ships with its `n_gaps` support
count; every significant-contig count ships with its `pct_contigs_sig` fraction.

## Column glossary

Throughout, **"sites" on its own means TESTED sites**. Significant ones always
carry `sig` in the name. `{d}` is one of `div`, `hf`, `lf`, `union`. Worked
examples use SGB **`SLG191_DASTool_bins_89`** and, within it, contig
**`k141_116087`** (89,547 bp).

### Contig level: `contigs/<MAG>.tsv` and `contig_level_all.tsv`

Same columns; the combined file adds only `MAG_ID`. There is no "number of
contigs" column, the row count *is* that number.

| column | meaning | how it is calculated |
|---|---|---|
| `MAG_ID` | which SGB the contig belongs to | from the MAG-to-contig mapping |
| `contig` | contig name, `<MAG>.fa_k141_<n>` | as written by the assembler |
| `contig_len` | length in bp | column 2 of the `.fai`; `-1` if absent from it |
| `{d}_n_sites` | tested sites on this contig | count of `p_value_summary` rows |
| `{d}_mean_gap` | mean spacing between them (Q3) | `np.diff(sorted(positions)).mean()` |
| `{d}_median_gap` | median spacing | `np.median` of the same gaps |
| `{d}_n_gaps` | gaps behind those means | `n_sites - 1`, or 0 when `n_sites < 2` |
| `{d}_expected_gap` | spacing if the same number of sites were scattered uniformly (Q4) | `(contig_len + 1) / (n_sites + 1)` |
| `{d}_n_sig` | significant sites on this contig | those with `q_value < 0.05` |
| `{d}_has_sig` | whether it carries any | `n_sig > 0`; the SGB level counts contigs by this |
| `{d}_sig_mean_gap` / `{d}_sig_median_gap` | spacing between significant sites (Q6/Q8/Q10) | same `np.diff` on the significant positions |
| `{d}_sig_n_gaps` | gaps behind them | `n_sig - 1`, or 0 |

Worked on `k141_116087`:

```
contig_len          = 89,547
div_n_sites        =    164   mean_gap =     23.69   n_gaps =    163
hf_n_sites         =    169   mean_gap =     10.34   n_gaps =    168
lf_n_sites         =    196   mean_gap =    318.32   n_gaps =    195
union_n_sites      =    264   mean_gap =    236.02   n_gaps =    263
```

`div_expected_gap` = (89,547 + 1) / (164 + 1) = **542.72**, against an observed
23.69: the tested sites on this contig are about 23x more tightly packed than
random placement predicts. It is `L + 1` over `n + 1`, not `L / n`: sites never
reach the contig ends, so the span they cover is under `L`. `METHODS.md` derives
this and checks it against exhaustive enumeration.

**Zeros and NaNs are not the same.** Of this SGB's 69 contigs, 66 have
`div_n_sites = 0` (divergence never tested them) and their `div_mean_gap` is
`NaN`, meaning no gap could be measured, not that the gap is zero. A contig with
exactly one site is also `NaN`, which is why every mean ships beside its `n_gaps`.

### SGB level: `sgb_level.tsv`

One row per SGB. Worked on `SLG191_DASTool_bins_89`, which has **102 reference
contigs** and 69 contig rows.

| column | meaning | how it is calculated |
|---|---|---|
| `MAG_ID` | the SGB | |
| `species` | GTDB species label | from the scores table |
| `n_ref_contigs` | **contigs in the reference** (Q2) | counted in the MAG-to-contig mapping. **Not** from the data, so it includes contigs no test reached |
| `ref_genome_len` | summed length of those | from the `.fai` |
| `{d}_n_contigs_with_sites` | **contigs carrying >=1 TESTED site** (Q1) | `(contig_table.{d}_n_sites > 0).sum()` |
| `{d}_n_sites` | tested sites | `contig_table.{d}_n_sites.sum()` |
| `{d}_mean_gap` | pooled spacing | gap-count-weighted: `sum(mean_gap * n_gaps) / sum(n_gaps)`, not a mean of contig means |
| `{d}_expected_gap` | pooled null | same weights, so the two stay comparable |
| `{d}_n_gaps` | total gaps | `contig_table.{d}_n_gaps.sum()` |
| `{d}_pct_genome_tested` | tested sites per bp of reference | `100 * n_sites / ref_genome_len` |
| `{d}_n_contigs_sig` | **contigs carrying >=1 SIGNIFICANT site** (Q5/Q7/Q9) | `contig_table.{d}_has_sig.sum()` |
| `{d}_pct_contigs_sig` | that, as a share of this test's own tested contigs | `100 * n_contigs_sig / n_contigs_with_sites` |
| `{d}_n_sig_sites` | significant sites | `contig_table.{d}_n_sig.sum()` |
| `{d}_sig_mean_gap` | pooled spacing between significant sites | weighted by `sig_n_gaps` |
| `{d}_sig_n_gaps` | gaps behind that | |

**The pooled mean is not an average of averages:**

```
contigs contributing gaps      : 3
mean of the per-contig means   : 40.4800   <- WRONG, each contig counts equally
gap-count-weighted (as coded)  : 25.1094   <- what SLG191_DASTool_bins_89 reports
total gap distance / total gaps: 4,469 / 178 = 25.1094
```

The last two lines agree exactly: pooling reconstructs the answer you would get
from the raw gaps, while averaging the per-contig means lets a 4-gap contig weigh
as much as a 163-gap one.

**The three contig counts nest**, `n_ref_contigs` >= `{d}_n_contigs_with_sites`
>= `{d}_n_contigs_sig`:

| | reference | with tested sites | with significant sites | % |
|---|---|---|---|---|
| `div` | 102 | 3 | 2 | 66.7% |
| `hf` | 102 | 4 | 2 | 50.0% |
| `lf` | 102 | 69 | 0 | 0.0% |
| `union` | 102 | 69 | — | — |

So `div_pct_contigs_sig` = 100 x 2 / 3 = **66.67%**. The denominator is
divergence's own tested contigs, not the reference count and not the contigs
some other test reached.

### Summary: `summary.tsv`

One row per definition, pooled over every SGB. Same meanings as the SGB level
with the `{d}_` prefix dropped, since the definition is now a value in the
`definition` column. Two additions: `definition`, and `n_SGBs` (SGBs with at
least one tested site under it). `n_ref_contigs` is identical on every row.

**Percentages are recomputed, not averaged:** as coded, 100 * 1,155 / 3,008 =
38.40%; averaging the 62 per-SGB percentages instead gives 34.31%. The same
weighting trap as the means.

**Zero and blank mean different things:**

```
definition  n_contigs_sig  pct_contigs_sig  n_sig_sites  sig_mean_gap
       div         1155.0            38.40      62490.0        161.74
        hf         1948.0            55.55      90904.0        234.55
        lf            0.0             0.00          0.0           NaN
     union            NaN              NaN          NaN           NaN
```

`lf` shows **0**, a measured absence: control-group sites were tested and none
reached `q < 0.05`. `union` shows **blank**: the question does not exist for it.
Filling union with 0 would falsely claim it was measured and came back empty.

## Headline results

Pooled over all 62 SGBs (`summary.tsv`):

| definition | contigs w/ sites | sites | mean gap | expected gap | ratio |
|---|---|---|---|---|---|
| `div` | 3,008 | 393,475 | 160 | 232 | 0.69 |
| `hf` | 3,507 | 332,435 | 207 | 308 | 0.67 |
| `lf` | 3,949 | 353,873 | 213 | 316 | 0.67 |
| `union` | 4,402 | 575,812 | 151 | 214 | 0.70 |

Significance, as a fraction of each test's own tested contigs:

| | overall | min | median | max |
|---|---|---|---|---|
| divergence | 38.4% | 0.0% | 35.6% | 100.0% |
| HF parallelism | 55.5% | 0.0% | 41.4% | 96.0% |
| LF parallelism | 0.0% | 0.0% | 0.0% | 0.0% |

Observed spacing is consistently **tighter** than the uniform null — the ratio
runs 0.001–0.999 (median 0.314) for `div`, i.e. tested sites are clustered
roughly threefold. That comparison only became informative once a variable site
was defined as a tested site; under the old dense definition the observed mean
and the null coincide arithmetically.

**LF parallelism is empty.** No control-group site reaches q < 0.05 anywhere; the
minimum q genome-wide is 0.195. Q9 is 0% for every SGB and Q10 undefined. LF is
not short of tested sites (353,873, more than divergence's 393,475): this is an
absence of signal, not of data.

## Known caveats

1. **FDR is pooled across all SGBs**, not within each one, so a per-SGB q<0.05
   count is not a per-SGB FDR. Recomputing per SGB moves q by up to 0.83.
2. **HF = `fat`, LF = `control`** is an assumption. The Fig 1 code only says
   `Fat`/`Control`; the sole HF/LF mapping is in
   `figures/Fig2/per_litter_violins.ipynb`.
3. **The uniform null is `(L+1)/(n+1)`, not `L/n`.** Sites never reach the contig
   ends, so the span they cover is under L. `METHODS.md` derives this and checks
   it against exhaustive enumeration. Do not simplify it away.
4. **`p_value_summary` needs `test_type` AND `group_analyzed` pinned.** The code
   asserts uniqueness rather than deduping; a dedupe would silently merge the two
   diet groups. The assertion was tested by forcing the collision (203,926 dupes).
5. **`position` is 0-based** and position 0 occurs, so coordinates occupy
   `[0, L-1]`; gap arithmetic is unaffected, anything compared against contig
   length is not. `min_p_value` is a minimum over four nucleotides with no
   selection correction, which inflates counts but not spacing geometry.
