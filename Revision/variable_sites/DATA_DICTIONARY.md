# Data dictionary

What each file holds, what every column means, and how it was calculated —
with the arithmetic worked through on real rows.

The numbers below were read from the tables directly. The generator that produced
this file has been removed, so it is now maintained by hand — update it alongside
any change to `variable_site_distribution.py`. For the derivations behind the
metrics (why the uniform null is `(L+1)/(n+1)`, why gaps are within-contig only)
see `METHODS.md`.

Worked examples use SGB **`SLG191_DASTool_bins_89`** and, within it, contig **`k141_116087`**.

---

## The files

| file | rows x cols | one row is | derived from |
|---|---|---|---|
| `contigs/<MAG>.tsv` | 62 files | one contig of that SGB | the run |
| `contig_level_all.tsv` | 4,402 x 38 | one contig, any SGB | the per-MAG files, concatenated |
| `sgb_level.tsv` | 62 x 43 | one SGB | that SGB's contig table |
| `summary.tsv` | 4 x 13 | one variable-site definition | the SGB table |

Each level is built from the level below it, never recomputed from the raw sites.
That is what guarantees the three cannot disagree: an SGB row is an aggregation
of exactly the contig rows written to disk beside it.

`figs/` holds the plots, written by `plots.py` from these tables. It contains no
numbers of its own.

## The four definitions

A **variable site is a site that was TESTED** — a row in `p_value_summary`,
significant or not. The between-group and within-group tests preprocess
separately and therefore test different sites, so the definition is not
single-valued. Every metric is reported under all four, as a column prefix:

| prefix | tested by | sites (all SGBs) | contigs |
|---|---|---|---|
| `div_` | `two_sample_paired_tTest` | 393,475 | 3,008 |
| `hf_` | `single_sample_tTest`, `group_analyzed == 'fat'` | 332,435 | 3,507 |
| `lf_` | `single_sample_tTest`, `group_analyzed == 'control'` | 353,873 | 3,949 |
| `union_` | any of the three above | 575,812 | 4,402 |

**Significance is always judged against the matching tested set** — `div_n_sig`
counts divergence-significant sites among divergence-tested ones. A significant
site can therefore never fall outside the denominator it is expressed against.
`union` has no `sig` columns at all: it is a union of tested sets, so it has no
test and no significance question of its own.

Throughout, **"sites" unqualified means tested sites**. Significant ones always
carry `sig` in the column name.

## Contig-level columns

`contigs/<MAG>.tsv` and `contig_level_all.tsv` carry the same columns; the
combined file adds nothing but the `MAG_ID` needed to tell the rows apart.
There is no "number of contigs" column here — the row count *is* that number.

Worked on contig `k141_116087` of `SLG191_DASTool_bins_89`, length **89,547 bp**.

| column | meaning | how it is calculated |
|---|---|---|
| `MAG_ID` | which SGB the contig belongs to | from the MAG-to-contig mapping |
| `contig` | contig name, `<MAG>.fa_k141_<n>` | as written by the assembler |
| `contig_len` | length in bp | column 2 of the `.fai`; `-1` if absent from it |
| `{d}_n_sites` | tested sites on this contig | count of `p_value_summary` rows |
| `{d}_mean_gap` | mean spacing between them | `np.diff(sorted(positions)).mean()` |
| `{d}_median_gap` | median spacing | `np.median` of the same gaps |
| `{d}_n_gaps` | gaps behind those means | `n_sites - 1`, or 0 when `n_sites < 2` |
| `{d}_expected_gap` | spacing if the same number of sites were scattered uniformly | `(contig_len + 1) / (n_sites + 1)` |
| `{d}_n_sig` | significant sites on this contig | those with `q_value < 0.05` |
| `{d}_has_sig` | whether it carries any | `n_sig > 0`; the SGB level counts contigs by this |
| `{d}_sig_mean_gap` | mean spacing between significant sites | same `np.diff` on the significant positions |
| `{d}_sig_median_gap` | median of those | |
| `{d}_sig_n_gaps` | gaps behind them | `n_sig - 1`, or 0 |

### Worked

```
contig_len          = 89,547
div_n_sites        =    164   mean_gap =     23.69   n_gaps =    163
hf_n_sites         =    169   mean_gap =     10.34   n_gaps =    168
lf_n_sites         =    196   mean_gap =    318.32   n_gaps =    195
union_n_sites      =    264   mean_gap =    236.02   n_gaps =    263
```

Two of those are worth checking by hand.

**`div_n_gaps`** is 163, exactly one fewer than the
164 sites — 164 points have 163 gaps between them, always.

**`div_expected_gap`** = (89,547 + 1) / (164 + 1) = **542.72**, against an observed 23.69. The observed value is
0.044x the null, i.e. the tested sites on this contig are
about 23x more tightly packed than random placement predicts.

It is `L + 1` over `n + 1`, not `L / n`: sites never reach the contig ends, so
the span they cover is under `L`. `METHODS.md` derives this and checks it against
exhaustive enumeration.

**Zeros and NaNs are not the same.** Of this SGB's 69 contigs, 66 have
`div_n_sites = 0` — divergence never tested them — and their `div_mean_gap` is
`NaN`, meaning no gap could be measured, not that the gap is zero. A contig with
exactly one site is also `NaN` for the same reason, which is why every mean ships
beside its `n_gaps`.

## SGB-level columns

One row per SGB. Worked on `SLG191_DASTool_bins_89`, which has **102 reference contigs** and 69 contig rows.

| column | meaning | how it is calculated |
|---|---|---|
| `MAG_ID` | the SGB | |
| `species` | GTDB species label | from the scores table |
| `n_ref_contigs` | contigs in the reference | counted in the MAG-to-contig mapping. **Not** from the data, so it includes contigs no test reached |
| `ref_genome_len` | summed length of those | from the `.fai` |
| `{d}_n_contigs_with_sites` | contigs carrying >=1 tested site | `(contig_table.{d}_n_sites > 0).sum()` |
| `{d}_n_sites` | tested sites | `contig_table.{d}_n_sites.sum()` |
| `{d}_mean_gap` | pooled spacing | gap-count-weighted: `sum(mean_gap * n_gaps) / sum(n_gaps)` |
| `{d}_expected_gap` | pooled null | same weights, so the two stay comparable |
| `{d}_n_gaps` | total gaps | `contig_table.{d}_n_gaps.sum()` |
| `{d}_pct_genome_tested` | tested sites per bp of reference | `100 * n_sites / ref_genome_len` |
| `{d}_n_contigs_sig` | contigs carrying >=1 significant site | `contig_table.{d}_has_sig.sum()` |
| `{d}_pct_contigs_sig` | that, as a share of the tested contigs | `100 * n_contigs_sig / n_contigs_with_sites` |
| `{d}_n_sig_sites` | significant sites | `contig_table.{d}_n_sig.sum()` |
| `{d}_sig_mean_gap` | pooled spacing between significant sites | weighted by `sig_n_gaps` |
| `{d}_sig_n_gaps` | gaps behind that | |

### Worked: the pooled mean is not an average of averages

```
contigs contributing gaps      : 3
mean of the per-contig means   : 40.4800   <- WRONG, each contig counts equally
gap-count-weighted (as coded)  : 25.1094   <- what SLG191_DASTool_bins_89 reports
total gap distance / total gaps: 4,469 / 178 = 25.1094
```

The last two lines agree exactly, which is the point: pooling reconstructs the
answer you would get from the raw gaps, while averaging the per-contig means
lets a 4-gap contig weigh as much as a 163-gap one.

### Worked: the three contig counts nest

| | reference | with tested sites | with significant sites | % |
|---|---|---|---|---|
| `div` | 102 | 3 | 2 | 66.7% |
| `hf` | 102 | 4 | 2 | 50.0% |
| `lf` | 102 | 69 | 0 | 0.0% |
| `union` | 102 | 69 | — | — |

So `div_pct_contigs_sig` = 100 x 2 / 3 = **66.67%**. The
denominator is divergence's own tested contigs — not the reference count, and
not the contigs some other test reached.

## summary.tsv columns

One row per definition, pooled over every SGB. Same meanings as the SGB level
with the `{d}_` prefix dropped, since the definition is now a value in the
`definition` column rather than part of the column name. Two additions:

| column | meaning |
|---|---|
| `definition` | which tested set this row describes |
| `n_SGBs` | SGBs with at least one tested site under it |

`n_ref_contigs` is identical on every row — the reference does not depend on
which test you ask about.

### Worked: percentages are recomputed, not averaged

```
as coded : 100 * 1,155 / 3,008 = 38.40%
if you averaged the 62 per-SGB percentages instead: 34.31%
```

4.09 points
apart. The same weighting trap as the means, wearing different clothes.

### Zero and blank mean different things

```
definition  n_contigs_sig  pct_contigs_sig  n_sig_sites  sig_mean_gap
       div         1155.0            38.40      62490.0        161.74
        hf         1948.0            55.55      90904.0        234.55
        lf            0.0             0.00          0.0           NaN
     union            NaN              NaN          NaN           NaN
```

`lf` shows **0** — a measured absence. Control-group sites were tested
(353,873 of them, more than
divergence's 393,475) and none
reached `q < 0.05`. `union` shows **blank** — the question does not exist for it.
Filling union with 0 would falsely claim it was measured and came back empty.

## Things that change the numbers

1. **FDR is pooled genome-wide across all SGBs**, not within each one, so a
   per-SGB `q < 0.05` count is not a per-SGB FDR. Recomputing per SGB moves
   individual q values by as much as 0.83.
2. **HF = `fat` and LF = `control` is an assumption.** The Fig 1 code only ever
   says `Fat` and `Control`.
3. **`position` is 0-based** and position 0 genuinely occurs, so coordinates
   occupy `[0, L-1]` while there are `L` of them. Gap arithmetic is unaffected;
   anything compared against contig length is not.
4. **`min_p_value` is a minimum over four nucleotide p-values** with no
   correction for having taken a minimum, which inflates significant-site counts
   throughout. It affects counts, not spacing geometry.

