# Variable sites and their distribution across the Fig 1 SGBs

Per-SGB and **per-contig** counts of variable sites, how far apart they sit, and
how that compares with BH-significant divergence and parallelism sites — comparing
the END and PRE timepoints, over the 62 SGBs tested in the divergence test.

## What a variable site is

**A site that was TESTED** — a row in `p_value_summary`, significant or not.
The between-group and within-group tests preprocess separately and so test
different sites (21.5% of HF-tested sites are not divergence-tested), which means
the definition is not single-valued. Four are reported side by side:

| key | tested by | sites | contigs |
|---|---|---|---|
| `div` | `two_sample_paired_tTest` | 393,475 | 3,008 |
| `hf` | `single_sample_tTest` / `fat` | 332,435 | 3,507 |
| `lf` | `single_sample_tTest` / `control` | 353,873 | 3,949 |
| `union` | any of the three | 575,812 | 4,402 |

Each significance question is judged against its **own** tested set — Q5/Q6 vs
`div`, Q7/Q8 vs `hf`, Q9/Q10 vs `lf` — so a hit can never fall outside the
denominator its fraction is expressed against.

This replaces an earlier definition based on AlleleFlux's zero-difference filter.
That filter removes almost nothing (the resulting set covered 1–60% of each
genome), so it measured coverage breadth rather than variability.

## Layout

```
variable_site_distribution.py   analysis -> four output tables
plots.py                        figures -> 6 PNGs + a 62-page per-MAG PDF
METHODS.md                      how each metric is derived, worked on real data
DATA_DICTIONARY.md              every file and column, with worked sums
```

Both markdown files are maintained by hand; keep them in step with
`variable_site_distribution.py` when the columns or the arithmetic change.

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
They are the two SGB-level spacing views with contigs as the unit, so a page answers
"how is spacing distributed inside this genome" rather than repeating the totals.
Output goes to `--outdir`,
which lives on scratch beside the tables (`$OUT/figs`) rather than in the repo — the
per-MAG PDF alone is ~360 KB and regenerates in 25 s, so it is output, not source.
It reads only the tables and computes no new statistics, so a figure can never disagree
with them. Colours are the Okabe-Ito colourblind-safe set, validated rather than eyeballed;
`union` is deliberately neutral grey because it is an aggregate of the other three, not a
peer category.

## Outputs

| file | shape | one row per |
|---|---|---|
| `contigs/<MAG>.tsv` | 62 files | contig, for that SGB |
| `contig_level_all.tsv` | 4,402 x 38 | contig, all SGBs, with `MAG_ID` |
| `sgb_level.tsv` | 62 x 43 | SGB |
| `summary.tsv` | 4 x 13 | variable-site definition, all SGBs pooled |

Each level is derived from the one below, so they cannot disagree. Column prefixes
are the definition keys: `div_`, `hf_`, `lf_`, `union_`. Every mean spacing ships
with its `n_gaps` support count; every significant-contig count ships with its
`pct_contigs_sig` fraction.

### Column glossary

Throughout, **"sites" on its own means TESTED sites** — the variable-site definition.
Significant ones always carry `sig` in the name. `{d}` is one of `div`, `hf`, `lf`,
`union`; the `union` prefix has no `sig` columns because it is a union of tested sets
and has no test, hence no significance question, of its own.

**`contigs/<MAG>.tsv` and `contig_level_all.tsv`** — one row per contig, so there is
no "number of contigs" column here; the row count *is* that number.

| column | meaning |
|---|---|
| `contig_len` | contig length in bp, from the `.fai`; `-1` if absent from it |
| `{d}_n_sites` | tested sites on this contig |
| `{d}_mean_gap` / `{d}_median_gap` | spacing between tested sites (Q3) |
| `{d}_n_gaps` | gaps behind that mean; always `n_sites - 1`, or 0 |
| `{d}_expected_gap` | uniform-null spacing, `(L+1)/(n+1)` (Q4) |
| `{d}_n_sig` | significant sites on this contig |
| `{d}_has_sig` | whether it carries any; what the SGB level counts contigs by |
| `{d}_sig_mean_gap` / `{d}_sig_median_gap` | spacing between significant sites (Q6/Q8/Q10) |
| `{d}_sig_n_gaps` | gaps behind that mean |

**`sgb_level.tsv`** — one row per SGB.

| column | meaning |
|---|---|
| `n_ref_contigs` | **contigs in the reference** (Q2). From the mapping file, not the data, so it counts contigs no test ever reached |
| `ref_genome_len` | summed length of those contigs |
| `{d}_n_contigs_with_sites` | **contigs carrying at least one TESTED site** (Q1) |
| `{d}_n_sites` | tested sites, summed over contigs |
| `{d}_mean_gap` | pooled spacing, weighted by `n_gaps` — not a mean of contig means |
| `{d}_expected_gap` | pooled uniform null, weighted the same way |
| `{d}_pct_genome_tested` | `100 * n_sites / ref_genome_len` |
| `{d}_n_contigs_sig` | **contigs carrying at least one SIGNIFICANT site** (Q5/Q7/Q9) |
| `{d}_pct_contigs_sig` | that count over `{d}_n_contigs_with_sites` — this test's own tested contigs, which is the denominator the questions imply |
| `{d}_n_sig_sites` | significant sites, summed over contigs |
| `{d}_sig_mean_gap` | pooled spacing between significant sites |

So for one SGB the three contig counts nest: `n_ref_contigs` >= `{d}_n_contigs_with_sites`
>= `{d}_n_contigs_sig`. Worked example, `SLG191_DASTool_bins_89`:

| | reference | with tested sites | with significant sites |
|---|---|---|---|
| `div` | 102 | 3 | 2 |
| `hf` | 102 | 4 | 2 |
| `lf` | 102 | 69 | 0 |
| `union` | 102 | 69 | — |

**`summary.tsv`** — one row per definition, pooled over all 62 SGBs. Same column
meanings without the `{d}_` prefix, plus `definition` and `n_SGBs` (SGBs with at
least one tested site under that definition).

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
runs 0.001–0.999 (median 0.314) for `div`, i.e. tested sites
are clustered roughly threefold. That comparison only became informative once a
variable site was defined as a tested site; under the old dense definition the
observed mean and the null coincide arithmetically.

**LF parallelism is empty.** No control-group site reaches q < 0.05 anywhere; the
minimum q genome-wide is 0.195. Q9 is 0% for every SGB and Q10 undefined. Note LF
is not short of tested sites (353,873, more than divergence's
393,475) — this is an absence of signal, not of data.

## Known caveats

1. **FDR is pooled across all SGBs**, not within each one, so a per-SGB q<0.05
   count is not a per-SGB FDR. Recomputing per SGB moves q by up to 0.83.
2. **HF = `fat`, LF = `control`** is an assumption. The Fig 1 code only says
   `Fat`/`Control`; the sole HF/LF mapping is in
   `Revision/AlleleFlux/notebooks/per_litter_boxplots.ipynb`.
3. **The uniform null is `(L+1)/(n+1)`, not `L/n`.** Sites never reach the contig
   ends, so the span they cover is under L. `METHODS.md` derives this and checks
   it against exhaustive enumeration. Do not simplify it away.
4. **`p_value_summary` needs `test_type` AND `group_analyzed` pinned.** The code
   asserts uniqueness rather than deduping; a dedupe would silently merge the two
   diet groups. The assertion was tested by forcing the collision (203,926 dupes).
5. **`position` is 0-based** and position 0 occurs; `min_p_value` is a minimum over
   four nucleotides with no selection correction, which inflates counts.

