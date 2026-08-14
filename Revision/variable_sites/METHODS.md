# How each reported number is calculated

Every figure below was computed directly from the run, so the arithmetic here is
the arithmetic the analysis performs. The generator that produced it has been
removed, so this file is now maintained by hand — update it alongside any change
to `variable_site_distribution.py`.

Worked example throughout: SGB **`SLG191_DASTool_bins_89`**.

---

## What is produced

```
  AlleleFlux run on /scratch
      |  p_value_summary/pre_end-fat_control/
      |      ..._two_sample_paired_pre_end.tsv   <- divergence, tested + q
      |      ..._single_sample_pre_end.tsv       <- HF and LF, tested + q
      |  redo_representative_megamag.fa.fai + megamag_mapping.tsv  <- geometry
      v
  variable_site_distribution.py
      v
  contigs/<MAG>.tsv      one file per SGB, one row per contig
  contig_level_all.tsv   all of the above concatenated, with a MAG_ID column
  sgb_level.tsv          one row per SGB
  summary.tsv            one row per variable-site definition, all SGBs pooled
```

Each level is derived from the one below it, so they cannot disagree: an SGB row
is built from that SGB's contig table, and the summary from the SGB table.

The analysis reads **only** `p_value_summary` and the reference geometry. It no
longer touches `allele_analysis/` at all, which is why it runs in seconds.

## Step 1 — what a 'variable site' is

**A variable site is a site that was TESTED.** If a statistical test evaluated
that position, it counts. Concretely, it is a row in `p_value_summary` for the
relevant test — irrespective of whether the result was significant.

That definition is not single-valued, because AlleleFlux preprocesses the
between-group and within-group tests separately and they end up testing
different sites. So four definitions are reported side by side:

| key | tested by | sites (62 SGBs) | contigs |
|---|---|---|---|
| `div` | `two_sample_paired_tTest` | 393,475 | 3,008 |
| `hf` | `single_sample_tTest` / `fat` | 332,435 | 3,507 |
| `lf` | `single_sample_tTest` / `control` | 353,873 | 3,949 |
| `union` | any of the three | 575,812 | 4,402 |

They genuinely differ: **71,441 HF-tested sites are not divergence-tested**
(21.5% of the HF set). That is why each significance question
is judged against its **own** tested set — Q5/Q6 against `div`, Q7/Q8 against
`hf`, Q9/Q10 against `lf`. A significant site can then never fall outside the
denominator its fraction is expressed against.

`union` is the fourth definition: tested by anything. A site tested by two
families appears in both source frames, so the union is deduplicated on
`(contig, position)`, keeping the smallest q_value.

### Why not the zero-difference filter?

An earlier version defined a variable site as one surviving AlleleFlux's
zero-difference filter, read from `allele_analysis`. That filter drops a position
only when the summed *signed* allele-frequency change, pooled over every mouse in
both groups, is exactly zero — so it removes almost nothing, and the resulting
set covered 1–60% of each genome. It measured coverage breadth, not variability.
The tested-site definition is roughly 240x smaller and is the set the
significance questions are actually drawn from.

## Q1 — how many contigs contain at least one variable site?

Count distinct `contig` values in each tested set.

```python
n_contigs_with_sites = tested.contig.nunique()
```

For `SLG191_DASTool_bins_89`, out of 102 reference contigs (Q2 below):

| definition | contigs with >=1 tested site | tested sites |
|---|---|---|
| `div` | 3 | 181 |
| `hf` | 4 | 194 |
| `lf` | 69 | 2,053 |
| `union` | 69 | 2,136 |

This can never exceed Q2. If it ever did, the mapping and the run would be
mismatched, which is worth treating as a hard error rather than a curiosity.

## Q2 — how many contigs is the SGB reference?

Pure bookkeeping, independent of any data: count the SGB's contigs in the
MAG-to-contig mapping, and sum their lengths from the FASTA index.

```
megamag_mapping.tsv   mag_id -> contig_id      (one row per contig)
*.fa.fai              contig -> length         (column 2)
```

For `SLG191_DASTool_bins_89`: **n_ref_contigs = 102**, ref_genome_len = 1,867,134 bp.

The join is exact — every mapped contig has a length — and the summed lengths
reproduce the `genome_size` column AlleleFlux writes in its own QC tables for
all 160 MAGs, which is the cross-check that mapping and reference correspond.

## Q3 — average distance between variable sites

Distances are computed **inside a contig only**. A gap spanning two contigs has
no meaning: they are separate DNA molecules of unknown relative position and
orientation. Sort the positions within a contig, take successive differences,
then pool those differences across the SGB's contigs.

```python
gaps = np.diff(np.sort(positions_on_this_contig))    # n sites -> n-1 gaps
mean_gap = pooled over contigs, weighted by gap count
```

Worked on `k141_116087`, the divergence-tested contig of `SLG191_DASTool_bins_89`
carrying the most sites (164). First ten sorted positions:

```
positions : 30449, 32603, 32604, 32657, 32683, 32716, 32721, 32724, 32726, 32770
gaps      :   2154,  1,  53,  26,  33,  5,  3,  2,  44
```

That contig contributes **163 gaps** (always one fewer than its site
count), mean 23.7 nt, median 5 nt.

Pooled across `SLG191_DASTool_bins_89`:

| definition | mean gap (nt) | n_gaps |
|---|---|---|
| `div` | 25.1 | 178 |
| `hf` | 12.3 | 190 |
| `lf` | 369.6 | 1,984 |
| `union` | 354.7 | 2,067 |

Two things that silently change the answer, so both are explicit:

**Contigs with fewer than two sites contribute nothing.** One site yields zero
gaps. Such contigs are neither counted as a gap of zero nor as a gap of the
contig length; they are simply absent from the pool. That is why every mean is
published next to its `n_gaps` — a mean over 3 gaps and a mean over 40,000 read
identically otherwise.

**Pooling is not averaging per-contig means.** Pooling weights each contig by
how many gaps it contributes, so a dense contig dominates a sparse one;
averaging contig means would weight them equally. The SGB figure is rebuilt from
the per-contig means and gap counts, so the two levels cannot disagree.

## Q4 — average distance between randomly drawn sites

The yardstick for Q3: if the same number of sites were scattered blindly along
the same DNA, how far apart would neighbours be?

### The sidewalk

Picture a sidewalk of 10 squares. Toss 3 pebbles onto 3 different squares, then
walk from your first pebble to your last.

- **How far do you walk?** From the first pebble to the last. You do *not* walk
  the whole sidewalk — the squares before your first pebble and after your last
  are never crossed. This is the step intuition skips.
- **How many hops?** Three pebbles, two hops between them: *n* - 1.
- **Average hop** = distance walked / 2.

Averaged over every possible throw, 3 pebbles on 10 squares span only **5.5**
squares — a little over half the sidewalk — so the average hop is 5.5/2 = 2.75.
The tempting `L/n` = 10/3 = 3.33 is wrong because it pretends you walked all ten.

### The formula

```
E[gap] = (L + 1) / (n + 1)
```

`L + 1`, not `L` and not `L - 1`. Nothing is added to the DNA — the contig is L
bases. The `+1`s are what a ratio of two binomial coefficients simplifies to:

1. The chance the first site sits at position *k* or later is the chance all *n*
   sites land in the `L-k+1` positions from *k* on: `C(L-k+1, n) / C(L, n)`.
2. For a positive whole number, E[X] = P(X>=1) + P(X>=2) + ... (an outcome of 3
   is counted in exactly the first three terms). Summing step 1 over *k* and
   applying the hockey-stick identity gives `E[first] = C(L+1,n+1) / C(L,n)`,
   whose factorials cancel to **(L+1)/(n+1)**.
3. Mirroring the contig (*p* -> *L+1-p*) turns first into last, so
   `E[last] = (L+1) - E[first]`, hence `E[span] = (L+1)(n-1)/(n+1)`.
4. The mean gap is always span/(n-1), so dividing cancels the `(n-1)`.

Step 4 answers the obvious objection: Q3 divides by *n*-1 while Q4 appears to
divide by *n*+1. There is only ever one denominator, *n*-1. It cancels.

### The proof: every arrangement enumerated

For a small contig you can simply list every possibility. All `C(L,n)` ways to
place *n* sites on a contig of L = 10, with span and mean gap averaged over all:

| n | avg span | hops (n-1) | span/(n-1) | (L+1)/(n+1) | L/n |
|---|---|---|---|---|---|
| 2 | 3.6667 | 1 | **3.6667** | **3.6667** | 5.000 |
| 3 | 5.5000 | 2 | **2.7500** | **2.7500** | 3.333 |
| 4 | 6.6000 | 3 | **2.2000** | **2.2000** | 2.500 |
| 5 | 7.3333 | 4 | **1.8333** | **1.8333** | 2.000 |
| 8 | 8.5556 | 7 | **1.2222** | **1.2222** | 1.250 |

The two bold columns agree exactly in every row — the same quantity written two
ways, not two competing methods. `L/n` matches neither.

Hockey-stick check for step 2 (L=10, n=3): sum of C(m,3) for m=3..10 = 330, and C(11,4) = 330.

### How the contigs are combined

```python
num += (n - 1) * (contig_len + 1) / (n + 1)   # gap count x that contig's expectation
den += (n - 1)                                 # accumulate gap counts
expected = num / den
```

The weighting is forced, not stylistic: pooling observed gaps already lets each
contig count in proportion to its gap count, so the null must match. Two 1,000 bp
contigs with 101 and 3 sites have expected gaps of 9.9 and 250.3; averaging those
plainly gives 130, weighting by gap count gives 14.5, and the observed pooled
mean is 14.6.

### What it says

For `SLG191_DASTool_bins_89`:

| definition | observed mean gap | expected under uniform | ratio |
|---|---|---|---|
| `div` | 25.1 | 931.7 | 0.027 |
| `hf` | 12.3 | 898.0 | 0.014 |
| `lf` | 369.6 | 583.5 | 0.633 |
| `union` | 354.7 | 562.0 | 0.631 |

Across all 62 SGBs the observed/expected ratio runs
0.001–0.999 (median **0.314**) for `div` and
0.150–1.006 (median **0.580**) for `union`.

A ratio below 1 means sites sit closer together than blind scattering would put
them — tested sites are clustered, by roughly threefold at the median. This is a
real signal, and it is worth noting that it only became visible once a variable
site was defined as a *tested* site. Under the old zero-difference definition the
sites were so dense that they spanned essentially the whole contig, which forces
the observed mean and the null to coincide arithmetically.

## Q5 and Q6 — divergence between the diet groups

**Which test.** Divergence is the between-group question: do fat and control
move apart? We use

    test_type == 'two_sample_paired_tTest'   and   q_value < 0.05

The paired t-test is the project's canonical divergence test, and the only one
whose label-permuted runs return no significant MAGs at all — which is what makes
a real signal interpretable.

**The file is not one row per site, and needs no deduplication.** It carries one
row per `(site x test_type)`, four variants in all, each with 393,475 rows:
`tTest`, `tTest_abs`, `Wilcoxon`, `Wilcoxon_abs`. Pinning `test_type` already
leaves exactly one row per site. The parallelism file is the same with one extra
key, `group_analyzed`.

A `drop_duplicates()` here would be worse than redundant: forget the
`group_analyzed` filter and the two diet groups collide on the same position
(721,676 rows for 517,704 sites over the unthresholded rows), and deduping would
quietly merge fat and control into one plausible-looking, wrong set. The code
**asserts** uniqueness instead and aborts. Forcing that collision confirms the
guard fires: it caught 203,926 duplicate rows.

**Q5 counts contigs against the divergence-tested denominator.** Of the contigs
carrying a divergence-*tested* site, how many carry a *significant* one? Because
both sets come from the same test, a significant site is always inside the tested
set and nothing can be dropped:

```python
denom = contigs where div_n_sites > 0
n_contigs_div_sig = contigs where div_n_sig > 0
pct_contigs_div_sig = 100 * n_contigs_div_sig / denom
```

For `SLG191_DASTool_bins_89`: **2 of 3 divergence-tested contigs = 66.67%**, over 70 significant sites.

Across the 62 SGBs that fraction runs 0.0%–100.0% (median 35.6%).

**Q6** applies exactly the Q3 procedure to the divergence-significant positions.
Significant sites are sparser than tested ones, so many contigs hold only one and
contribute no gap; `div_sig_n_gaps` is how you tell a well-supported mean from a
fragile one.

## Q7 to Q10 — parallelism within each diet group

Parallelism is the within-group question: do allele frequencies move
*consistently* inside one diet group over time, irrespective of the other group?
Both come from the same file, separated by one column:

```
p_value_summary_single_sample_pre_end.tsv
  test_type      == 'single_sample_tTest'
  group_analyzed == 'fat'      -> HF   (Q7, Q8)
  group_analyzed == 'control'  -> LF   (Q9, Q10)
```

The contig counting and gap arithmetic are identical to Q5/Q6, each against its
own tested denominator.

For `SLG191_DASTool_bins_89`: **HF 2 of 4 HF-tested contigs = 50.00%**, over 90 sites.

Across the 62 SGBs: 0.0%–96.0% (median 41.4%).

### The LF result is real, not a bug

**LF: 0 significant sites in the entire dataset.** The smallest q value
genome-wide is 0.195 for the t-test (0.509 for Wilcoxon), so
nothing reaches q < 0.05. Q9 is 0% for every SGB and Q10 is undefined for every
SGB; the code emits 0 and NaN respectively.

Read biologically that is a finding: the high-fat group shows parallel
allele-frequency change and the control group does not. Anyone re-running this
and getting non-zero LF numbers has changed the test, the threshold, or the FDR
scope. Note LF is not short of *tested* sites — it has 353,873 of them,
more than divergence's 393,475 — so this is an absence of signal, not
an absence of data.

**HF and LF map to `fat` and `control` by assumption.** The Fig 1 code only says
`Fat` and `Control`; the only HF/LF mapping in the repo is in
`Revision/AlleleFlux/notebooks/per_litter_boxplots.ipynb`. If the control diet is
not a low-fat diet, Q9 and Q10 were asking about a group this run does not hold.

## Two properties to state when presenting the numbers

### FDR is pooled across all SGBs, not within each one

Benjamini-Hochberg correction is applied once over every site of a given test
across all MAGs together. So 'BH-significant for this SGB' means 'significant in
a genome-wide ranking', and a strongly-affected SGB shifts everyone else's q
values. Recomputing FDR per SGB moves individual q values by as much as 0.83.

### Smaller cautions

`position` is **0-based** and position 0 genuinely occurs, so coordinates occupy
[0, L-1] while there are L of them; Q4's formula is written in terms of the
*count* L, so it already accounts for this, and differences between positions are
unaffected either way. Separately, `min_p_value` is the minimum of four nucleotide
p-values with no correction for having taken a minimum, which inflates
significant-site counts throughout — it affects counts, not spacing geometry.

## Cross-check against the published tables

Recomputed here, independently of the analysis script, for `SLG191_DASTool_bins_89`:

| column | recomputed | sgb_level.tsv |
|---|---|---|
| n_ref_contigs | 102 | 102  |
| ref_genome_len | 1,867,134 | 1,867,134  |
| div_n_sites | 181 | 181  |
| div_n_contigs_with_sites | 3 | 3  |
| hf_n_sites | 194 | 194  |
| hf_n_contigs_with_sites | 4 | 4  |
| lf_n_sites | 2,053 | 2,053  |
| lf_n_contigs_with_sites | 69 | 69  |
| union_n_sites | 2,136 | 2,136  |
| union_n_contigs_with_sites | 69 | 69  |
| div_n_sig_sites | 70 | 70  |
| hf_n_sig_sites | 90 | 90  |
| lf_n_sig_sites | 0 | 0  |

All 13 values match. If they ever stop
matching, one of the two code paths has changed and the discrepancy is the bug.

