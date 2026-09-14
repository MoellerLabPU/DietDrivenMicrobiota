# `alleleflux-baseline-presence`

**Question:** for every site the statistics called significant, was the winning
allele already in the mice at the baseline timepoint (standing variation), or did
it appear afterwards (de novo candidate)?

One run = one comparison directory (`pre_end-fat_control`) and one test. It
reads only AlleleFlux's own outputs; the strain-turnover step is optional.

**Naming.** Every column name and label that refers to a timepoint uses the
comparison's own timepoint names, exactly as the metadata spells them: `pre`
and `end` for `pre_end-fat_control` (so `n_pre_samples_covered`,
`allele_absent_at_end`), `5mo` and `22mo` for a DRiDO comparison. Nothing is
hard-coded; "earlier" and "later" below mean the first and second timepoint of
the comparison.

```
alleleflux-baseline-presence \
  --run_dir  <run>/longitudinal  --comparison pre_end-fat_control \
  --summary  two_sample_paired   --test_type two_sample_paired_tTest \
  --profiles_dir <profiles>  --metadata <metadata.tsv> \
  --fasta <ref.fa>  --mag_mapping <mapping.tsv>  --output_dir <out> \
  [--threshold_column q_value] [--threshold 0.05] [--min_cov 5] [--min_freq 0.05] \
  [--turnover_dir <strain_turnover outputs>] [--mags MAG ...] [--cpus N]
```

## Flow

1. **Roster.** From the metadata (the run's original input file, not a per-MAG
   `inputMetadata` file, whose rosters differ per MAG), keep the two groups and
   two timepoints of the comparison.
2. **Sites.** Open `p_value_summary_{summary}_*.tsv`, keep rows whose
   `test_type` equals `--test_type` and whose `--threshold_column` is
   `<= --threshold`.
3. **Alleles.** For each site, open the MAG's per-base test file and keep the
   base(s) whose p-value equals the site's `min_p_value`. At a biallelic site
   the two alleles are mirror images and both tie; both are reported, with
   `n_alleles_tied_at_min_p = 2`.
4. **Reads.** One job per (MAG, sample), for **every** sample of the comparison,
   both timepoints alike: load the profile once, and for every candidate (site,
   allele) count that allele's reads and give a status.
5. **Tables.** Join site and sample context, label each `end` row against its
   own mouse's `pre`, write the long table and the per-site summary.

## Filters, in the order they apply

| filter | flag / rule | what it removes |
|---|---|---|
| test family | `--summary` (two_sample_paired, two_sample_unpaired, single_sample, lmm, lmm_across_time) | other families' files. CMH is not offered: it has no per-base p-value, so it cannot name the allele |
| test statistic | `--test_type`, exactly as the summary spells it | other statistics in the same file (Wilcoxon, `_abs`, …) |
| significance | `--threshold_column` `<=` `--threshold` | non-significant sites |
| allele | per-base p `==` site `min_p_value` | bases that were not the significant one |
| depth gate | `--min_cov` (default 5, 1 = off) total reads at the position in that sample | samples too thin to say anything: status `not_covered`, out of every denominator |
| presence | reads `>=` null-model bar for that depth **and** `>= --min_freq` | reads that could be sequencing error |

**Per-sample status** (`allele_status`), after the depth gate:

| status | meaning |
|---|---|
| `present` | clears the bar and the frequency floor |
| `below_detection` | 1+ reads but under the bar or under 5 % |
| `absent` | zero reads of the allele |
| `not_covered` | fewer than `min_cov` reads at the position, **or** the MAG has no profile for that sample at all (no reads mapped; the run warns per MAG) |

"Covered" is about the **position** (enough reads to look); "present" / "absent"
are about the **allele** at a covered position. "Covered and absent" therefore
means "we looked properly and it was not there".

## Output 1: long table, `{comparison}_{summary}_{stat}_baseline_presence.tsv.gz`

One row per site × allele × sample. Real rows, MAG bin.012, allele **G** at
position 13509 (paired tTest, q = 0.044):

| sample | mouse | group | time | allele_reads | total_reads | bar | status | origin_in_own_mouse |
|---|---|---|---|---|---|---|---|---|
| SLG194 | 533 | fat | pre | 14 | 19 | 3 | present | |
| SLG1120 | 546 | control | end | 3 | 5 | 3 | present | standing_variation |
| SLG1121 | 547 | control | end | 7 | 9 | 3 | present | pre_not_covered |
| SLG1102 | 542 | control | end | 0 | 2 | 2 | not_covered | end_not_covered |

How the per-sample columns are computed (filters in **bold**):

| column | computed as | filter |
|---|---|---|
| allele_reads | reads of this base at this position in this sample | none |
| total_reads | A+C+G+T at the position | none |
| detection_threshold_reads | null-model bar at `total_reads` (`--fdr`, `--min_base_quality`) | none, it is the bar itself |
| allele_frequency | allele_reads / total_reads; NaN when not covered | **depth gate** |
| allele_status | `not_covered` if total_reads < min_cov; else `present` if allele_reads >= bar **and** frequency >= min_freq; else `absent` if 0 reads; else `below_detection` | **depth gate, bar, 5 % floor** |
| allele_present | allele_status == present | same |
| origin_in_own_mouse | this row's status + the same mouse's `pre` status, table below | same, on both timepoints |

`origin_in_own_mouse` is filled on `end` rows only and reads **both** timepoints
of the same mouse:

| end status | own pre status | label |
|---|---|---|
| present | present | `standing_variation` |
| present | below_detection | `de_novo_candidate_below_detection_at_pre` |
| present | absent | `de_novo_candidate` |
| present | not_covered | `pre_not_covered` |
| present | no pre sample | `no_pre_sample` |
| below_detection | any | `allele_below_detection_at_end` |
| absent | any | `allele_absent_at_end` |
| not_covered | any | `end_not_covered` |

Other columns: mag_id, contig, position, gene_id, test_type, group_analyzed,
min_p_value, q_value, n_alleles_tied_at_min_p, replicate, time,
allele_frequency, allele_present, strain_background (only with
`--turnover_dir`), min_cov, min_freq.

## Output 2: summary, `..._baseline_presence_summary.tsv`

One row per site × allele, deliberately lean: the numbers the baseline question
asks for, nothing else. Everything further is in the long table. The columns
come in two kinds and it matters which is which:

* **Sample counts** (`n_*_samples_*`, `n_replicates_*`, `n_mice_*`,
  `origin_*`) apply the presence rule: a sample counts only if it passed the
  depth gate, and "present" means the allele cleared the bar and the 5 % floor.
* **Read counts** (`total_reads_*`, `allele_reads_*`, `allele_frequency_*`)
  apply **no filter at all**: every read at the position from every sample at
  that timepoint, thin samples and profile-less samples (0 reads) included.
  That is what makes "the allele was never seen in 1,000 reads at pre" a
  frequency bound of < 1/1,000 rather than a statement about the filtered
  subset.

The same G at 13509:

| column | value | computed as |
|---|---|---|
| mag_id, contig, position, gene_id, group_analyzed, allele | … | site and allele identity (`group_analyzed` is blank for two-sample tests; it separates the per-group rows of single-sample tests) |
| n_alleles_tied_at_min_p, q_value | 2, 0.044 | copied from the site |
| origin_any_mouse | standing_variation | allele_not_present_at_end if no covered end sample shows the allele; else standing_variation if any covered pre sample has it present; else de_novo_candidate_below_detection_at_pre if any pre has it below detection; else de_novo_candidate if any covered pre is absent; else pre_not_covered |
| n_pre_samples_allele_present / n_pre_samples_covered | 12 / 14 | present pre samples over pre samples deep enough to judge (filtered) |
| n_replicates_with_allele_at_pre | 6 | distinct replicates among the present pre samples (filtered) |
| pre_mice_allele_present | 530,532,533,… | their subjectIDs |
| n_mice_standing_variation | 6 | end samples whose own pre had the allele present |
| n_mice_de_novo_candidate | 1 | end present, own pre covered and absent |
| n_mice_de_novo_candidate_below_detection_at_pre | 1 | end present, own pre had reads under the bar |
| total_reads_pre | 168 | A+C+G+T at the position, summed over **all 30 pre samples** (unfiltered; 16 of them were too thin to count above) |
| allele_reads_pre | 125 | G reads among them (unfiltered) |
| allele_frequency_pre | 0.744 | 125 / 168; NaN when total is 0 |
| total_reads_end | 245 | same, over all 30 end samples |
| allele_reads_end | 179 | |
| allele_frequency_end | 0.731 | 179 / 245 |

Read as a sentence: G was already present in 12 of the 14 baseline mice we
could see, in 6 of 8 replicates; of the mice we can check against their own
baseline, 6 had it, 1 did not, 1 had a trace. Pooling every read, 125 of 168
at pre carried G. Standing variation.

For a de novo candidate the read columns give the requested upper bound: if
`allele_reads_pre` is 0 and `total_reads_pre` is 1,000, the allele was absent
at pre or below 1/1,000. Note the two kinds can disagree on purpose: a site can
be `de_novo_candidate` (no covered pre sample had it present) while
`allele_reads_pre` is 3, because those reads sat in samples under the bar.

The three framings, any mouse / same replicate / same mouse, are all here. Which
one is the headline is a design question: littermates sharing a colony justify
"any mouse"; outbred DRiDO mice, half without an earlier sample, justify "same
mouse". The per-mouse verdicts for every `end` sample, including the ones with
no origin to assign (pre not covered, no pre sample, allele not seen at end),
are in the long table's `origin_in_own_mouse` column.

## Answering the original question

> How many BH-corrected significant divergence alleles at END were detected in
> any PRE mouse? For each position, which replicates did or did not have the
> allele at PRE, and how many reads per mouse?

| part of the question | where |
|---|---|
| BH-corrected significant divergence alleles at END | run with `--summary two_sample_paired --test_type two_sample_paired_tTest --threshold_column q_value --threshold 0.05`; one summary row per site × allele. Count sites with `n_alleles_tied_at_min_p` in mind: a biallelic site contributes two rows |
| detected in **any** PRE mouse | `origin_any_mouse == standing_variation` (equivalently `n_pre_samples_allele_present >= 1`); the headline number is the share of rows |
| detected in the **same** mouse | `n_mice_standing_variation` vs `n_mice_de_novo_candidate` (+ `_below_detection_at_pre`), per site; the stricter framing |
| which replicates had it at PRE | `n_replicates_with_allele_at_pre` and `pre_mice_allele_present` (mouse → replicate via the metadata); the mice that lacked it or had only a trace are in the long table's pre rows by `allele_status` |
| reads per mouse | long table, pre rows: `allele_reads`, `total_reads`, `detection_threshold_reads`, `allele_status` |
| total reads per position across all mice, and reads carrying the allele, at PRE and END (follow-up ask) | summary: `total_reads_pre`, `allele_reads_pre`, `allele_frequency_pre` and the `_end` trio; unfiltered sums. Upper bound on a de novo allele's pre frequency = 1 / `total_reads_pre` when `allele_reads_pre` is 0 |
| MAG, contig | on every row of both tables (`mag_id`, `contig`, `position`, `gene_id`) |
| taxonomy | not written by the command (dropped by design); join `gtdbtk.bac120.summary.tsv` on `user_genome == mag_id` |

Standing variation vs mutation, per site, is `origin_any_mouse` under the
loosest framing and the `n_mice_*` split under the strictest; both are reported
so the choice of framing is made in the analysis, not baked into the file.

## Undoing a filter without rerunning

The long table keeps the raw numbers behind every verdict:

| to relax | look at |
|---|---|
| the depth gate | `total_reads` (kept on `not_covered` rows too) |
| the null-model bar | `allele_reads` against `detection_threshold_reads`; "any read" is `allele_reads > 0` |
| the 5 % floor | `allele_frequency` |
| the significance cutoff | rerun with `--threshold 1` or `--threshold_column min_p_value`; unread sites have no rows |

The summary is the roll-up under the filters as run and cannot be undone this way.

## Scale

Profiles load once per (MAG, sample). The pre_end divergence run, 57 MAGs × 60
samples, took 23 min and 6.4 GB on 8 cores (SLURM job 13697897); run it on
SLURM, not the login node. Zero significant sites (e.g. Wilcoxon at n = 8)
writes header-only files and exits 0.
