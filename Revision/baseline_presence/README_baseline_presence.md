# `alleleflux-baseline-presence`

**Question:** for every site the statistics called significant, was the winning
allele already in the mice at the baseline timepoint (standing variation), or did
it appear afterwards (de novo candidate)?

One run = one comparison directory (`pre_end-fat_control`) and one test. It
reads only AlleleFlux's own outputs; the strain-turnover step is optional.

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
   two timepoints of the comparison; label each sample `t0` (earlier) or `t1`
   (later).
2. **Sites.** Open `p_value_summary_{summary}_*.tsv`, keep rows whose
   `test_type` equals `--test_type` and whose `--threshold_column` is
   `<= --threshold`.
3. **Alleles.** For each site, open the MAG's per-base test file and keep the
   base(s) whose p-value equals the site's `min_p_value`. At a biallelic site
   the two alleles are mirror images and both tie; both are reported, with
   `n_alleles_tied_at_min_p = 2`.
4. **Reads.** One job per (MAG, sample), for **every** sample of the comparison,
   t0 and t1 alike: load the profile once, and for every candidate (site,
   allele) count that allele's reads and give a status.
5. **Tables.** Join site and sample context, label each t1 row against its own
   mouse's t0, write the long table and the per-site summary.

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

| sample | mouse | group | role | allele_reads | total_reads | bar | status | origin_in_own_mouse |
|---|---|---|---|---|---|---|---|---|
| SLG194 | 533 | fat | t0 | 14 | 19 | 3 | present | |
| SLG1120 | 546 | control | t1 | 3 | 5 | 3 | present | standing_variation |
| SLG1121 | 547 | control | t1 | 7 | 9 | 3 | present | t0_not_covered |
| SLG1102 | 542 | control | t1 | 0 | 2 | 2 | not_covered | t1_not_covered |

How the per-sample columns are computed (filters in **bold**):

| column | computed as | filter |
|---|---|---|
| allele_reads | reads of this base at this position in this sample | none |
| total_reads | A+C+G+T at the position | none |
| detection_threshold_reads | null-model bar at `total_reads` (`--fdr`, `--min_base_quality`) | none, it is the bar itself |
| allele_frequency | allele_reads / total_reads; NaN when not covered | **depth gate** |
| allele_status | `not_covered` if total_reads < min_cov; else `present` if allele_reads >= bar **and** frequency >= min_freq; else `absent` if 0 reads; else `below_detection` | **depth gate, bar, 5 % floor** |
| allele_present | allele_status == present | same |
| origin_in_own_mouse | this row's status + the same mouse's t0 status, table below | same, on both timepoints |

`origin_in_own_mouse` is filled on t1 rows only and reads **both** timepoints of
the same mouse:

| t1 status | own t0 status | label |
|---|---|---|
| present | present | `standing_variation` |
| present | below_detection | `de_novo_candidate_below_detection_at_t0` |
| present | absent | `de_novo_candidate` |
| present | not_covered | `t0_not_covered` |
| present | no t0 sample | `no_t0_sample` |
| below_detection | any | `allele_below_detection_at_t1` |
| absent | any | `allele_absent_at_t1` |
| not_covered | any | `t1_not_covered` |

Other columns: mag_id, contig, position, gene_id, test_type, group_analyzed,
min_p_value, q_value, n_alleles_tied_at_min_p, replicate, time,
allele_frequency, allele_present, strain_background (only with
`--turnover_dir`), min_cov, min_freq.

## Output 2: summary, `..._baseline_presence_summary.tsv`

One row per site × allele, deliberately lean: the numbers the baseline question
asks for, nothing else. Everything further is in the long table. "covered" =
passed the depth gate; "present" = allele_status is `present`. The same G at 13509:

| column | value | computed as |
|---|---|---|
| mag_id, contig, position, gene_id, group_analyzed, allele | … | site and allele identity (`group_analyzed` is blank for two-sample tests; it separates the per-group rows of single-sample tests) |
| n_alleles_tied_at_min_p, q_value | 2, 0.044 | copied from the site |
| origin_any_mouse | standing_variation | allele_not_present_at_t1 if no covered t1 sample shows the allele; else standing_variation if any covered t0 sample has it present; else de_novo_candidate_below_detection_at_t0 if any t0 has it below detection; else de_novo_candidate if any covered t0 is absent; else t0_not_covered |
| n_t0_samples_allele_present / n_t0_samples_covered | 12 / 14 | present t0 samples over t0 samples deep enough to judge |
| n_replicates_with_allele_at_t0 | 6 | distinct replicates among the present t0 samples |
| t0_mice_allele_present | 530,532,533,… | their subjectIDs |
| n_mice_standing_variation | 6 | t1 samples whose own t0 had the allele present |
| n_mice_de_novo_candidate | 1 | t1 present, own t0 covered and absent |
| n_mice_de_novo_candidate_below_detection_at_t0 | 1 | t1 present, own t0 had reads under the bar |

Read as a sentence: G was already present in 12 of the 14 baseline mice we
could see, in 6 of 8 replicates; of the mice we can check against their own
baseline, 6 had it, 1 did not, 1 had a trace. Standing variation.

The three framings, any mouse / same replicate / same mouse, are all here. Which
one is the headline is a design question: littermates sharing a colony justify
"any mouse"; outbred DRiDO mice, half without a t0 sample, justify "same mouse".
The per-mouse verdicts for every t1 sample, including the ones with no origin to
assign (t0 not covered, no t0 sample, allele not seen at t1), are in the long
table's `origin_in_own_mouse` column.

## Answering the original question

> How many BH-corrected significant divergence alleles at END were detected in
> any PRE mouse? For each position, which replicates did or did not have the
> allele at PRE, and how many reads per mouse?

| part of the question | where |
|---|---|
| BH-corrected significant divergence alleles at END | run with `--summary two_sample_paired --test_type two_sample_paired_tTest --threshold_column q_value --threshold 0.05`; one summary row per site × allele. Count sites with `n_alleles_tied_at_min_p` in mind: a biallelic site contributes two rows |
| detected in **any** PRE mouse | `origin_any_mouse == standing_variation` (equivalently `n_t0_samples_allele_present >= 1`); the headline number is the share of rows |
| detected in the **same** mouse | `n_mice_standing_variation` vs `n_mice_de_novo_candidate` (+ `_below_detection_at_t0`), per site; the stricter framing |
| which replicates had it at PRE | `n_replicates_with_allele_at_t0` and `t0_mice_allele_present` (mouse → replicate via the metadata); the mice that lacked it or had only a trace are in the long table's t0 rows by `allele_status` |
| reads per mouse | long table, t0 rows: `allele_reads`, `total_reads`, `detection_threshold_reads`, `allele_status` |
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

Profiles load once per (MAG, sample), ~3 s each. Sam's pre_end divergence:
57 MAGs × 60 samples ≈ 3,400 loads, ~3 h on 8 cores; run it on SLURM, not the
login node. Zero significant sites (e.g. Wilcoxon at n = 8) writes header-only
files and exits 0.
