# AlleleFlux MIDAS-filter runs (October 2026): `site_ratio` 3 and 10

Two reruns of the diet-manipulation data with AlleleFlux's per-site depth filter
(the MIDAS2 `merge_snps` rule, applied once per sample when the allele-frequency
cache is written). The two runs differ only in the filter's `site_ratio`:

| run | `site_ratio` | config | working directory (under `midas_filter_Oct_2026/`) |
|---|---|---|---|
| first run, submitted 2026-10-23 | 3.0 | `alleleflux_config_midas_filter.yaml` | `AlleleFlux_midas_filter` |
| follow-up, submitted 2026-10-27 | 10.0 | `alleleflux_config_midas_filter_ratio10.yaml` | `AlleleFlux_midas_filter_ratio10` |

The follow-up has its own section at the end. What follows applies to both runs.
Inputs (BAMs, reference, metadata sheet) and resources are
the same as the July 2026 revision run described in `README.md`; the analysis
is not. Per-MAG QC uses breadth ≥ 0.1 with no minimum mean depth (the values in
`figures/Fig1/alleleflux_config_mapq20.yml`, not the revision run's 0.5 and 1),
only the `pre_end` comparison is tested, and regional contrast, gene scores and
outlier detection are off. The full list of differences is under Files.

The `site_ratio` 3 run: set up 2026-10-23; submitted 2026-10-23 (driver job id in `slurm_scripts/logs/`; child jobs carry the WCKey below).

## Files

- **`alleleflux_config_midas_filter.yaml`** — the run configuration, derived from
  `alleleflux_config.yaml`. The diff against it is exactly:
  - `run_name` and `output.root_dir`;
  - `quality_control.breadth_threshold` 0.1 and `coverage_threshold` 0
    (0.5 and 1 in the revision run);
  - `use_regional_contrast`, `use_gene_scores` and `use_outlier_detection`
    switched to `False`;
  - the `pre_post` entry of `timepoints_combinations` commented out, so only
    `pre_end` runs;
  - the three strain-level switches written out as `False`
    (`use_pairwise_ani`, `use_strain_turnover`, `use_baseline_presence`; off in
    the revision run too, they just did not exist yet);
  - the `quality_control.site_depth_filter` block.
- **`slurm_scripts/run_alleleflux_midas_filter.sh`** — the driver job. Submit it
  with `sbatch` from `slurm_scripts/` (its logs go to `slurm_scripts/logs/`).
- **`alleleflux_config_midas_filter_ratio10.yaml`** and
  **`slurm_scripts/run_alleleflux_midas_filter_ratio10.sh`** — config and driver
  for the `site_ratio` 10 follow-up. The config differs from the ratio-3 config
  only in `run_name`, `output.root_dir`, `input.profiles_path` and `site_ratio`.

## Filter settings

Pipeline defaults, written out explicitly (decision 2026-10-23). The follow-up
run changes only `site_ratio`, from 3 to 10:

| key | value | meaning |
|---|---|---|
| `enabled` | `True` | filter on |
| `site_depth` | `0` | no minimum depth per position (MIDAS2 would use 5) |
| `site_ratio` | `3.0` (first run) / `10.0` (follow-up) | drop a position whose depth exceeds 3× (or 10×) the sample's mean coverage over covered positions |
| `min_allele_reads` | `0` | every read counts toward depth (MIDAS2 would use 2) |
| `keep_original_counts` | `True` | kept positions keep their original base counts (MIDAS2 would zero counts ≤ `min_allele_reads`) |

What the filter removed per sample is written next to each cache file
(`*_site_filter.json`) and rolled up by the `site_depth_filter_summary` rule
under `site_depth_filter/` in the output directory.

## Code: a frozen clone, not the development tree

The pipeline runs from a clone of AlleleFlux pinned to one commit, with its own
conda environment, so that development work in `~/AlleleFlux-dev` (Phase 2 of
the de novo rules, in progress at the time) can never reach a running job.

| what | where |
|---|---|
| repository | https://github.com/MoellerLabPU/AlleleFlux |
| clone | `/home/su2806/AlleleFlux-run-2026-10` |
| commit | `48e3db6` — "cache-first (Phase 2 Tasks 1-4): every per-site consumer reads the allele-frequency cache", 2026-10-19; full test suite (721) and both pipeline golden scenarios green at that commit |
| local branch in the clone | `frozen-run-2026-10` (points at that commit; nothing else on it) |
| conda env | `alleleflux-run-dm`, a clone of `alleleflux-dev` (same Snakemake, SLURM plugin and R versions), pinned to the clone by an editable install |
| SLURM profile | `<clone>/alleleflux/smk_workflow/slurm_profile_native/` |
| working directory | `/scratch/gpfs/AMOELLER/sidd/diet_manip/midas_filter_Oct_2026/AlleleFlux_midas_filter` |
| outputs | `<working directory>/AlleleFlux/longitudinal/` |

## Steps followed (reproducible)

```bash
module load anaconda3/2025.12

# 1. frozen clone at the chosen commit, on its own branch
git clone https://github.com/MoellerLabPU/AlleleFlux /home/su2806/AlleleFlux-run-2026-10
cd /home/su2806/AlleleFlux-run-2026-10
git checkout -b frozen-run-2026-10 48e3db6

# 2. dedicated env with the same dependency versions as the development env
conda create -n alleleflux-run-dm --clone alleleflux-dev

# 3. pin the env to the clone (run from the clone directory and nowhere else:
#    an editable install records the directory it was run from, last writer wins)
conda run -n alleleflux-run-dm pip install -e . --no-deps

# 4. verify the pin: must print /home/su2806/AlleleFlux-run-2026-10/alleleflux
grep -ohP "MAPPING.*?\{'alleleflux': '\K[^']+" \
  /home/su2806/.conda/envs/alleleflux-run-dm/lib/python*/site-packages/__editable___alleleflux_*_finder.py

# 5. dry run from the frozen env (done 2026-10-23: exit 0, 97 jobs before the
#    eligibility checkpoint = 90 profiling + 2 metadata + 2 QC + 2 eligibility + all;
#    the cache, statistics and scoring jobs appear only after the checkpoint)
conda run -n alleleflux-run-dm alleleflux run \
  -w /scratch/gpfs/AMOELLER/sidd/diet_manip/midas_filter_Oct_2026/AlleleFlux_midas_filter \
  -c /home/su2806/DietDrivenMicrobiota/additional_analyses/AlleleFlux/alleleflux_config_midas_filter.yaml \
  -p /home/su2806/AlleleFlux-run-2026-10/alleleflux/smk_workflow/slurm_profile_native/ --dry-run

# 6. submit
cd /home/su2806/DietDrivenMicrobiota/additional_analyses/AlleleFlux/slurm_scripts
mkdir -p logs && sbatch run_alleleflux_midas_filter.sh
```

The driver repeats step 4 at start-up and exits if the env points anywhere but
the clone. Every child job carries the WCKey `alleleflux_AlleleFlux_midas_filter`,
so `sacct --wckeys=alleleflux_AlleleFlux_midas_filter` isolates the run.

## Notes

- Profiles are regenerated from the BAMs: the filter needs the nine-column
  profile format introduced with it, so `input.reuse_from` stays empty.
- The metadata sheet is the corrected copy from the revision run
  (`metadata_md_bam_fixed.tsv`; 90 samples, 15 per group × timepoint).
- To reproduce later, redo steps 1 to 3 with the same commit; nothing in
  `~/AlleleFlux-dev` is part of this run.

## Follow-up: `site_ratio` 10 run (2026-10-27)

Same frozen code (`48e3db6`), same env (`alleleflux-run-dm`), same inputs and
tests; the only analytical change is `site_ratio: 10.0` instead of `3.0`.

| what | where |
|---|---|
| config | `alleleflux_config_midas_filter_ratio10.yaml` |
| driver | `slurm_scripts/run_alleleflux_midas_filter_ratio10.sh` (job 16324300, submitted 2026-10-27; a first submission, 16323991, died at start-up on an unbound `PS1` under `set -u`, fixed in the driver); FINISHED 2026-10-30 09:11, 1899/1899 steps, driver wall time 2 d 20 h. 157 child jobs failed and were retried successfully: all CMH tests OOM-killed at the shared 50G `statistical_tests` allowance (second attempt 100G passed for most, 16 needed 150G); the ratio-3 run showed the same 155-failure pattern. Verified after the run: every cache job's command carried `--site_ratio 10.0`; the filter removed 0.097% of covered positions vs 0.994% at ratio 3 (same 4584 sample x MAG rows, same covered total) |
| working directory | `/scratch/gpfs/AMOELLER/sidd/diet_manip/midas_filter_Oct_2026/AlleleFlux_midas_filter_ratio10` |
| WCKey | `alleleflux_AlleleFlux_midas_filter_ratio10` |

Profiles are **reused** from the ratio-3 run through `input.profiles_path`
(the 90 per-sample directories with their `.profile_done` sentinels), because
the filter is applied when the allele-frequency cache is written, not during
profiling. `input.reuse_from` is deliberately left empty: it would also reuse
the ratio-3 QC and cache, and the new ratio would never take effect. The dry
run from the frozen env listed 4 jobs before the eligibility checkpoint
(metadata, QC, eligibility, all) against 97 for the ratio-3 run, i.e. all 90
profiling jobs dropped out.
