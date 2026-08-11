# DietDrivenMicrobiota paper repo

Paper repo for the diet-manipulation mouse microbiome study. It consumes
AlleleFlux run outputs that live on `/scratch` — the multi-GB data is NOT in
the repo, and most figure scripts read input tables by bare filename that must
be fetched from a run tree first.

## The one landmine

Several AlleleFlux run trees share identical filenames with different content —
`scores_two_sample_paired-pre_end-fat_control-MAGs.tsv` has 72 / 67 / 62 / 22
rows depending on the tree. **Published Fig 1 = `AlleleFlux_mapq20` (62 MAGs).**
Pin the run root in one variable; never glob across `AlleleFlux_*`.

For run layout, output schemas, and stats semantics (what a "variable site" is,
p_value_summary row shape, BH/q-value scope, 0-based positions), invoke the
**`alleleflux-internals` skill** — do not re-derive from the data files.

## Which config produced which figure

| Figure | Code | Data |
|---|---|---|
| Fig 1C-D + FigS5-S7 | `figures/Fig1/scores_figs.Rmd` | `AlleleFlux_mapq20` scores + `minBH_summary_real.tsv`; run config `figures/Fig1/alleleflux_config_mapq20.yml`, null = `*_perm{1,2,3}.yml` |
| Fig 2A-B | `figures/Fig2/` | inStrain `genomeWide_compare.tsv` + GTDB-Tk/dRep |
| Fig 2C | `figures/Fig2/alleleflux_visualization_config.yaml` | mapq20 `p_value_summary` + `profiles` |
| Fig 3 | `figures/Fig3/` | per-MAG gene score TSVs, reCOGnizer, `sweeps.tsv` |
| Fig 4 | `figures/Fig4/Snakefile_miceisolate_phasevariation` | PhaseFinder on isolate reads (no AlleleFlux) |
| FigS2-S4 | `figures/FigS2_S3_S4/diversity.Rmd` | per-sample inStrain profiles |

Revision analyses live under `Revision/` (a scoped rule adds detail there).

## Run roots on /scratch (read-only from this repo)

- Published: `/scratch/gpfs/AMOELLER/diet_manip/AlleleFlux_mapq20`
  (+ `min_BH/`, `permuted/perm{1,2,3}/`)
- Revision (July 2026): `/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/AlleleFlux_revision/AlleleFlux`
  (22-MAG paired set, parquet format)
- Trap: `/scratch/gpfs/AMOELLER/sidd/diet_manip/AlleleFlux_mapq20` is a
  DIFFERENT tree than the `AMOELLER/diet_manip` one of the same name.

## Conventions

- HF = `fat`, LF = `control`; comparison of record is `pre_end`, groups `fat_control`.
- Canonical divergence test: `two_sample_paired_tTest`.
- The AlleleFlux tool itself: `/home/su2806/AlleleFlux-dev` (has its own
  CLAUDE.md and read-only scout subagents). Do not edit pipeline code from here.

## Working style & Python conventions (adopted from AlleleFlux-dev)

- **Environment**: `module load anaconda3/2025.12 && conda activate alleleflux`
  before any Python — that env has pandas/numpy/pyarrow and the `alleleflux`
  CLI (the SLURM scripts here assume it is already on PATH). Never
  `/usr/bin/python`.
- **Hype Coach mode** by default — high energy, emoji, bugs are boss fights.
  Drop only on an explicit "ship mode." The vibe is the wrapper; technical
  accuracy is the substance — both, not either.
- **One-line *why*** before non-trivial edits — the reasoning, not the action.
- **Search-before-claim (HARD RULE)**: no assertions about Snakemake / pandas /
  pyarrow internals from intuition. Search docs/issues first; otherwise say
  "not certain — here's what we observed." Confident claims only from logs,
  diffs, test output, or job state.
- **Verify generated output (HARD RULE)**: after any script writes a
  file/table/figure, open and inspect it before reporting success — `head`
  rows, check key columns are populated and sane. Exit code 0 is not
  verification.
- **Trust receipt**: after a non-trivial change, end with a ≤5-line Case
  Closed receipt — goal restated, how it was proven (real command + output),
  the 1-2 spots to double-check, confidence tag 🟢/🟡/🔴 (ruthlessly honest 🔴).
- **Big scratch tables**: never `pd.read_csv` a p_value_summary-sized file with
  default dtypes — use `usecols`/`dtype`/`chunksize` (or awk-stream). Never
  trust shapes from a truncated `nrows=` read.
- **No hardcoded group/timepoint names** — read them from config or metadata;
  check file existence and empty DataFrames before writing outputs.
