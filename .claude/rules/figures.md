---
paths:
  - "figures/**"
---

# Figure-code gotchas (verified 2026-08-10)

## Fig 1 (`figures/Fig1/scores_figs.Rmd`)

- The Rmd reads every input table by **bare filename** (no directory); none of
  those files are committed. Copy them from
  `/scratch/gpfs/AMOELLER/diet_manip/AlleleFlux_mapq20/longitudinal/scores/processed/combined/MAG/`
  (score TSVs) and `.../AlleleFlux_mapq20/min_BH/` (`minBH_summary_real.tsv`)
  into the working dir before knitting.
- The MAG set is gated on `two_samp_pre_end` — used 19×, **assigned nowhere**
  in the repo (`scores_figs.Rmd:39` first use). It is almost certainly the
  62-row `scores_two_sample_paired-pre_end-fat_control-MAGs.tsv` (read into
  `div_mapq20` at lines 5/200); that identification is inferred, not stated.
- `mouse_MAGs.contree` (the Fig-1 tree) exists **nowhere on this filesystem**
  (scoping 2026-08) — Fig 1 cannot currently be re-rendered from scratch.
- `alleleflux_config_mapq20.yml:10` has a **broken** `fasta_path` (missing
  `sg4230/`). The perm configs use the working path. Working copies:
  `/scratch/gpfs/AMOELLER/diet_manip/redo_representative_megamag.fa` and
  `.../copy_sg4230_scratch/sg4230/popgentoolkit/redo_representative_megamag.fa`.
- The three perm configs differ from the base in FOUR ways: `run_name`,
  `metadata_path` (→ `permuted/permute_diet_every_other_mouse_md/perm_group_swap_set{1,2,3}.tsv`),
  `output.root_dir` (→ `permuted/perm{1,2,3}`), and the fixed `fasta_path`.
  They are full reruns on hand-made swap sheets — NOT seed-generated
  `permutation:`-block runs (that newer mechanism is used in `Revision/`).

## Fig 2C (`figures/Fig2/alleleflux_visualization_config.yaml`)

- `test_type: "two_sample_paired_tTest"` (line 54) and
  `significant_sites_file` → mapq20 `p_value_summary_two_sample_paired_pre_end.tsv`
  (line 40). This config is the base that
  `Revision/AlleleFlux/alleleflux_visualization_config_per_replicate.yaml`
  was copied from.

Schemas / stats semantics for anything these figures read: invoke the
`alleleflux-internals` skill.
