# *Phocaeicola sartorii* isolate ANI, PRE vs END

Did the same *P. sartorii* strain persist in a mouse from pre-treatment to the
end of the experiment, and does that hold under both diets? Answered with
average nucleotide identity (ANI) between isolate genome assemblies from the
Moeller Lab isolate library. Not part of the paper.

## Files

- **`copy_sartorii_genomes.py`** — pulls the *P. sartorii* rows out of the isolate
  library spreadsheet, copies their assemblies out of the shared folder, and
  writes `sartorii_metadata.tsv` (isolate, mouse, diet, timepoint) plus a list of
  any isolates with no assembly on disk.
- **`Snakefile`** — the workflow, from assemblies to annotated pair tables.
- **`annotate_ani.py`** — attaches mouse metadata to every genome pair and labels
  it, so the within-mouse PRE-vs-END pairs can be judged against the
  between-mouse pairs rather than against a fixed ANI threshold.

## Workflow

1. **Copy assemblies and curate metadata.** One diet label is corrected
   (`DIET_FIXES` in the Snakefile). Only pre-treatment and end isolates are used.
2. **Pair table.** Counts of PRE × END isolates per mouse, before any filtering.
3. **Assembly QC.** seqkit stats (report only) and CheckM2. Genomes are kept at
   the MIMAG high-quality bar: > 90 % complete, < 5 % contamination. The pair
   table is recomputed after the filter.
4. **ANI.** skani, at two sensitivities (`--slow` as the primary setting for these
   fragmented assemblies, default as a sensitivity check), and fastANI at two
   fragment lengths (3000 and 1000 bp). fastANI is included because skani reports
   ANI to two decimals and every comparison here saturates at 99.98 to 100.00;
   fastANI resolves differences skani cannot. Each mode is run and annotated
   identically so the modes are directly comparable.
5. **Annotate and classify.** Every pair is labelled `within_host_across_time`
   (the question), `within_host_same_time`, `between_host_same_diet` or
   `between_host_diff_diet`. The between-host pairs are the null: if unrelated
   isolates score as high as within-mouse PRE-vs-END pairs, high ANI proves
   nothing. The script logs the within-host and between-host medians and ranges.
6. **Summary table** per mode, with the QC bounds of the genomes that went in.

Paths at the top of the Snakefile and scripts refer to the workstation where
this was run (`/workdir1/sidd/...`); run with `--use-conda` so the CheckM2 rule
gets its environment.
