# Fig 2 — Strain Replacement and Allele Frequency Visualization

## Fig 2A-B — Strain Replacement Analysis

Panels A and B assess whether diet-driven changes are due to strain replacement using popANI clustering from inStrain.

- **`strain_replacement.Rmd`** — R Markdown that builds popANI distance matrices per species-group bin, performs permutation t-tests to test within-diet clustering vs. across-diet comparisons, and applies BH correction.

### Required Inputs

- `instrainComparer_breadth0.05_new_genomeWide_compare.tsv` — inStrain genome-wide comparison output
- `gtdbtk_representative_taxonomy.txt` — GTDB-Tk taxonomy with dRep secondary cluster assignments
- Mouse metadata table

## Fig 2D — Per-litter Anchor-allele Frequency Violins

Panel D shows, for *Phocaeicola sartorii* (`SLG443_DASTool_bins_SLG443_bin.96`), the distribution of anchor-allele frequencies across the top 2,000 divergent sites, one violin per litter, in a 2 × 3 grid (HF / LF × PRE / POST / END). Each value is the mean frequency across a litter's mice for that diet and timepoint; violin width scales with the number of sites behind it (`density_norm="count"`).

- **`per_litter_violins.ipynb`** — Jupyter notebook (pandas + seaborn) that streams the TruSeq rows out of the per-site long frequency table written by the visualization workflow's `track_freqs` step, averages each litter's mice per site, and draws the figure.

### Required Inputs

- `<MAG>_frequency_table.long.tsv` — long-format allele frequency table from the AlleleFlux visualization workflow (`track_freqs/`)

## Fig 2E — Allele Frequency Trajectories, Five Mice Excluded

All three files live in `Fig2E_excl5mice/`. Panel E repeats the Fig 2C trajectory plot for `SLG443_DASTool_bins_SLG443_bin.96` with mice 534, 538, 539, 540 and 541 removed from the Hackflex metadata. The anchor alleles are unchanged, since those mice have no rows in the TruSeq metadata the anchor step uses.

- **`Fig2E_excl5mice/make_filtered_hackflex_metadata.py`** — Drops the five mice from the Hackflex metadata sheet and writes the filtered sheet the config points at. Run once before the workflow.
- **`Fig2E_excl5mice/alleleflux_visualization_config_excl5mice.yaml`** — Configuration for the visualization workflow: the filtered Hackflex sheet, the single MAG, one combined line plot per replicate, and bin widths of 10, 20, 30, 40 and 50 days in a single run. The paper uses the **50-day** bin width; the other widths were produced for comparison and the config is kept as run.
- **`Fig2E_excl5mice/run_viz_excl5mice.sh`** — SLURM submission script for the run (paths are those of the original run and are cluster-specific).

### Required Inputs

- Hackflex and TruSeq sample metadata sheets
- AlleleFlux `two_sample_paired_tTest` p-value summary for `pre_end`, `fat_control`
- Nucleotide frequency profiles for both library types
