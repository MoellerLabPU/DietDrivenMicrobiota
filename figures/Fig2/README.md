# Fig 2 — Strain Replacement and Allele Frequency Visualization

## Fig 2A-B — Strain Replacement Analysis

Panels A and B assess whether diet-driven changes are due to strain replacement using popANI clustering from inStrain.

- **`strain_replacement.Rmd`** — R Markdown that builds popANI distance matrices per species-group bin, performs permutation t-tests to test within-diet clustering vs. across-diet comparisons, and applies BH correction.

### Required Inputs

- `instrainComparer_breadth0.05_new_genomeWide_compare.tsv` — inStrain genome-wide comparison output
- `gtdbtk_representative_taxonomy.txt` — GTDB-Tk taxonomy with dRep secondary cluster assignments
- Mouse metadata table

## Fig 2C — Allele Frequency Trajectories

Panel C visualizes allele frequency trajectories at significant sites using the [AlleleFlux visualization workflow](https://github.com/MoellerLabPU/AlleleFlux).

- **`alleleflux_visualization_config.yaml`** — Example configuration file for the AlleleFlux visualization Snakemake workflow.

### Required Inputs

- AlleleFlux significance test outputs (p-value summary tables)
- Nucleotide frequency profiles
- Sample metadata
