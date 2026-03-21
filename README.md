# DietDrivenMicrobiota

Analysis code and workflows for studying diet-driven allele frequency changes in gut microbiota.

This repository contains the computational analysis pipeline for investigating how dietary interventions (high-fat vs. control diet) shape within-species genetic variation in the mouse gut microbiome. Using longitudinal metagenomic data, we track allele frequency shifts in metagenome-assembled genomes (MAGs) across experimental timepoints and perform statistical testing to identify positions under diet-driven selection.

## Repository Structure

```
DietDrivenMicrobiota/
│
├── figures/                               # Code for generating paper figures
│   ├── Fig1/                              # Fig 1C-D: Phylogeny + score heatmaps
│   ├── Fig2/                              # Fig 2: Strain replacement & allele freq visualization
│   ├── Fig3/                              # Fig 3: Gene-level functional enrichment
│   ├── Fig4/                              # Fig 4: Phase variation (B. muris)
│   ├── FigS2_S3_S4/                       # Fig S2-S4: Community diversity analysis
│   └── FigS5_S6_S7/                       # Fig S5-S7: Supplementary score comparisons
│
├── miscellaneous scripts/                 # Helper Python scripts
│   ├── bh_list_sam.py                     # BH-corrected p-value summary across tests
│   └── combine_files.py                   # Concatenate per-MAG significance files
│
├── notebooks/                             # Jupyter notebooks for analysis
│   ├── add_protein_description.ipynb      # Merge p-value summaries with protein annotations
│   ├── isolate_analysis.ipynb             # Isolate-level analysis and visualization
│   └── tested_sites.ipynb                 # Tested sites analysis
│
├── megaTable/                             # MegaTable Snakemake workflow
│   ├── README.md                          # Detailed workflow documentation
│   ├── MEGATABLE_COLUMN_DOCUMENTATION.md  # Column reference for output tables
│   ├── create_megatable.smk               # Snakemake workflow
│   ├── extract_mag_positions.py           # Extract (MAG, contig, position) tuples
│   └── config*.yml                        # Workflow configurations
│
├── isolate_Bacteroides_muris/             # B. muris isolate processing pipeline
│   ├── qc.smk                             # QC workflow (FastQC → fastp → MultiQC)
│   ├── align.smk                          # Alignment workflow (Bowtie2 → SAMtools)
│   ├── copy_reads.py                      # Subset & copy isolate reads
│   ├── prepare_metadat.py                 # Generate AlleleFlux metadata
│   └── alleleflux_config.yaml             # AlleleFlux config for isolate analysis
│
├── processing_MAGs/                       # MAG processing commands
│   ├── QC.txt                             # Quality control (GUNC, CheckM2)
│   └── MAG_processing.txt                 # Taxonomy, phylogeny, dereplication, mapping
│
├── LICENSE                                # GNU General Public License v3.0
└── README.md                              # This file
```

## Figures

Each figure's code is organized in its own subdirectory under [`figures/`](figures/). See the [figures README](figures/README.md) for a quick reference table.

| Figure | Directory | Description |
|--------|-----------|-------------|
| Fig 1C-D | [`figures/Fig1/`](figures/Fig1/) | Phylogeny with divergence/parallelism score heatmaps |
| Fig 2A-C | [`figures/Fig2/`](figures/Fig2/) | Strain replacement (popANI) and allele frequency trajectories |
| Fig 3 | [`figures/Fig3/`](figures/Fig3/) | Gene-level COG functional enrichment (hypergeometric tests) |
| Fig 4 | [`figures/Fig4/`](figures/Fig4/) | Phase variation analysis using PhaseFinder |
| Fig S2-S4 | [`figures/FigS2_S3_S4/`](figures/FigS2_S3_S4/) | Community diversity (relative abundance, Shannon, ANCOM-BC2) |
| Fig S5-S7 | [`figures/FigS5_S6_S7/`](figures/FigS5_S6_S7/) | Supplementary score comparisons (uses code from Fig 1) |

## Components

### Miscellaneous Scripts (`miscellaneous scripts/`)

| Script | Description |
|--------|-------------|
| `bh_list_sam.py` | Computes BH-corrected p-value summaries across multiple statistical tests for two comparison periods |
| `combine_files.py` | Concatenates per-MAG AlleleFlux significance test result files into unified tables |

### MegaTable Workflow (`megaTable/`)

Snakemake workflow that consolidates AlleleFlux analysis outputs into comprehensive summary tables combining BH-corrected p-values, per-position coverage and allele frequency statistics, and quality control metrics. See [`megaTable/README.md`](megaTable/README.md) for full documentation.

### Isolate Analysis (`isolate_Bacteroides_muris/`)

Pipeline for processing *Bacteroides muris* isolate sequencing data: read subsetting, quality control (FastQC → fastp → MultiQC), alignment (Bowtie2 → SAMtools), and AlleleFlux allele frequency profiling.

### MAG Processing (`processing_MAGs/`)

Reference commands for MAG quality control (GUNC, CheckM2), taxonomy assignment (GTDB-Tk), phylogenetic tree construction (IQ-TREE), dereplication (dRep), read mapping (Bowtie2), gene prediction (Prodigal), functional annotation (reCOGnizer), and strain-level profiling (inStrain).

## License

This project is licensed under the GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.
