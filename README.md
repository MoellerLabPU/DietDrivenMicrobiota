# DietDrivenMicrobiota

Analysis code and workflows for studying diet-driven allele frequency changes in gut microbiota using [AlleleFlux](https://github.com/MoellerLabPU/AlleleFlux).

This repository contains the computational analysis pipeline for investigating how dietary interventions (high-fat vs. control diet) shape within-species genetic variation in the mouse gut microbiome. Using longitudinal metagenomic data, we track allele frequency shifts in metagenome-assembled genomes (MAGs) across experimental timepoints and perform statistical testing to identify positions under diet-driven selection.

## Repository Structure

```
DietDrivenMicrobiota/
│
├── scripts/                           # Helper Python scripts
│   ├── bh_list_sam.py                 # BH-corrected p-value summary across tests
│   ├── combine_files.py               # Concatenate per-MAG significance files
│   └── extract_mag_positions.py       # Extract (MAG, contig, position) tuples
│
├── notebooks/                         # Jupyter notebooks for analysis
│   ├── add_protein_description.ipynb  # Merge p-value summaries with protein annotations
│   └── isolate_analysis.ipynb         # Isolate-level analysis and visualization
│
├── configs/                           # Example configuration files
│   └── visualization_config.yaml      # AlleleFlux visualization workflow config
│
├── megaTable/                         # MegaTable Snakemake workflow
│   ├── README.md                      # Detailed workflow documentation
│   ├── MEGATABLE_COLUMN_DOCUMENTATION.md  # Column reference for output tables
│   ├── create_megatable.smk           # Snakemake workflow
│   ├── config.yml                     # Workflow configuration
│   ├── config_mapq2.yml               # Config for MAPQ ≥ 2 analysis
│   ├── config_mapq20.yml              # Config for MAPQ ≥ 20 analysis
│   └── profile/                       # SLURM cluster profile
│
├── isolate_Bacteroides_muris/         # Bacteroides muris isolate analysis
│   ├── qc.smk                         # QC workflow (FastQC → fastp → MultiQC)
│   ├── align.smk                      # Alignment workflow (Bowtie2 → SAMtools)
│   ├── copy_reads.py                  # Subset & copy isolate reads
│   ├── prepare_metadat.py             # Generate AlleleFlux metadata from isolate sheet
│   ├── alleleflux_config.yaml         # AlleleFlux config for isolate analysis
│   └── slurm_profile/                 # SLURM cluster profile
│
├── LICENSE                            # GNU General Public License v3.0
└── README.md                          # This file
```

## Components

### Helper Scripts (`scripts/`)

| Script | Description |
|--------|-------------|
| `bh_list_sam.py` | Computes Benjamini–Hochberg-corrected p-value summaries across multiple statistical tests (paired/unpaired t-tests, Wilcoxon, LMM, across-time models) for two comparison periods. Outputs both wide summary tables and detailed significant-row lists. |
| `combine_files.py` | Concatenates per-MAG AlleleFlux significance test result files (paired-sample and single-sample) into unified tables, ensuring MAG ID consistency across file types. |
| `extract_mag_positions.py` | Extracts and combines unique (MAG, contig, position) tuples from per-MAG AlleleFlux result files. Supports both `two_sample_paired` and `single_sample` test types with group-level partitioning. |

### MegaTable Workflow (`megaTable/`)

A Snakemake workflow that consolidates AlleleFlux analysis outputs into comprehensive "megatables" combining:

- BH-corrected p-values from multiple statistical tests
- Per-position coverage and allele frequency statistics
- Quality control metrics per MAG

See [`megaTable/README.md`](megaTable/README.md) for full documentation, and [`megaTable/MEGATABLE_COLUMN_DOCUMENTATION.md`](megaTable/MEGATABLE_COLUMN_DOCUMENTATION.md) for a complete column reference.

### Isolate Analysis (`isolate_Bacteroides_muris/`)

Pipeline for processing *Bacteroides muris* isolate sequencing data:

1. **Read subsetting** — Identify and copy reads for *B. muris* isolates from the strain library
2. **Quality control** — FastQC → fastp trimming → FastQC → MultiQC
3. **Alignment** — Bowtie2 alignment → SAMtools BAM conversion, sorting, and indexing
4. **AlleleFlux analysis** — Allele frequency profiling and significance testing on isolate data

### Notebooks (`notebooks/`)

- **`add_protein_description.ipynb`** — Merges p-value summary tables with Prodigal protein annotations for functional context
- **`isolate_analysis.ipynb`** — Exploratory analysis and visualization of isolate-level AlleleFlux results

### Visualization Config (`configs/`)

- **`visualization_config.yaml`** — Example configuration for the AlleleFlux visualization workflow, including terminal nucleotide analysis, allele frequency tracking, and plotting parameters
