# MegaTable Workflow

## Overview

The MegaTable workflow consolidates [AlleleFlux](https://github.com/MoellerLabPU/AlleleFlux) analysis results into comprehensive tables for downstream analysis. It combines statistical test results (p-values), coverage statistics, and quality control metrics into unified "megatables" for both paired and single-sample analyses.

The workflow processes allele frequency data from metagenomic samples, filtering for statistically significant positions (q-value < 0.05) and enriching them with coverage statistics, allele frequencies, and QC metrics.

## Workflow Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         INPUT DATA SOURCES                                  │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Significance directories (two_sample_paired & single_sample)             │
│  • P-value summary files (CMH, LMM, paired/unpaired t-tests)                │
│  • MAG metadata files (per-sample profiles with allele counts)              │
│  • Reference FASTA & MAG-to-contig mapping                                  │
│  • QC directory (for sample filtering by breadth threshold)                 │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                         STEP 1: EXTRACT POSITIONS                           │
│  Rule: extract_mag_positions                                                │
│  Script: extract_mag_positions.py                                           │
├─────────────────────────────────────────────────────────────────────────────┤
│  Extracts unique (MAG, contig, position) tuples from significance tests     │
│  → tested_positions/two_sample_paired_mag_positions.tsv.gz                  │
│  → tested_positions/single_sample_{group}_mag_positions.tsv.gz              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STEP 2: FILTER & PIVOT P-VALUES                          │
│  Rules: bh_pValues_paired / bh_pValues_single                               │
│  (Inline Python using _process_and_pivot_pvalue_files helper)               │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Filters positions with q_value < 0.05 from primary test                  │
│  • Merges p-values from multiple tests (CMH, LMM, paired/unpaired)          │
│  • Pivots test results to wide format: min_p_value_{test}, q_value_{test}   │
│  → bh_pValues/two_sample_paired_bh_pValues.tsv.gz                           │
│  → bh_pValues/single_sample_{group}_bh_pValues.tsv.gz                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                      STEP 3: QUALITY CONTROL                                │
│  Rules: quality_control_paired / quality_control_single                     │
│  CLI: alleleflux-qc-positions                                               │
├─────────────────────────────────────────────────────────────────────────────┤
│  Runs position-based QC, generates three summary files:                     │
│  • ALL_MAGs_positions_QC_overall_summary.tsv (breadth, coverage per MAG)    │
│  • ALL_MAGs_positions_QC_group_summary.tsv (metrics stratified by group)    │
│  • ALL_MAGs_positions_QC_group_time_summary.tsv (group×time, paired only)   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                   STEP 4: COMPUTE MAG STATISTICS                            │
│  Rule: compute_stats (parallelized per MAG)                                 │
│  CLI: alleleflux-coverage-allele-stats                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│  Calculates per-position coverage and allele frequency statistics:          │
│  • Mean/std coverage (present-only & with-zeros variants)                   │
│  • Allele frequencies (A, C, G, T) with mean/std                            │
│  • Group-stratified statistics (if group/time metadata available)           │
│  → coverage_stats/{mag}_stats.tsv (one per MAG)                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                STEP 5: MERGE P-VALUES WITH COVERAGE STATS                   │
│  Rules: create_merged_table_paired / create_merged_table_single             │
│  (Inline Python using _merge_coverage_with_pvalues helper)                  │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Loads MAG stats for MAGs present in p-value tables                       │
│  • Parallel processing via multiprocessing.Pool                             │
│  • Merges on (MAG, contig, position) keys                                   │
│  → pValue_stats_merged/two_sample_paired_merged_table.tsv.gz                │
│  → pValue_stats_merged/single_sample_{group}_merged_table.tsv.gz            │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                  STEP 6: ADD QC SUMMARIES (FINAL MEGATABLES)                │
│  Rules: create_megatable_paired / create_megatable_single                   │
│  (Inline Python using _load_and_flatten_qc_summaries helper)                │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Flattens QC summaries (group/time metrics become separate columns)       │
│    - Example: breadth_mean_group_fat, breadth_mean_group_fat_time_pre       │
│  • Merges QC data with merged tables on MAG_ID                              │
│  → megatables/two_sample_paired_megatable.tsv.gz                            │
│  → megatables/single_sample_{group}_megatable.tsv.gz                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
                         ┌────────────────────┐
                         │  FINAL MEGATABLES  │
                         │  (Comprehensive)   │
                         └────────────────────┘
```

## Input Data Requirements

### Input File Formats

| File Type | Format | Required Columns | Source |
|-----------|--------|------------------|--------|
| **Metadata TSV** | Tab-separated | `sample_id, file_path` (optional: `group, time`) | AlleleFlux Step 1 output |
| **MAG Mapping** | Tab-separated | `mag_id, contig_id` | Created by `alleleflux-create-mag-mapping` |
| **Profile Files** | Tab-separated | `contig, position, total_coverage, A, C, G, T` | AlleleFlux profiling |
| **Significance Results** | Gzipped TSV | `contig, position, gene, min_p_value, q_value` | AlleleFlux Step 2 statistical tests |
| **P-value Summary** | Tab-separated | `MAG, contig, position, gene, min_p_value, q_value` | AlleleFlux Step 2 summary files |
| **Positions File** | Gzipped TSV | `MAG, contig, position` | Generated by Step 1 of this workflow |

### Required AlleleFlux Outputs

This workflow consumes outputs from a completed AlleleFlux analysis run:

```
AlleleFlux_output/
├── significance_tests/
│   ├── two_sample_paired_{timepoints}/     # Per-MAG .tsv.gz files
│   │   └── {MAG}_two_sample_paired.tsv.gz
│   └── single_sample_{timepoints}/         # Per-MAG per-group files
│       └── {MAG}_{group}_single_sample.tsv.gz
├── p_value_summaries/
│   ├── p_value_summary_cmh_{timepoints}.tsv
│   ├── p_value_summary_two_sample_paired_{timepoints}.tsv
│   ├── p_value_summary_two_sample_unpaired_{timepoints}.tsv
│   ├── p_value_summary_lmm_{timepoints}.tsv
│   ├── p_value_summary_single_sample_{timepoints}.tsv
│   ├── p_value_summary_lmm_across_time_{timepoints}.tsv
│   └── p_value_summary_cmh_across_time_{timepoints}.tsv
├── metadata/
│   └── {MAG}_metadata.tsv                  # Per-sample file paths and group/time info
├── QC/
│   └── {MAG}_QC.tsv                        # Per-sample breadth and coverage
└── reference.fasta                         # Concatenated MAG sequences
```

## Output Structure

### Final Megatables Column Reference

The megatables contain all information needed for downstream analysis:

#### 1. Position Identifiers

| Column | Description |
|--------|-------------|
| `MAG` | MAG identifier |
| `contig` | Contig name within the MAG |
| `position` | 1-based nucleotide position |
| `gene` | Gene ID if position falls within a gene (from Prodigal annotations) |

#### 2. P-values & Q-values (BH-corrected)

**For Paired Analysis (`two_sample_paired_megatable.tsv.gz`):**

| Column | Description |
|--------|-------------|
| `min_p_value_CMH` | Minimum p-value from Cochran-Mantel-Haenszel test |
| `q_value_CMH` | BH-corrected q-value for CMH |
| `min_p_value_two_sample_paired_tTest` | Minimum p-value from paired t-test |
| `q_value_two_sample_paired_tTest` | BH-corrected q-value for paired t-test |
| `min_p_value_LMM` | Minimum p-value from Linear Mixed Model |
| `q_value_LMM` | BH-corrected q-value for LMM |
| `min_p_value_two_sample_unpaired_tTest` | Minimum p-value from unpaired t-test |
| `q_value_two_sample_unpaired_tTest` | BH-corrected q-value for unpaired t-test |

**For Single-Sample Analysis (`single_sample_{group}_megatable.tsv.gz`):**

| Column | Description |
|--------|-------------|
| `min_p_value_single_sample_tTest_{group}` | P-value from single-sample t-test for group |
| `q_value_single_sample_tTest_{group}` | BH-corrected q-value |
| `min_p_value_LMM_across_time` | P-value from LMM across time |
| `q_value_LMM_across_time` | BH-corrected q-value |
| `min_p_value_CMH_across_time` | P-value from CMH across time |
| `q_value_CMH_across_time` | BH-corrected q-value |

#### 3. Coverage Statistics (per position)

| Column | Description |
|--------|-------------|
| `total_coverage` | Sum of coverage across all samples |
| `mean_coverage` | Mean coverage (samples with coverage > 0 only) |
| `std_coverage` | Standard deviation of coverage (present samples) |
| `mean_coverage_with_zeros` | Mean coverage including samples with 0 coverage |
| `std_coverage_with_zeros` | Std dev including samples with 0 coverage |
| `n_present` | Number of samples with coverage > 0 at this position |
| `n_samples` | Total number of samples |

**Group-stratified Coverage:**

| Column Pattern | Description |
|----------------|-------------|
| `mean_coverage_group_{group}` | Mean coverage for samples in group |
| `std_coverage_group_{group}` | Std dev of coverage for samples in group |
| `mean_coverage_with_zeros_group_{group}` | Mean coverage including zeros for group |

**Time-stratified Coverage (longitudinal data only):**

| Column Pattern | Description |
|----------------|-------------|
| `mean_coverage_group_{group}_time_{time}` | Mean coverage for group at timepoint |
| `std_coverage_group_{group}_time_{time}` | Std dev for group at timepoint |

#### 4. Allele Frequencies (A, C, G, T)

| Column Pattern | Description |
|----------------|-------------|
| `mean_freq_{nucleotide}` | Mean allele frequency across samples |
| `std_freq_{nucleotide}` | Standard deviation of allele frequency |
| `mean_freq_{nucleotide}_group_{group}` | Mean frequency for group |
| `std_freq_{nucleotide}_group_{group}` | Std dev of frequency for group |
| `mean_freq_{nucleotide}_group_{group}_time_{time}` | Frequency at group×time |

*Where `{nucleotide}` is one of: A, C, G, T*

#### 5. QC Metrics (per MAG)

| Column | Description |
|--------|-------------|
| `num_samples` | Total samples with data for this MAG |
| `breadth_mean` | Mean genome breadth (fraction covered) across samples |
| `breadth_std` | Std dev of genome breadth |
| `average_coverage_mean` | Mean of mean coverage across samples |
| `median_coverage_mean` | Mean of median coverage across samples |

**Group-stratified QC:**

| Column Pattern | Description |
|----------------|-------------|
| `breadth_mean_group_{group}` | Mean breadth for samples in group |
| `average_coverage_mean_group_{group}` | Mean coverage for group |

**Time-stratified QC (paired data only):**

| Column Pattern | Description |
|----------------|-------------|
| `breadth_mean_group_{group}_time_{time}` | Mean breadth at group×time |
| `average_coverage_mean_group_{group}_time_{time}` | Mean coverage at group×time |

## Configuration

### Required Config Parameters (`config.yml`)

```yaml
# Run identifier for logging
run_name: "my_analysis"

# Output directory
output_dir: "/path/to/output"

# Input directories (from AlleleFlux analysis)
significance_dir_two_sample_paired: "/path/to/significance/two_sample_paired_{timepoints}"
significance_dir_single_sample: "/path/to/significance/single_sample_{timepoints}"
p_value_summary_dir: "/path/to/p_value_summaries"
metadata_dir: "/path/to/mag_metadata"
qc_dir: "/path/to/qc_results"

# Reference files
fasta: "/path/to/reference.fasta"
mag_mapping: "/path/to/mag_to_contig_mapping.tsv"

# Analysis parameters
single_sample_groups: ["fat", "control"]  # List of groups for single-sample analysis
breadth_threshold: 0.5                     # Minimum breadth for sample inclusion
cpus: 16                                   # Parallel workers for compute_stats

# Optional
mag_list: "/path/to/mag_list.txt"  # Auto-generated if not provided
min_samples: 1                     # Minimum samples for MAG inclusion
```

### Configuration Parameters Explained

| Parameter | Required | Description |
|-----------|----------|-------------|
| `run_name` | Yes | Identifier for log file organization |
| `output_dir` | Yes | Base directory for all workflow outputs |
| `significance_dir_two_sample_paired` | Yes | Directory with per-MAG paired test `.tsv.gz` files |
| `significance_dir_single_sample` | Yes | Directory with per-MAG single-sample test files |
| `p_value_summary_dir` | Yes | Directory with `p_value_summary_*.tsv` files |
| `metadata_dir` | Yes | Directory with `{MAG}_metadata.tsv` files |
| `qc_dir` | Yes | Directory with `{MAG}_QC.tsv` files |
| `fasta` | Yes | Reference FASTA with concatenated MAG sequences |
| `mag_mapping` | Yes | TSV file mapping contig IDs to MAG IDs |
| `single_sample_groups` | Yes | List of group names for single-sample tests |
| `breadth_threshold` | Yes | Minimum genome breadth (0-1) for sample inclusion |
| `cpus` | Yes | Number of CPUs for parallel processing |
| `mag_list` | No | File with MAG IDs (one per line); auto-generated if missing |
| `min_samples` | No | Minimum samples required for MAG inclusion (default: 1) |

## Usage

### Running the Complete Workflow

```bash
cd /path/to/megaTable

# Dry run to check workflow
snakemake -s create_megatable.smk --configfile config.yml -n

# Execute workflow
snakemake -s create_megatable.smk --configfile config.yml --cores 16

# With cluster profile (recommended for large datasets)
snakemake -s create_megatable.smk --configfile config.yml --profile profile/
```

### Running Specific Steps

```bash
# Generate only BH p-value tables
snakemake -s create_megatable.smk --configfile config.yml \
    --cores 16 --until bh_pValues_paired bh_pValues_single

# Generate coverage stats only
snakemake -s create_megatable.smk --configfile config.yml \
   --cores 16 --until compute_stats

# Generate merged tables without QC summaries
snakemake -s create_megatable.smk --configfile config.yml \
    --cores 16 --until create_merged_table_paired create_merged_table_single
```

## Output Directory Structure

```
{output_dir}/
├── mag_list.txt                           # Auto-generated list of MAG IDs
├── tested_positions/
│   ├── two_sample_paired_mag_positions.tsv.gz
│   └── single_sample_{group}_mag_positions.tsv.gz
├── bh_pValues/
│   ├── two_sample_paired_bh_pValues.tsv.gz
│   └── single_sample_{group}_bh_pValues.tsv.gz
├── qc_results_paired/
│   ├── {MAG}_positions_QC.tsv             # Per-MAG per-sample QC
│   ├── ALL_MAGs_positions_QC_overall_summary.tsv
│   ├── ALL_MAGs_positions_QC_group_summary.tsv
│   └── ALL_MAGs_positions_QC_group_time_summary.tsv
├── qc_results_single_{group}/
│   ├── {MAG}_positions_QC.tsv
│   ├── ALL_MAGs_positions_QC_overall_summary.tsv
│   └── ALL_MAGs_positions_QC_group_summary.tsv
├── coverage_stats/
│   └── {MAG}_stats.tsv                    # Per-MAG coverage and allele stats
├── pValue_stats_merged/
│   ├── two_sample_paired_merged_table.tsv.gz
│   └── single_sample_{group}_merged_table.tsv.gz
└── megatables/
    ├── two_sample_paired_megatable.tsv.gz     # FINAL OUTPUT
    └── single_sample_{group}_megatable.tsv.gz # FINAL OUTPUT
```

## AlleleFlux CLI Commands Used

This workflow uses the following AlleleFlux command-line tools:

| Command | Script Source | Purpose |
|---------|---------------|---------|
| `alleleflux-coverage-allele-stats` | `alleleflux.scripts.accessory.coverage_and_allele_stats` | Compute per-position coverage and allele frequency statistics |
| `alleleflux-list-mags` | `alleleflux.scripts.accessory.list_mags` | Generate list of MAG IDs from metadata directory |
| `alleleflux-qc-positions` | `alleleflux.scripts.preprocessing.quality_control_positions` | Position-specific quality control |

### CLI Command Details

#### `alleleflux-coverage-allele-stats`

Computes per-position coverage and allele frequency statistics for a single MAG.

```bash
alleleflux-coverage-allele-stats \
    --rootDir /path/to/metadata_dir \
    --output_dir /path/to/coverage_stats \
    --mag_id MAG001 \
    --qc_dir /path/to/qc_results
```

| Argument | Description |
|----------|-------------|
| `--rootDir` | Directory containing `{MAG}_metadata.tsv` files |
| `--output_dir` | Output directory for stats files |
| `--mag_id` | Single MAG ID to process (memory-efficient mode) |
| `--qc_dir` | QC directory to filter samples by breadth threshold |

#### `alleleflux-qc-positions`

Runs position-based quality control, calculating breadth and coverage metrics.

```bash
alleleflux-qc-positions \
    --metadata_dir /path/to/metadata \
    --fasta /path/to/reference.fasta \
    --mag_mapping_file /path/to/mapping.tsv \
    --positions_file /path/to/positions.tsv.gz \
    --positions_denominator positions \
    --output_dir /path/to/qc_output \
    --breadth_threshold 0.5 \
    --coverage_threshold 0 \
    --cpus 16 \
    --data_type longitudinal
```

#### `alleleflux-list-mags`

Generates a list of MAG IDs from metadata files.

```bash
alleleflux-list-mags \
    --metadata_dir /path/to/metadata \
    --min_samples 1 \
    --output mag_list.txt
```

## Helper Functions

The workflow defines several helper functions for data processing:

| Function | Purpose |
|----------|---------|
| `_process_and_pivot_pvalue_files()` | Loads p-value files, filters by q-value threshold, pivots to wide format |
| `_load_and_flatten_qc_summaries()` | Loads QC summary files and flattens group/time columns |
| `_load_and_merge_single_mag()` | Worker for parallel processing: loads coverage data and merges with p-values |
| `_merge_coverage_with_pvalues()` | Orchestrates parallel merge of coverage stats with p-values |

## Logging

The workflow uses a custom logger (`workflow_logger`) that:

- Does not interfere with Snakemake's internal logging
- Formats numbers with thousand separators for readability
- Tracks progress through each step with informative messages

Example log output:

```
2025-10-08 14:30:15 [INFO] Loaded BH p-values with 1,234,567 positions
2025-10-08 14:30:20 [INFO] Found 273 unique MAGs in p-value table
Loading MAG coverage files: 100%|██████████| 273/273 MAGs [02:15<00:00, 2.02 MAG/s]
2025-10-08 14:32:35 [INFO] Loaded coverage stats with 45,678,901 total positions across all MAGs
2025-10-08 14:35:42 [INFO] Created merged_table with 1,234,567 positions and 156 columns
```

## Troubleshooting

### Common Issues

1. **Memory errors during merge operations**
   - Increase `mem_mb` in `resources` section of merge rules
   - Process fewer MAGs per batch by filtering input

2. **Missing coverage files**
   - Check that `compute_stats` completed for all MAGs
   - Verify MAG list file contains correct MAG IDs
   - Check log files in `logs/` directory

3. **QC summary files not found**
   - Ensure `quality_control_paired`/`quality_control_single` rules completed
   - Check that input positions files exist and are not empty

4. **Column name conflicts**
   - QC merge uses `_qc` suffix to avoid conflicts
   - Check for duplicate column names in intermediate files

## Example Usage

### Complete Workflow Example

```bash
# 1. Set up configuration
cd /path/to/megaTable
cp config.template.yml config.yml
# Edit config.yml with your paths

# 2. Dry run to verify DAG
snakemake -s create_megatable.smk --configfile config.yml -n

# 3. Run the full workflow
snakemake -s create_megatable.smk --configfile config.yml --cores 16

# 4. Check outputs
ls -la output_dir/megatables/
```

### Using Cluster Profile (SLURM)

```bash
# With cluster profile for large datasets
snakemake -s create_megatable.smk \
    --configfile config.yml \
    --profile profile/ \
    --jobs 100
```

### Rerunning Specific Steps

```bash
# Force rerun of coverage stats
snakemake -s create_megatable.smk --configfile config.yml \
    --cores 16 --forcerun compute_stats

# Rerun from merged tables onward
snakemake -s create_megatable.smk --configfile config.yml \
    --cores 16 --forcerun create_merged_table_paired create_merged_table_single
```

## Dependencies

- Python 3.8+
- pandas
- snakemake (≥7.0)
- tqdm
- AlleleFlux package (provides CLI commands: `alleleflux-coverage-allele-stats`, `alleleflux-list-mags`, `alleleflux-qc-positions`)

### Installing AlleleFlux

```bash
# From source (development mode)
cd /path/to/AlleleFlux
pip install -e .

# Or from conda
conda install -c bioconda alleleflux
```

## Workflow Rules Summary

| Rule | Step | Inputs | Outputs | Parallelization |
|------|------|--------|---------|-----------------|
| `extract_mag_positions` | 1 | Significance test directories | Position files | Sequential |
| `bh_pValues_paired` | 2 | P-value summary files | Filtered BH p-values | Sequential |
| `bh_pValues_single` | 2 | P-value summary files | Filtered BH p-values (per group) | Sequential |
| `quality_control_paired` | 3 | Positions + metadata | QC summaries | Multi-threaded |
| `quality_control_single` | 3 | Positions + metadata | QC summaries (per group) | Multi-threaded |
| `compute_stats` | 4 | Metadata + QC | Coverage stats (per MAG) | Per-MAG parallel |
| `create_merged_table_paired` | 5 | BH p-values + coverage stats | Merged table | Pool workers |
| `create_merged_table_single` | 5 | BH p-values + coverage stats | Merged table (per group) | Pool workers |
| `create_megatable_paired` | 6 | Merged table + QC summaries | Final megatable | Sequential |
| `create_megatable_single` | 6 | Merged table + QC summaries | Final megatable (per group) | Sequential |
