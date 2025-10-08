# MegaTable Workflow

## Overview

The MegaTable workflow consolidates AlleleFlux analysis results into comprehensive tables for downstream analysis. It combines statistical test results (p-values), coverage statistics, and quality control metrics into unified "megatables" for both paired and single-sample analyses.

## Workflow Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         INPUT DATA SOURCES                                  │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Significance directories (two_sample_paired & single_sample)             │
│  • P-value summary files (CMH, LMM, paired/unpaired t-tests)                │
│  • MAG metadata files (per-sample profiles)                                 │
│  • Reference FASTA & MAG mapping                                            │
│  • QC directory (optional, for filtering)                                   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                         STEP 1: EXTRACT POSITIONS                           │
│  extract_mag_positions                                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│  Extracts tested positions from significance test results                   │
│  → tested_positions/two_sample_paired_mag_positions.tsv.gz                  │
│  → tested_positions/single_sample_{group}_mag_positions.tsv.gz              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STEP 2: FILTER & PIVOT P-VALUES                          │
│  bh_pValues_paired / bh_pValues_single                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Filters positions with q_value < 0.05                                    │
│  • Merges p-values from multiple tests (CMH, LMM, paired/unpaired)          │
│  • Pivots test results to wide format with test-specific columns            │
│  → bh_pValues/two_sample_paired_bh_pValues.tsv.gz                           │
│  → bh_pValues/single_sample_{group}_bh_pValues.tsv.gz                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                      STEP 3: QUALITY CONTROL                                │
│  quality_control_paired / quality_control_single                            │
├─────────────────────────────────────────────────────────────────────────────┤
│  Runs QC on extracted positions, generates three summary files:             │
│  • ALL_MAGs_QC_overall_summary.tsv (breadth, coverage per MAG)              │
│  • ALL_MAGs_QC_group_summary.tsv (metrics stratified by group)              │
│  • ALL_MAGs_QC_group_time_summary.tsv (metrics by group×time, paired only)  │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                   STEP 4: COMPUTE MAG STATISTICS                            │
│  compute_stats (per MAG)                                                    │
├─────────────────────────────────────────────────────────────────────────────┤
│  Calculates per-position coverage and allele frequency statistics:          │
│  • Mean/std coverage (present-only & with-zeros)                            │
│  • Allele frequencies (A, C, G, T) with mean/std                            │
│  • Group-stratified statistics (if group/time metadata available)           │
│  → coverage_stats/{mag}_stats.tsv (one per MAG)                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                STEP 5: MERGE P-VALUES WITH MAG STATS                   │
│  create_merged_table_paired / create_merged_table_single                    │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Loads MAG stats for MAGs present in p-value tables                  │
│  • Merges on (MAG, contig, position) keys                                   │
│  • Subsets to only significant positions (q < 0.05)                         │
│  → pValue_stats_merged/two_sample_paired_merged_table.tsv.gz                │
│  → pValue_stats_merged/single_sample_{group}_merged_table.tsv.gz            │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│                  STEP 6: ADD QC SUMMARIES (FINAL MEGATABLES)                │
│  create_megatable_paired / create_megatable_single                          │
├─────────────────────────────────────────────────────────────────────────────┤
│  • Flattens QC summaries (group/time metrics become separate columns)       │
│    - Example: breadth_mean_group_fat, breadth_mean_group_fat_time_pre       │
│  • Merges QC data with merged tables on MAG                                 │
│  → megatables/two_sample_paired_megatable.tsv.gz                            │
│  → megatables/single_sample_{group}_megatable.tsv.gz                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                    ↓
                         ┌────────────────────┐
                         │  FINAL MEGATABLES  │
                         │  (Comprehensive)   │
                         └────────────────────┘
```

## Output Structure

### Final Megatables
The megatables contain all information needed for downstream analysis:

#### Columns Included:
1. **Position Identifiers**: MAG, contig, position, gene
2. **P-values & Q-values**: 
   - `min_p_value_CMH`, `q_value_CMH`
   - `min_p_value_two_sample_paired_tTest`, `q_value_two_sample_paired_tTest`
   - `min_p_value_LMM`, `q_value_LMM`
   - `min_p_value_two_sample_unpaired_tTest`, `q_value_two_sample_unpaired_tTest`
   - Additional test-specific columns (LMM_across_time, CMH_across_time, single_sample_tTest)

3. **Coverage Statistics** (per position):
   - `total_coverage`, `mean_coverage`, `std_coverage`
   - `mean_coverage_with_zeros`, `std_coverage_with_zeros`
   - `n_present`, `n_samples`
   - Group-stratified: `mean_coverage_group_{group}`, `std_coverage_group_{group}`
   - Time-stratified: `mean_coverage_group_{group}_time_{time}`

4. **Allele Frequencies** (A, C, G, T):
   - Overall: `mean_freq_A`, `std_freq_A`, etc.
   - Group-stratified: `mean_freq_A_group_{group}`, etc.

5. **QC Metrics** (per MAG):
   - Overall: `num_samples`, `breadth_mean`, `average_coverage_mean`, `median_coverage_mean`
   - Group-stratified: `breadth_mean_group_{group}`, `average_coverage_mean_group_{group}`
   - Time-stratified: `breadth_mean_group_{group}_time_{time}` (paired data only)

## Configuration

### Required Config Parameters (`config.yml`):

```yaml
# Output directory
output_dir: "/path/to/output"

# Input directories
significance_dir_two_sample_paired: "/path/to/significance/two_sample_paired"
significance_dir_single_sample: "/path/to/significance/single_sample"
p_value_summary_dir: "/path/to/p_value_summaries"
metadata_dir: "/path/to/mag_metadata"
qc_dir: "/path/to/qc_results"

# Reference files
fasta: "/path/to/reference.fasta"
mag_mapping: "/path/to/mag_to_contig_mapping.tsv"

# Analysis parameters
single_sample_groups: ["fat", "control"]  # List of groups for single-sample analysis
breadth_threshold: 0.5
cpus: 16

# Optional
mag_list: "/path/to/mag_list.txt"  # Auto-generated if not provided
min_samples: 1
```

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
output_dir/
├── tested_positions/
│   ├── two_sample_paired_mag_positions.tsv.gz
│   └── single_sample_{group}_mag_positions.tsv.gz
├── bh_pValues/
│   ├── two_sample_paired_bh_pValues.tsv.gz
│   └── single_sample_{group}_bh_pValues.tsv.gz
├── qc_results_paired/
│   ├── ALL_MAGs_QC_overall_summary.tsv
│   ├── ALL_MAGs_QC_group_summary.tsv
│   └── ALL_MAGs_QC_group_time_summary.tsv
├── qc_results_single_{group}/
│   ├── ALL_MAGs_QC_overall_summary.tsv
│   └── ALL_MAGs_QC_group_summary.tsv
├── coverage_stats/
│   └── {mag}_stats.tsv (one per MAG)
├── pValue_stats_merged/
│   ├── two_sample_paired_merged_table.tsv.gz
│   └── single_sample_{group}_merged_table.tsv.gz
└── megatables/
    ├── two_sample_paired_megatable.tsv.gz
    └── single_sample_{group}_megatable.tsv.gz
```

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

## Dependencies

- Python 3.8+
- pandas
- snakemake
- tqdm
- AlleleFlux package (for `alleleflux-coverage-allele-stats`, `alleleflux-list-mags`)

