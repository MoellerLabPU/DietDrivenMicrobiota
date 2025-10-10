# Megatable Column Documentation

This document defines every column that appears in the final megatable(s) produced by `megaTable/create_megatable.smk`, based on the exact behavior of the workflow and its referenced scripts in the AlleleFlux repository.

## Table of Contents
1. [Overview](#overview)
2. [Position Identifier Columns](#position-identifier-columns)
3. [Statistical Test Results Columns](#statistical-test-results-columns)
4. [Coverage Statistics Columns](#coverage-statistics-columns)
5. [Allele Frequency Columns](#allele-frequency-columns)
6. [Quality Control Metrics Columns](#quality-control-metrics-columns)
7. [Quick Reference Tables](#quick-reference-tables)
8. [Interpretation Guidelines](#interpretation-guidelines)
9. [Data Quality Notes](#data-quality-notes)

---

## Overview

The megatable is the comprehensive output of the AlleleFlux analysis pipeline, consolidating statistical test results, coverage metrics, allele frequencies, and quality control data into a single table for downstream analysis.

### Outputs

- `megatables/two_sample_paired_megatable.tsv.gz`
- `megatables/single_sample_{group}_megatable.tsv.gz`

Each megatable is built by merging, per-position within a MAG:
- BH-corrected test results (pivoted to wide) from `alleleflux.scripts.preprocessing.p_value_summary`
- Per-position coverage and allele-frequency statistics from `alleleflux.scripts.accessory.coverage_and_allele_stats`
- Per-MAG quality-control summaries from `alleleflux.scripts.preprocessing.quality_control`


### Two Analysis Types

**1. Paired Analysis (`two_sample_paired_megatable.tsv.gz`)**
- Compares pre/end timepoints between two groups
- Includes time-stratified metrics
- Tests: CMH, paired t-test, LMM, unpaired t-test

**2. Single-Sample Analysis (`single_sample_{group}_megatable.tsv.gz`)**
- Analyzes temporal changes within a single group
- One file per group (e.g., "fat", "control")
- Tests: single-sample t-test, LMM across time, CMH across time

### Row Selection Criteria

**Critical:** Only positions meeting statistical significance thresholds are included:
- **Paired analysis:** Positions with `q_value < 0.05` for `two_sample_paired_tTest`
- **Single-sample analysis:** Positions with `q_value < 0.05` for `single_sample_tTest` within the specified group

This filtering ensures the megatable contains only positions with evidence of significant change.

---

## Position Identifier Columns

These columns uniquely identify each genomic position in the analysis.

| Column | Data Type | Source | Description |
|--------|-----------|--------|-------------|
| `MAG` | string | Extracted from filename | Metagenome-Assembled Genome identifier (e.g., "SLG779_DASTool_bins_SLG779_bin.32") |
| `contig` | string | Significance test files | Contig name within the MAG |
| `position` | integer | Significance test files | 1-based position on the contig |
| `gene` | string | Significance test files | Gene identifier if position is within a gene, otherwise NA |

**Uniqueness:** The combination of (MAG, contig, position) is guaranteed to be unique within each megatable.

---

## Statistical Test Results Columns

P-values and Q-values (FDR-corrected p-values) from multiple statistical tests are included for each significant position.

### Column Naming Pattern

```
{statistic}_{test_label}
```

Where:
- **statistic**: `min_p_value` or `q_value`
- **test_label**: Test name, optionally with group suffix

### Paired Analysis Tests

| Test Label | Full Column Names | Description |
|------------|-------------------|-------------|
| `two_sample_paired_tTest` | `min_p_value_two_sample_paired_tTest`<br>`q_value_two_sample_paired_tTest` | Paired t-test comparing pre vs. end within each group |
| `CMH` | `min_p_value_CMH`<br>`q_value_CMH` | Cochran-Mantel-Haenszel test for allele frequency changes |
| `LMM` | `min_p_value_LMM`<br>`q_value_LMM` | Linear Mixed Model for longitudinal changes |
| `two_sample_unpaired_tTest` | `min_p_value_two_sample_unpaired_tTest`<br>`q_value_two_sample_unpaired_tTest` | Unpaired t-test between groups |

**Note:** All positions in the paired megatable have `q_value_two_sample_paired_tTest < 0.05` by design.

### Single-Sample Analysis Tests

For group "fat" as an example:

| Test Label | Full Column Names | Description |
|------------|-------------------|-------------|
| `single_sample_tTest_fat` | `min_p_value_single_sample_tTest_fat`<br>`q_value_single_sample_tTest_fat` | Single-sample t-test for temporal change within "fat" group |
| `LMM_across_time` | `min_p_value_LMM_across_time`<br>`q_value_LMM_across_time` | Linear Mixed Model across timepoints |
| `CMH_across_time` | `min_p_value_CMH_across_time`<br>`q_value_CMH_across_time` | Cochran-Mantel-Haenszel test across timepoints |

**Note:** The test label includes the group name (e.g., `_fat`, `_control`) for single-sample t-test results.

### Understanding P-values and Q-values

- **min_p_value**: The minimum p-value from the test across all 4 alleles (most significant result)
- **q_value**: False Discovery Rate-corrected minimum p-value (controls for multiple testing)
- **Interpretation**: `q_value < 0.05` indicates statistical significance after multiple testing correction

### NaN in Statistical Results

- If a given test was not run for a position/context, the corresponding p-value/q-value columns are NaN.
- In single-sample megatables, the pivot may include columns for multiple groups; for positions selected by `{group}`, other groups' columns are typically NaN unless the same position also appeared in their tests.

---

## Coverage Statistics Columns

**IMPORTANT: These are PER-POSITION statistics** - each row represents coverage at a single genomic position.

Coverage metrics quantify sequencing depth at each position, with stratification by group and time.

### ⚠️ CRITICAL: Sample Filtering in Coverage Statistics

**Coverage statistics are calculated using only QC-passing samples.** The workflow applies breadth threshold filtering through the following mechanism. This is different than the QC metrics described below:

1. **QC Step**: `quality_control.py`generates QC files (`*_QC.tsv`) with a `breadth_threshold_passed` column indicating which samples passed the breadth threshold for each MAG.

2. **Coverage Stats Step**: `create_megatable.smk` passes `--qc_dir` to `alleleflux-coverage-allele-stats`, which:
   - Reads the QC files to identify samples that passed breadth threshold
   - Filters the metadata to include only QC-passing samples
   - Computes coverage statistics using this filtered sample set

**Result:** `n_samples` reflects the number of QC-passing samples for that MAG. All coverage and allele frequency statistics are calculated using only samples that passed the breadth threshold.

**Code Citations:**
```python
# coverage_and_allele_stats.py lines 717-729
# Apply QC filtering if requested
qc_filtered_samples = None
if qc_dir:
    qc_filtered_samples = read_qc_filtered_samples(qc_dir, mag_id)

# Filter metadata to only include QC-passed samples
if qc_filtered_samples is not None:
    meta = meta[meta_sample_ids.isin(qc_filtered_samples)].copy()
```

### Dual Statistics Approach

Coverage statistics are computed in **TWO ways** to handle samples with zero coverage:

#### **Present-Only Statistics (EXCLUDES ZEROS)**
**Columns:** `mean_coverage`, `std_coverage`

**Calculation:**
```
mean_coverage = sum(coverage for samples with coverage > 0) / n_present

where n_present = count of QC-passing samples with coverage > 0 at this position
```

- **Denominator:** Only QC-passing samples with non-zero coverage at that position (`n_present`)
- **Zero-coverage samples:** Completely excluded from calculations
- **Use case:** Understanding actual sequencing depth when the MAG is detected at this position
- **Biological meaning:** "When this position is covered, what depth do we typically see?"

**Example:**
```
Among QC-passing samples at this position:
Sample A: 10x coverage
Sample B: 15x coverage
Sample C: 0x coverage (MAG not detected - EXCLUDED)
Sample D: 0x coverage (MAG not detected - EXCLUDED)

mean_coverage = (10 + 15) / 2 = 12.5x
n_present = 2
```

#### **Include-Zeros Statistics (INCLUDES ZEROS)**
**Columns:** `mean_coverage_with_zeros`, `std_coverage_with_zeros`

**Calculation:**
```
mean_coverage_with_zeros = sum(coverage across ALL QC-passing samples) / n_samples

where n_samples = total number of QC-passing samples for this MAG
```

- **Denominator:** All QC-passing samples for the MAG (`n_samples`)
- **Zero-coverage samples:** Included as 0 in calculations
- **Use case:** Population-level prevalence and overall detection patterns among quality samples
- **Biological meaning:** "What's the average coverage across the QC-passing population?"

**Example (same data, assume 4 QC-passing samples total):**
```
Sample A: 10x coverage
Sample B: 15x coverage
Sample C: 0x coverage (MAG not present - INCLUDED as 0)
Sample D: 0x coverage (MAG not present - INCLUDED as 0)

mean_coverage_with_zeros = (10 + 15 + 0 + 0) / 4 = 6.25x
n_samples = 4
```

### Why Both Metrics Exist

1. **Exclude-zeros metrics** (`mean_coverage`): Show actual depth when position IS covered
2. **Include-zeros metrics** (`mean_coverage_with_zeros`): Show overall prevalence and population-level coverage
3. **Together**: Large differences indicate sparse but well-covered positions; similar values indicate ubiquitous coverage

**Interpretation Example:**
```
If mean_coverage = 20x but mean_coverage_with_zeros = 4x:
→ Position is well-covered (20x) when present
→ But only detected in ~20% of QC-passing samples (4/20 = 0.2)
→ This is a sparsely-distributed but high-quality position
```

### Base Coverage Columns

**Per-position coverage statistics for all QC-passing samples:**

| Column | Zero Handling | Denominator | Description |
|--------|---------------|-------------|-------------|
| `total_coverage` | N/A | N/A | Sum of coverage across all QC-passing samples at this position |
| `mean_coverage` | ⚠️ **EXCLUDES** zeros | n_present | Mean coverage for QC-passing samples with coverage > 0 |
| `std_coverage` | ⚠️ **EXCLUDES** zeros | n_present | Standard deviation for present QC-passing samples |
| `mean_coverage_with_zeros` | ✅ **INCLUDES** zeros | n_samples | Mean coverage across ALL QC-passing samples (zeros included) |
| `std_coverage_with_zeros` | ✅ **INCLUDES** zeros | n_samples | Standard deviation across ALL QC-passing samples |
| `n_present` | N/A | N/A | Number of QC-passing samples with non-zero coverage |
| `n_samples` | N/A | N/A | Total number of QC-passing samples for this MAG |

### Group-Stratified Coverage

**Per-position coverage statistics stratified by experimental group:**

**Pattern:** `{metric}_group_{group}`

**Examples:**
- `mean_coverage_group_fat`
- `std_coverage_group_fat`
- `mean_coverage_with_zeros_group_fat`
- `std_coverage_with_zeros_group_fat`

**Zero Handling:** Group-stratified metrics follow the same pattern as base metrics:
- `mean_coverage_group_{group}` = **EXCLUDES** zeros (denominator = n_present in that group)
- `mean_coverage_with_zeros_group_{group}` = **INCLUDES** zeros (denominator = n_samples in that group)

### Time-Stratified Coverage (Paired Analysis Only)

**Pattern:** `{metric}_group_{group}_time_{time}`

**Examples:**
- `mean_coverage_group_fat_time_pre`
- `std_coverage_group_fat_time_end`
- `mean_coverage_with_zeros_group_control_time_pre`

**Zero Handling:** Time-stratified metrics follow the same exclude/include pattern.

---

## Allele Frequency Columns

Allele frequencies represent the proportion of each nucleotide (A, C, G, T) at a position.

### ⚠️ Note: Allele Frequencies Are ALWAYS Present-Only

Unlike coverage metrics, allele frequencies are **only computed for QC-passing samples where the position has coverage > 0**. There are no "include-zeros" variants for allele frequencies, as it would be meaningless to average in zero-coverage samples.

### Base Allele Frequency Columns

**Pattern:** `{statistic}_freq_{nucleotide}`

| Column | Description |
|--------|-------------|
| `mean_freq_A` | Mean frequency of adenine across QC-passing samples with coverage |
| `std_freq_A` | Standard deviation of adenine frequency |
| `mean_freq_C` | Mean frequency of cytosine across QC-passing samples with coverage |
| `std_freq_C` | Standard deviation of cytosine frequency |
| `mean_freq_G` | Mean frequency of guanine across QC-passing samples with coverage |
| `std_freq_G` | Standard deviation of guanine frequency |
| `mean_freq_T` | Mean frequency of thymine across QC-passing samples with coverage |
| `std_freq_T` | Standard deviation of thymine frequency |

### Group-Stratified Allele Frequencies

**Pattern:** `{statistic}_freq_{nucleotide}_group_{group}`

**Examples:**
- `mean_freq_A_group_fat`
- `std_freq_A_group_fat`
- `mean_freq_G_group_control`

### Time-Stratified Allele Frequencies (Paired Analysis Only)

**Pattern:** `{statistic}_freq_{nucleotide}_group_{group}_time_{time}`

**Examples:**
- `mean_freq_A_group_fat_time_pre`
- `std_freq_C_group_control_time_end`

### Understanding Allele Frequencies

```
Example at position X (among QC-passing samples):
- Sample 1: 10x coverage, all reads show 'A' → freq_A=1.0, others=0.0
- Sample 2: 20x coverage, 15 'A' + 5 'G' → freq_A=0.75, freq_G=0.25
- Sample 3: 0x coverage → excluded from calculation

mean_freq_A = (1.0 + 0.75) / 2 = 0.875
mean_freq_G = (0.0 + 0.25) / 2 = 0.125
```

**Note:** At each position in each sample, frequencies sum to 1.0 (or close to 1.0 accounting for rounding).

---

## Quality Control Metrics Columns

QC metrics are computed per-MAG (not per-position) and provide information about sequencing quality and coverage breadth. These metrics are merged into the megatable based on the MAG column, so all positions from the same MAG share the same QC values.

### QC Metrics in Megatable Workflow Context

**The megatable workflow implements dual breadth calculation to enable both position-specific and genome-wide quality assessment.** When a positions file is provided with `positions_denominator="positions"`, the workflow calculates two distinct breadth metrics at the per-sample level:

**Dual Breadth Mode (Per-Sample Calculation):**
1. **`breadth`: Coverage at the **tested positions only**
   - Captured: After position filtering
   - Count: Positions with coverage ≥ 1 in the FILTERED dataframe
   - Denominator: Number of positions in the positions file
   - Used for position-specific quality metrics

2. **`breadth_genome`: Coverage across the **entire genome**
   - Captured: **BEFORE** position filtering
   - Count: Positions with coverage ≥ 1 in the FULL MAG profile
   - Denominator: Full genome size
   - Used for breadth threshold filtering

**Quality Filtering Behavior:**
- Threshold checking uses `breadth_genome`
- Sample-level `breadth_threshold_passed` correctly reflects genome-wide breadth
- Only samples passing `breadth_genome` threshold are included in aggregated QC summaries

### Denominator Considerations

When QC is run with a positions file (as in the megatable workflow):
- **Denominator mode**: Set to `"positions"` (see [`create_megatable.smk:524`](create_megatable.smk:524))
- **`breadth` calculation**: Proportion of tested positions with coverage ≥ 1 (denominator = number of specified positions)
- **`breadth_genome` calculation**: Proportion of ALL genome positions with coverage ≥ 1 (denominator = genome size, calculated BEFORE position filtering)
- **Average coverage**: Sum of coverage at tested positions / number of tested positions
- **All other coverage metrics**: Calculated from FILTERED positions only, not the entire genome

This means most QC metrics describe quality at the specified positions only, while `breadth_genome` provides genome-wide context.

### Overall QC Metrics

These metrics summarize quality across all samples for each MAG.

**Source:** `ALL_MAGs_QC_overall_summary.tsv`

| Column | Description | Calculation Details (the mean is across the samples considered) |
|--------|-------------|-----------------------------------------------------------------|
| `num_samples` | Total number of samples passing breadth threshold | Samples with breadth_genome ≥ threshold |
| `breadth_mean` | Mean breadth of coverage at tested positions | Mean proportion of tested positions with coverage ≥ 1 (calculated AFTER position filtering) |
| `average_coverage_mean` | Mean of average coverage per sample | Mean of (sum coverage at tested positions / number of tested positions) per sample |
| `median_coverage_mean` | Mean of median coverage per sample | Mean of median coverage at tested positions (zeros excluded per sample) |
| `median_coverage_including_zeros_mean` | Median including zeros for absent positions | Mean of median across tested positions with zeros for absent positions |
| `coverage_std_mean` | Mean of coverage standard deviation | Mean of per-sample std at tested positions (zeros excluded) |
| `coverage_std_including_zeros_mean` | Std including zeros | Mean of per-sample std including zeros for absent positions |

**Note:** `length_weighted_coverage` is omitted when a positions file is used, as it is only meaningful for full genome analysis.

### Group-Stratified QC Metrics

**Source:** `ALL_MAGs_QC_group_summary.tsv` (pivoted)

**Pattern:** `{metric}_group_{group}`

**Examples:**
- `breadth_mean_group_fat` (position-specific breadth)
- `average_coverage_mean_group_control`
- `median_coverage_mean_group_control`

**Note:** Group-stratified `breadth_genome` columns will NOT appear due to the same aggregation bug.

### Time-Stratified QC Metrics (Paired Analysis Only)

**Source:** `ALL_MAGs_QC_group_time_summary.tsv` (pivoted)

**Pattern:** `{metric}_group_{group}_time_{time}`

**Examples:**
- `breadth_mean_group_fat_time_pre` (position-specific breadth)
- `average_coverage_mean_group_fat_time_end`
- `median_coverage_mean_group_control_time_pre`


### QC Metric Interpretation

```
Example MAG at tested positions:
- breadth_mean = 0.85 → 85% of tested positions are covered on average
- breadth_mean_group_fat = 0.90 → 90% of tested positions covered in fat group
- average_coverage_mean = 45.2x → Average depth of 45x across tested positions
```

---

## Quick Reference Tables

### Column Families Summary

| Family | Columns/Patterns | Denominator | Zero Handling | Notes |
|--------|------------------|-------------|---------------|-------|
| Identity | `MAG`, `contig`, `position`, `gene` | n/a | n/a | Join keys + optional gene context |
| P-values | `min_p_value_*`, `q_value_*` | n/a | n/a | Labels reflect test + optional `_across_time` + optional `_GROUP` |
| Coverage overall | `total_coverage`, `n_present`, `n_samples` | see below | see below | Base counts/sums (QC-passing samples) |
| Coverage present-only | `mean_coverage`, `std_coverage` | `n_present` | exclude zeros | QC-passing samples with coverage > 0 only |
| Coverage include-zeros | `mean_coverage_with_zeros`, `std_coverage_with_zeros` | `n_samples` | include zeros | All QC-passing samples (zeros as 0) |
| Coverage grouped | `…_group_<G>`, `…_group_<G>_time_<T>` | group-wise | per variant | Mirrors overall logic within slices |
| Allele freq overall | `mean_freq_[ACGT]`, `std_freq_[ACGT]` | `n_present` | exclude zeros | Frequencies computed only where covered |
| Allele freq grouped | `mean_freq_[ACGT]_group_*…` | group-wise | exclude zeros | Same as above within slices |
| QC overall | `num_samples`, `breadth_mean`, etc. | samples passing breadth_genome threshold | metric-specific | Per-MAG means|
| QC grouped | `<field>_group_<G>[ _time_<T>]` | samples passing breadth_genome threshold | metric-specific | Flattened pivots by group/time; no `breadth_genome` variants |

### Zero-Handling Summary Table

| Metric Category | EXCLUDES Zeros (Present-Only) | INCLUDES Zeros |
|----------------|-------------------------------|----------------|
| **Coverage - Base** | `mean_coverage`<br>`std_coverage` | `mean_coverage_with_zeros`<br>`std_coverage_with_zeros` |
| **Coverage - Group** | `mean_coverage_group_{group}`<br>`std_coverage_group_{group}` | `mean_coverage_with_zeros_group_{group}`<br>`std_coverage_with_zeros_group_{group}` |
| **Coverage - Time** | `mean_coverage_group_{group}_time_{time}`<br>`std_coverage_group_{group}_time_{time}` | `mean_coverage_with_zeros_group_{group}_time_{time}`<br>`std_coverage_with_zeros_group_{group}_time_{time}` |
| **Allele Frequencies** | All allele freq columns<br>*(always present-only)* | *(not applicable)* |
| **QC Metrics** | `breadth_mean` (position-specific, AFTER filtering) | Not applicable<br>*(Breadth metrics distinguish timing: post-filter vs pre-filter)* |

### Test Result Column Finder

**Paired Analysis:**
```
CMH test:           min_p_value_CMH, q_value_CMH
Paired t-test:      min_p_value_two_sample_paired_tTest, q_value_two_sample_paired_tTest
LMM:                min_p_value_LMM, q_value_LMM
Unpaired t-test:    min_p_value_two_sample_unpaired_tTest, q_value_two_sample_unpaired_tTest
```

**Single-Sample Analysis (example group "fat"):**
```
T-test:             min_p_value_single_sample_tTest_fat, q_value_single_sample_tTest_fat
LMM across time:    min_p_value_LMM_across_time, q_value_LMM_across_time
CMH across time:    min_p_value_CMH_across_time, q_value_CMH_across_time
```

---

## Interpretation Guidelines

### Reading Zero-Handling Patterns

**Scenario 1: Sparse MAG with good coverage when present**
```
n_present = 20 (out of 100 QC-passing samples)
mean_coverage = 20.0x                    ← EXCLUDES zeros (only 20 present samples)
mean_coverage_with_zeros = 4.0x          ← INCLUDES zeros (all 100 QC-passing samples)
```
**Interpretation:** The MAG appears in only 20% of QC-passing samples, but when it does appear, it's well-covered at 20x depth.

**Scenario 2: Ubiquitous MAG**
```
n_present = 98 (out of 100 QC-passing samples)
mean_coverage = 15.1x                    ← EXCLUDES zeros
mean_coverage_with_zeros = 14.8x         ← INCLUDES zeros (nearly all samples)
```
**Interpretation:** The MAG is nearly universal among QC-passing samples. Both metrics converge because almost all samples have coverage.

**Scenario 3: Group-specific MAG**
```
mean_coverage_group_fat = 0.5x
mean_coverage_group_control = 15.0x
```
**Interpretation:** This MAG is predominantly found in the control group, not the fat group.

---

## Data Quality Notes

### QC Filtering Effects

The workflow applies quality control at multiple stages:

1. **Position filtering:** Only positions passing significance thresholds (q < 0.05) are included in the megatable
2. **Sample filtering in coverage stats:** Samples not meeting breadth threshold are excluded from coverage/allele calculations via `--qc_dir` parameter
3. **Breadth threshold filtering in QC metrics:** Samples are filtered based on genome-wide breadth (`breadth_genome`) calculated BEFORE position filtering

### Understanding Coverage Calculation Timing

**Critical Distinction - Two Types of Metrics:**

1. **Per-Position Coverage/Allele Statistics** (calculated AFTER position filtering):
   - Source: `coverage_and_allele_stats.py`
   - Calculates per-position statistics from specified positions specified from the positions file
   - All coverage metrics (`mean_coverage`, `std_coverage`, etc.) are calculated from the FILTERED position set
   - Uses `--qc_dir` parameter to identify QC-passing samples
   - Filters to samples passing breadth threshold (via genome-wide `breadth_genome`)
   - `n_samples` = number of QC-passing samples (those with breadth_genome ≥ threshold)
   - These appear as per-position rows in the megatable

2. **Per-MAG QC Summary Metrics** (with dual breadth calculation):
   - Source: `quality_control.py`
   - Uses `--positions_file` with `--positions_denominator positions`
   - Implements **dual breadth calculation** at per-sample level:
      - **`breadth_genome`**: Count captured BEFORE filtering → used for threshold checking
      - **`breadth`**: Count captured AFTER filtering → used for position-specific quality
      - **All other metrics**: Calculated AFTER filtering → all use filtered dataframe
   - Threshold filtering uses genome-wide `breadth_genome`
   - `num_samples` = number of samples passing breadth_genome threshold
   - These appear as MAG-level columns in the megatable (all positions from same MAG share same values)

### Missing Data Handling

**NaN values may appear in columns when:**

1. **Insufficient data for statistics:**
   - `std_coverage` = NaN when only 1 sample has coverage
   - Group-specific metrics = NaN when group has no samples with coverage

2. **Test not applicable:**
   - Time-stratified columns = NaN in single-sample analysis when that timepoint doesn't exist
   - Group-specific p-values = NaN when position wasn't tested for that group

3. **MAG not in analysis:**
   - QC metric columns = NaN when MAG wasn't included in quality control

### Column Availability by Analysis Type

**Paired Analysis ONLY:**
- Time-stratified coverage: `mean_coverage_group_{group}_time_{time}`
- Time-stratified allele freq: `mean_freq_A_group_{group}_time_{time}`
- Time-stratified QC: `breadth_mean_group_{group}_time_{time}`
- Paired t-test p-values: `min_p_value_two_sample_paired_tTest`

**Single-Sample Analysis ONLY:**
- Single-sample t-test p-values: `min_p_value_single_sample_tTest_{group}`
- Group-specific time columns (if multiple timepoints exist within group)


---

## Additional Notes

### Column Order

Columns appear in this general order:
1. Position identifiers (MAG, contig, position, gene)
2. Statistical test results (p-values, q-values)
3. Coverage statistics (base, then group, then time)
4. Allele frequencies (base, then group, then time)
5. QC metrics (overall, then group, then time)

The exact order may vary slightly based on the specific analysis configuration.

### Dynamic Column Generation

Some columns are dynamically generated based on:
- **Groups:** Defined in `config.yml` under `single_sample_groups`
- **Timepoints:** Extracted from metadata files
- **Tests performed:** Based on available p-value summary files

### Naming Conventions

- Group/time labels are sanitized to `[A-Za-z0-9_.-]`; other characters replaced by `_`
- Column suffixes follow pattern: `_group_{group}` or `_group_{group}_time_{time}`
- All numeric columns are rounded to 6 decimal places

---

## Provenance (for maintainers)

- Workflow driving merges and naming: `megaTable/create_megatable.smk`
- P-value summarization and BH: `alleleflux.scripts.preprocessing.p_value_summary`
- Coverage and allele stats: `alleleflux.scripts.accessory.coverage_and_allele_stats`
- QC sample metrics and aggregates: `alleleflux.scripts.preprocessing.quality_control`

