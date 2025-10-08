import os
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import logging

# Load configuration
configfile: "config.yml"

# Set up logging
workflow_logger = logging.getLogger("__name__")   # any non-root name
workflow_logger.propagate = False                   # don't bubble to root
workflow_logger.setLevel(logging.INFO)
_h = logging.StreamHandler()
_h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
workflow_logger.addHandler(_h)


def _process_and_pivot_pvalue_files(data_directory, p_value_files, merge_keys, final_df, q_threshold=0.05, group=None):
    """Process a list of p-value summary files and pivot them to a wide format.

    Parameters
    - data_directory: base dir (used for logging only)
    - p_value_files: list of file paths
    - merge_keys: list of columns to join on
    - final_df: DataFrame with base positions to keep (will be left-joined)
    - q_threshold: numeric threshold (not used for selection here but kept for API)
    - group: optional group name to filter rows in files that contain a 'group' column

    Returns a DataFrame (final_df merged with pivoted p-value columns).
    """
    processed_data_chunks = []
    for p_value_file_path in p_value_files:
        if not os.path.exists(p_value_file_path):
            workflow_logger.warning(f"File not found, skipping: {p_value_file_path}")
            continue

        p_value_df = pd.read_csv(p_value_file_path, sep="\t")
        workflow_logger.info(f"Processing '{p_value_file_path}' with shape {p_value_df.shape}")

        # If a group filter was requested, apply it when possible
        if group is not None and 'group' in p_value_df.columns:
            p_value_df = p_value_df[p_value_df['group'] == group]

        # Validate that the file has the necessary 'test_type' column
        if 'test_type' not in p_value_df.columns or p_value_df['test_type'].dropna().nunique() == 0:
            workflow_logger.warning(f"SKIP '{p_value_file_path}': missing/empty 'test_type'")
            continue

        # Normalize the data by ensuring a 'group_analyzed' column exists
        if 'group_analyzed' not in p_value_df.columns:
            p_value_df['group_analyzed'] = pd.NA

        # Reduce the dataframe to only the columns needed for the pivot
        required_columns = merge_keys + ['test_type', 'group_analyzed', 'min_p_value', 'q_value']
        p_value_df = p_value_df[[c for c in required_columns if c in p_value_df.columns]].copy()

        # Build label from test_type and optionally filename-based special cases
        filename = os.path.basename(p_value_file_path).lower()
        p_value_df['label'] = p_value_df['test_type'].astype(str).str.strip()

        # For files that indicate "across_time" in the filename adjust labels
        is_lmm_across_time = 'lmm_across_time' in filename
        is_cmh_across_time = 'cmh_across_time' in filename
        lmm_rows_mask = (p_value_df['label'].str.lower() == 'lmm') & is_lmm_across_time
        cmh_rows_mask = (p_value_df['label'].str.lower() == 'cmh') & is_cmh_across_time
        p_value_df.loc[lmm_rows_mask, 'label'] = 'LMM_across_time'
        p_value_df.loc[cmh_rows_mask, 'label'] = 'CMH_across_time'

        # Append group suffix when group_analyzed present
        has_valid_group_mask = p_value_df['group_analyzed'].notna() & (p_value_df['group_analyzed'].astype(str).str.strip() != "")
        p_value_df.loc[has_valid_group_mask, 'label'] = (
            p_value_df.loc[has_valid_group_mask, 'label'] + '_' + p_value_df.loc[has_valid_group_mask, 'group_analyzed'].astype(str)
        )

        processed_data_chunks.append(p_value_df[merge_keys + ['label', 'min_p_value', 'q_value']])

    if not processed_data_chunks:
        # Nothing to merge; return final_df unchanged
        return final_df

    combined_tidy_df = pd.concat(processed_data_chunks, ignore_index=True)

    dup_mask = combined_tidy_df.duplicated(subset=merge_keys + ['label'], keep=False)
    if dup_mask.any():
        raise ValueError(f"Duplicated rows across p-value files: {combined_tidy_df.loc[dup_mask]}")

    pivoted_wide_df = (
        combined_tidy_df.set_index(merge_keys + ['label'])[['min_p_value', 'q_value']]
            .unstack('label')
    )

    # Flatten MultiIndex columns
    pivoted_wide_df.columns = [f"{value}_{key}" for value, key in pivoted_wide_df.columns]
    pivoted_wide_df = pivoted_wide_df.reset_index()

    merged = final_df.merge(pivoted_wide_df, on=merge_keys, how='left')
    return merged


def _load_and_flatten_qc_summaries(qc_dir):
    """Load QC summary files and flatten group/time columns.
    
    Parameters
    - qc_dir: Directory containing QC summary files
    
    Returns a flattened DataFrame with MAG_ID as the key.
    """
    overall_file = os.path.join(qc_dir, "ALL_MAGs_QC_overall_summary.tsv")
    group_file = os.path.join(qc_dir, "ALL_MAGs_QC_group_summary.tsv")
    group_time_file = os.path.join(qc_dir, "ALL_MAGs_QC_group_time_summary.tsv")
    
    # Load overall summary
    if not os.path.exists(overall_file):
        raise FileNotFoundError(f"Overall QC summary not found: {overall_file}")
    
    df_overall = pd.read_csv(overall_file, sep='\t')
    workflow_logger.info(f"Loaded overall QC summary with {df_overall.shape[0]:,} MAGs")
    
    # Start with overall summary
    result = df_overall.copy()
    
    # Load and pivot group summary if it exists
    if os.path.exists(group_file):
        df_group = pd.read_csv(group_file, sep='\t')
        workflow_logger.info(f"Loaded group QC summary with {df_group.shape[0]:,} rows")
        
        # Get columns to pivot (exclude MAG_ID and group)
        value_cols = [col for col in df_group.columns if col not in ['MAG_ID', 'group']]
        
        # Pivot: one row per MAG_ID, columns for each group
        df_group_pivot = df_group.pivot(index='MAG_ID', columns='group', values=value_cols)
        
        # Flatten multi-index columns: (metric, group) -> metric_group_<group>
        df_group_pivot.columns = [f"{metric}_group_{group}" for metric, group in df_group_pivot.columns]
        df_group_pivot = df_group_pivot.reset_index()
        
        # Merge with result
        result = result.merge(df_group_pivot, on='MAG_ID', how='left')
        workflow_logger.info(f"Merged group QC data, now {result.shape[1]:,} columns")
    
    # Load and pivot group_time summary if it exists (only for paired data)
    if os.path.exists(group_time_file):
        df_group_time = pd.read_csv(group_time_file, sep='\t')
        workflow_logger.info(f"Loaded group_time QC summary with {df_group_time.shape[0]:,} rows")
        
        # Get columns to pivot (exclude MAG_ID, group, and time)
        value_cols = [col for col in df_group_time.columns if col not in ['MAG_ID', 'group', 'time']]
        
        # Pivot: one row per MAG_ID, columns for each (group, time) combination
        df_group_time_pivot = df_group_time.pivot(
            index='MAG_ID', 
            columns=['group', 'time'], 
            values=value_cols
        )
        
        # Flatten multi-index columns: (metric, group, time) -> metric_group_<group>_time_<time>
        df_group_time_pivot.columns = [
            f"{metric}_group_{group}_time_{time}" 
            for metric, group, time in df_group_time_pivot.columns
        ]
        df_group_time_pivot = df_group_time_pivot.reset_index()
        
        # Merge with result
        result = result.merge(df_group_time_pivot, on='MAG_ID', how='left')
        workflow_logger.info(f"Merged group_time QC data, now {result.shape[1]:,} columns")
    
    return result


def _merge_coverage_with_pvalues(bh_pvalues_file, coverage_stats_dir, output_file, group_label=None):
    """Merge coverage statistics with BH-corrected p-values.
    
    Parameters
    - bh_pvalues_file: Path to BH p-values file (gzipped TSV)
    - coverage_stats_dir: Directory containing MAG coverage stats files
    - output_file: Path for output merged_table file
    - group_label: Optional group label for logging (e.g., 'fat', 'control')
    
    Returns the merged DataFrame.
    """
    group_info = f" for group {group_label}" if group_label else ""
    
    # Load the BH p-value table
    df_pvalues = pd.read_csv(bh_pvalues_file, sep='\t', compression='gzip')
    workflow_logger.info(f"Loaded BH p-values{group_info} with {df_pvalues.shape[0]:,} positions")
    
    # Get unique MAGs from the p-value table
    mags_needed = df_pvalues['MAG'].unique()
    workflow_logger.info(f"Found {len(mags_needed):,} unique MAGs in p-value table")
    
    # Load and concatenate coverage stats for relevant MAGs
    coverage_dfs = []
    for mag_id in tqdm(mags_needed, desc="Loading MAG coverage files", unit="MAG"):
        coverage_file = os.path.join(coverage_stats_dir, f"{mag_id}_mean_coverage.tsv")
        if os.path.exists(coverage_file):
            df_cov = pd.read_csv(coverage_file, sep='\t')
            df_cov['MAG'] = mag_id
            coverage_dfs.append(df_cov)
        else:
            workflow_logger.warning(f"Coverage file not found for MAG {mag_id}")
    
    if not coverage_dfs:
        raise ValueError("No coverage stats files found for any MAGs in p-value table")
    
    df_coverage = pd.concat(coverage_dfs, ignore_index=True)
    workflow_logger.info(f"Loaded coverage stats with {df_coverage.shape[0]:,} total positions across all MAGs")
    
    # Merge on MAG, contig, and position
    merged_table = df_pvalues.merge(
        df_coverage,
        on=['MAG', 'contig', 'position'],
        how='left'
    )
    
    workflow_logger.info(f"Created merged_table{group_info} with {merged_table.shape[0]:,} positions and {merged_table.shape[1]:,} columns")
    
    # Write output
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    merged_table.to_csv(output_file, sep='\t', index=False, compression='gzip')
    workflow_logger.info(f"Wrote merged_table to {output_file}")
    
    return merged_table


# Base output directory for this workflow
BASE_OUTPUT_DIR = config["output_dir"]

# --- Target Files ---
# Paired sample outputs
PAIRED_POSITIONS_FILE = os.path.join(BASE_OUTPUT_DIR, "tested_positions", "two_sample_paired_mag_positions.tsv.gz")
PAIRED_BH_pValues = os.path.join(BASE_OUTPUT_DIR, "bh_pValues", "two_sample_paired_bh_pValues.tsv.gz")
PAIRED_QC_SUMMARY = os.path.join(BASE_OUTPUT_DIR, "qc_results_paired", "ALL_MAGs_QC_overall_summary.tsv")

# Single sample outputs, expanded for each group defined in the config
SINGLE_SAMPLE_GROUPS = config.get("single_sample_groups", [])
SINGLE_POSITIONS_FILES = expand(os.path.join(BASE_OUTPUT_DIR, "tested_positions", "single_sample_{group}_mag_positions.tsv.gz"), group=SINGLE_SAMPLE_GROUPS)
SINGLE_BH_pValues = expand(os.path.join(BASE_OUTPUT_DIR, "bh_pValues", "single_sample_{group}_bh_pValues.tsv.gz"), group=SINGLE_SAMPLE_GROUPS)
SINGLE_QC_SUMMARY = expand(os.path.join(BASE_OUTPUT_DIR, "qc_results_single_{group}", "ALL_MAGs_QC_overall_summary.tsv"), group=SINGLE_SAMPLE_GROUPS)

# Coverage and allele stats outputs
# Determine MAG list path (either from config or will be generated)
MAG_LIST_FILE = config.get("mag_list", os.path.join(BASE_OUTPUT_DIR, "mag_list.txt"))

# If mag_list is not provided in config, generate it now
if not config.get("mag_list"):
    import subprocess
    mag_list_output = MAG_LIST_FILE
    os.makedirs(os.path.dirname(mag_list_output), exist_ok=True)
    min_samples = config.get("min_samples", 1)
    
    subprocess.run([
        "alleleflux-list-mags",
        "--metadata_dir", config["metadata_dir"],
        "--output", mag_list_output,
    ], check=True)
    
    workflow_logger.info(f"Generated MAG list at {mag_list_output}")

# Read the MAG list (either provided or just generated)
with open(MAG_LIST_FILE) as f:
    MAGS = [line.strip() for line in f if line.strip()]

COVERAGE_STATS_FILES = expand(os.path.join(BASE_OUTPUT_DIR, "coverage_stats", "{mag}_mean_coverage.tsv"), mag=MAGS)


PAIRED_MERGED_TABLE_1 = os.path.join(BASE_OUTPUT_DIR, "pValue_stats_merged", "two_sample_paired_merged_table.tsv.gz")
SINGLE_MERGED_TABLES_1 = expand(os.path.join(BASE_OUTPUT_DIR, "pValue_stats_merged", "single_sample_{group}_merged_table.tsv.gz"), group=SINGLE_SAMPLE_GROUPS)

# Final megatables with QC data merged
PAIRED_MEGATABLE = os.path.join(BASE_OUTPUT_DIR, "megatables", "two_sample_paired_megatable.tsv.gz")
SINGLE_MEGATABLES = expand(os.path.join(BASE_OUTPUT_DIR, "megatables", "single_sample_{group}_megatable.tsv.gz"), group=SINGLE_SAMPLE_GROUPS)

rule all:
    """
    Specifies the final files the workflow should produce.
    """
    input:
        PAIRED_POSITIONS_FILE,
        SINGLE_POSITIONS_FILES,
        PAIRED_BH_pValues,
        SINGLE_BH_pValues,
        PAIRED_QC_SUMMARY,
        SINGLE_QC_SUMMARY,
        COVERAGE_STATS_FILES,
        PAIRED_MERGED_TABLE_1,
        SINGLE_MERGED_TABLES_1,
        PAIRED_MEGATABLE,
        SINGLE_MEGATABLES,

# --- Step 1: Extract MAG Positions ---

rule extract_mag_positions:
    """
    Runs extract_mag_positions.py for both 'two_sample_paired' and 'single_sample'
    test types. This rule uses a single script execution to generate all required
    position files for both test types, which is more efficient than separate rules.
    """
    input:
        significance_dir_two_sample_paired=config["significance_dir_two_sample_paired"],
        significance_dir_single_sample=config["significance_dir_single_sample"],
    output:
        # The script will create these files, so we list them as outputs
        paired=PAIRED_POSITIONS_FILE,
        single=SINGLE_POSITIONS_FILES,
    params:
        outdir=os.path.join(BASE_OUTPUT_DIR, "tested_positions"),
        script="../extract_mag_positions.py",
    threads: 1
    shell:
        """
        echo "Running for two_sample_paired..."
        python {params.script} \
            --input-dir {input.significance_dir_two_sample_paired} \
            --test-type two_sample_paired \
            --out {params.outdir}

        echo "Running for single_sample..."
        python {params.script} \
            --input-dir {input.significance_dir_single_sample} \
            --test-type single_sample \
            --out {params.outdir}
        """

# --- Step 1.5: Filter p_summary for significant positions ---

rule bh_pValues_paired:
    """
    Filter p_summary for two_sample_paired_tTest with q_value < 0.05,
    then merge p-values from CMH, LMM, and unpaired t-test for those positions.
    """
    input:
        p_value_dir=config["p_value_summary_dir"],
    output:
        filtered=PAIRED_BH_pValues,
    params:
        q_threshold=0.05,
    threads: 1
    run:
        data_directory = input.p_value_dir
        
        # Define p-value files for paired analysis
        p_value_files = [
            os.path.join(data_directory, 'p_value_summary_cmh_pre_end.tsv'),
            os.path.join(data_directory, 'p_value_summary_two_sample_paired_pre_end.tsv'),
            os.path.join(data_directory, 'p_value_summary_lmm_pre_end.tsv'),
            os.path.join(data_directory, 'p_value_summary_two_sample_unpaired_pre_end.tsv')
        ]
        
        # Define the core columns used to join the dataframes
        merge_keys = ['mag_id', 'contig', 'position', 'gene_id']
        
        # --- Step 1: Filter for two_sample_paired_tTest with q_value < 0.05 ---
        paired_file = os.path.join(data_directory, 'p_value_summary_two_sample_paired_pre_end.tsv')
        if not os.path.exists(paired_file):
            raise FileNotFoundError(f"Required file not found: {paired_file}")
        
        df_paired = pd.read_csv(paired_file, sep="\t")
        
        # Filter by test_type and q_value
        df_filtered = df_paired[
            (df_paired['test_type'] == 'two_sample_paired_tTest') &
            (df_paired['q_value'] < params.q_threshold)
        ].copy()
        
        workflow_logger.info(f"Filtered {df_filtered.shape[0]:,} positions with two_sample_paired_tTest q_value < {params.q_threshold}")
        
        # Select base columns and keep the filtered positions
        final_df = df_filtered[merge_keys].copy()
        
        # Process remaining p-value files, pivot to wide format, and merge with filtered positions
        merged = _process_and_pivot_pvalue_files(
            data_directory,
            p_value_files,
            merge_keys,
            final_df,
            q_threshold=params.q_threshold,
            group=None,
        )

        # Rename columns as specified and write output
        merged.rename(columns={'mag_id': 'MAG', 'gene_id': 'gene'}, inplace=True)
        Path(output.filtered).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output.filtered, sep='\t', index=False, compression='gzip')
        workflow_logger.info(f"Wrote merged table with {merged.shape[0]:,} positions and {merged.shape[1]:,} columns at {output.filtered}")

rule bh_pValues_single:
    """
    Filter p_summary for single_sample_tTest with q_value < 0.05 for each group,
    then merge p-values from LMM across time and CMH across time for those positions.
    """
    input:
        p_value_dir=config["p_value_summary_dir"],
    output:
        filtered=os.path.join(BASE_OUTPUT_DIR, "bh_pValues", "single_sample_{group}_bh_pValues.tsv.gz"),
    params:
        q_threshold=0.05,
        group="{group}",
    threads: 1
    run:
        data_directory = input.p_value_dir
        
        # Define p-value files for single sample analysis
        p_value_files = [
            os.path.join(data_directory, 'p_value_summary_lmm_across_time_pre_end.tsv'),
            os.path.join(data_directory, 'p_value_summary_single_sample_pre_end.tsv'),
            os.path.join(data_directory, 'p_value_summary_cmh_across_time_pre_end.tsv'),
        ]
        
        # Define the core columns used to join the dataframes
        merge_keys = ['mag_id', 'contig', 'position', 'gene_id']
        
        # --- Step 1: Filter for single_sample_tTest with q_value < 0.05 for this group ---
        single_file = os.path.join(data_directory, 'p_value_summary_single_sample_pre_end.tsv')
        if not os.path.exists(single_file):
            raise FileNotFoundError(f"Required file not found: {single_file}")
        
        df_single = pd.read_csv(single_file, sep="\t")
        
        # Filter by test_type, group, and q_value
        df_filtered = df_single[
            (df_single['test_type'] == 'single_sample_tTest') &
            (df_single['group_analyzed'] == params.group) &
            (df_single['q_value'] < params.q_threshold)
        ].copy()
        
        workflow_logger.info(f"Filtered {df_filtered.shape[0]:,} positions for group {params.group} with single_sample_tTest q_value < {params.q_threshold}")
        
        # Select base columns and keep the filtered positions
        final_df = df_filtered[merge_keys].copy()
        
        # Process remaining p-value files, pivot to wide format, and merge with filtered positions
        merged = _process_and_pivot_pvalue_files(
            data_directory,
            p_value_files,
            merge_keys,
            final_df,
            q_threshold=params.q_threshold,
            group=params.group,
        )

        # Rename columns as specified and write output
        merged.rename(columns={'mag_id': 'MAG', 'gene_id': 'gene'}, inplace=True)
        Path(output.filtered).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output.filtered, sep='\t', index=False, compression='gzip')
        workflow_logger.info(f"Wrote merged table for group {params.group} with {merged.shape[0]:,} positions and {merged.shape[1]:,} columns at {output.filtered}")

# --- Step 2: Run Quality Control ---

rule quality_control_paired:
    """
    Runs quality_control.py on the positions extracted from two-sample paired tests.
    """
    input:
        script=os.path.join(os.environ["HOME"], "AlleleFlux/alleleflux/scripts/preprocessing/quality_control.py"),
        positions=PAIRED_POSITIONS_FILE,
        metadata_dir=config["metadata_dir"],
        fasta=config["fasta"],
        mag_mapping=config["mag_mapping"],
    output:
        PAIRED_QC_SUMMARY,
    params:
        outdir=os.path.dirname(PAIRED_QC_SUMMARY),
        breadth=config["breadth_threshold"],
    threads: config["cpus"]
    shell:
        """
        python {input.script} \
            --rootDir {input.metadata_dir} \
            --fasta {input.fasta} \
            --mag_mapping_file {input.mag_mapping} \
            --positions_file {input.positions} \
            --output_dir {params.outdir} \
            --breadth_threshold {params.breadth} \
            --cpus {threads} \
            --positions_denominator positions \
            --data_type longitudinal
        """

rule quality_control_single:
    """
    Runs quality_control.py on the positions extracted from single-sample tests,
    one for each group.
    """
    input:
        script=os.path.join(os.environ["HOME"], "AlleleFlux/alleleflux/scripts/preprocessing/quality_control.py"),
        positions=os.path.join(BASE_OUTPUT_DIR, "tested_positions", "single_sample_{group}_mag_positions.tsv.gz"),
        metadata_dir=config["metadata_dir"],
        fasta=config["fasta"],
        mag_mapping=config["mag_mapping"],
    output:
        os.path.join(BASE_OUTPUT_DIR, "qc_results_single_{group}", "ALL_MAGs_QC_overall_summary.tsv"),
    params:
        group="{group}",
        outdir=os.path.join(BASE_OUTPUT_DIR, "qc_results_single_{group}"),
        breadth=config["breadth_threshold"],
    threads: config["cpus"]
    shell:
        """
        python {input.script} \
            --rootDir {input.metadata_dir} \
            --fasta {input.fasta} \
            --mag_mapping_file {input.mag_mapping} \
            --positions_file {input.positions} \
            --output_dir {params.outdir} \
            --breadth_threshold {params.breadth} \
            --cpus {threads} \
            --data_type single
        """

# --- Step 3: Compute Coverage and Allele Statistics per MAG ---

rule compute_mean_coverage:
    """
    Compute coverage and allele statistics for a single MAG.
    Requires metadata_dir and optionally uses qc_dir for filtering.
    """
    input:
        metadata_dir=config["metadata_dir"],
        qc_dir=config["qc_dir"]
    output:
        mean_coverage_out=os.path.join(BASE_OUTPUT_DIR, "coverage_stats", "{mag}_mean_coverage.tsv"),

    threads: 1
    shell:
        """
        alleleflux-coverage-allele-stats \
            --rootDir {input.metadata_dir} \
            --output_dir {BASE_OUTPUT_DIR}/coverage_stats \
            --mag_id {wildcards.mag} \
            --qc_dir {input.qc_dir}
        """

# --- Step 4: Create Megatables by Merging Coverage Stats with BH P-values ---

rule create_merged_table_paired:
    """
    Merge coverage statistics with BH-corrected p-values for paired test positions.
    Subsets coverage stats to only positions present in the p-value table.
    """
    input:
        bh_pvalues=PAIRED_BH_pValues,
        coverage_stats=COVERAGE_STATS_FILES,
    output:
        merged_table=PAIRED_MERGED_TABLE_1,
    threads: 1
    resources:
        mem_mb=400000
    run:
        _merge_coverage_with_pvalues(
            bh_pvalues_file=input.bh_pvalues,
            coverage_stats_dir=os.path.join(BASE_OUTPUT_DIR, "coverage_stats"),
            output_file=output.merged_table,
            group_label=None
        )

rule create_merged_table_single:
    """
    Merge coverage statistics with BH-corrected p-values for single-sample test positions.
    Subsets coverage stats to only positions present in the p-value table for each group.
    """
    input:
        bh_pvalues=os.path.join(BASE_OUTPUT_DIR, "bh_pValues", "single_sample_{group}_bh_pValues.tsv.gz"),
        coverage_stats=COVERAGE_STATS_FILES,
    output:
        merged_table=os.path.join(BASE_OUTPUT_DIR, "pValue_stats_merged", "single_sample_{group}_merged_table.tsv.gz"),
    params:
        group="{group}",
    resources:
        mem_mb=400000
    threads: 1
    run:
        _merge_coverage_with_pvalues(
            bh_pvalues_file=input.bh_pvalues,
            coverage_stats_dir=os.path.join(BASE_OUTPUT_DIR, "coverage_stats"),
            output_file=output.merged_table,
            group_label=params.group
        )

# --- Step 5: Create Final Megatables with QC Summary Data ---

rule create_megatable_paired:
    """
    Merge QC summary data (overall, group, and group_time) with the paired merged table.
    Creates the final comprehensive megatable for paired analysis.
    """
    input:
        merged_table=PAIRED_MERGED_TABLE_1,
        qc_summary=PAIRED_QC_SUMMARY,
    output:
        megatable=PAIRED_MEGATABLE,
    params:
        qc_dir=os.path.dirname(PAIRED_QC_SUMMARY),
    threads: 1
    resources:
        mem_mb=400000
    run:
        # Load the merged table (p-values + coverage stats)
        df_merged = pd.read_csv(input.merged_table, sep='\t')
        workflow_logger.info(f"Loaded merged table with {df_merged.shape[0]:,} positions and {df_merged.shape[1]:,} columns")
        
        # Load and flatten QC summaries
        df_qc = _load_and_flatten_qc_summaries(params.qc_dir)
        workflow_logger.info(f"Loaded and flattened QC summaries with {df_qc.shape[0]:,} MAGs and {df_qc.shape[1]:,} columns")
        
        # Rename MAG_ID to MAG in QC data for consistent merging
        if 'MAG_ID' in df_qc.columns:
            df_qc.rename(columns={'MAG_ID': 'MAG'}, inplace=True)
        
        # Merge on MAG
        megatable = df_merged.merge(df_qc, on='MAG', how='left', suffixes=('', '_qc'))
        
        workflow_logger.info(f"Created final megatable with {megatable.shape[0]:,} positions and {megatable.shape[1]:,} columns")
        
        # Write output
        Path(output.megatable).parent.mkdir(parents=True, exist_ok=True)
        megatable.to_csv(output.megatable, sep='\t', index=False, compression='gzip')
        workflow_logger.info(f"Wrote final megatable to {output.megatable}")

rule create_megatable_single:
    """
    Merge QC summary data (overall and group) with the single-sample merged table.
    Creates the final comprehensive megatable for single-sample analysis per group.
    """
    input:
        merged_table=os.path.join(BASE_OUTPUT_DIR, "pValue_stats_merged", "single_sample_{group}_merged_table.tsv.gz"),
        qc_summary=os.path.join(BASE_OUTPUT_DIR, "qc_results_single_{group}", "ALL_MAGs_QC_overall_summary.tsv"),
    output:
        megatable=os.path.join(BASE_OUTPUT_DIR, "megatables", "single_sample_{group}_megatable.tsv.gz"),
    params:
        group="{group}",
        qc_dir=os.path.join(BASE_OUTPUT_DIR, "qc_results_single_{group}"),
    threads: 1
    resources:
        mem_mb=400000
    run:
        # Load the merged table (p-values + coverage stats)
        df_merged = pd.read_csv(input.merged_table, sep='\t', compression='gzip')
        workflow_logger.info(f"Loaded merged table for group {params.group} with {df_merged.shape[0]:,} positions and {df_merged.shape[1]:,} columns")
        
        # Load and flatten QC summaries
        df_qc = _load_and_flatten_qc_summaries(params.qc_dir)
        workflow_logger.info(f"Loaded and flattened QC summaries for group {params.group} with {df_qc.shape[0]:,} MAGs and {df_qc.shape[1]:,} columns")
        
        # Rename MAG_ID to MAG in QC data for consistent merging
        if 'MAG_ID' in df_qc.columns:
            df_qc.rename(columns={'MAG_ID': 'MAG'}, inplace=True)
        
        # Merge on MAG
        megatable = df_merged.merge(df_qc, on='MAG', how='left', suffixes=('', '_qc'))
        
        workflow_logger.info(f"Created final megatable for group {params.group} with {megatable.shape[0]:,} positions and {megatable.shape[1]:,} columns")
        
        # Write output
        Path(output.megatable).parent.mkdir(parents=True, exist_ok=True)
        megatable.to_csv(output.megatable, sep='\t', index=False, compression='gzip')
        workflow_logger.info(f"Wrote final megatable to {output.megatable}")
