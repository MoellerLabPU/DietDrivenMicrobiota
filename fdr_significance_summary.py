#!/usr/bin/env python3
"""
Compute BH-corrected minima and counts for each test-type across all files for a
given pair of timepoints, but only for mag_ids present in the 'pre_end'
paired_tTest set.

This improved script takes a single pair of timepoints (e.g., 'pre' and 'post')
as input, and produces a wide summary TSV and a detailed TSV listing every
significant row for that specific comparison. It is more efficient by reading
files only once and includes detailed logging.
"""

import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# Centralized configuration for all statistical tests to be processed.
# This makes it easier to add, remove, or modify tests without changing the logic.
TEST_CONFIGS = [
    {
        "name": "paired_tTest",
        "dir_template": "two_sample_paired_{period}-fat_control",
        "suffix": "_two_sample_paired.tsv.gz",
        "p_col_regex": r"_p_value_tTest$",
        "is_master_for_mag_ids": True,  # Use this test to define the initial set of mag_ids
    },
    {
        "name": "paired_Wilcoxon",
        "dir_template": "two_sample_paired_{period}-fat_control",
        "suffix": "_two_sample_paired.tsv.gz",
        "p_col_regex": r"_p_value_Wilcoxon$",
    },
    {
        "name": "control_tTest",
        "dir_template": "single_sample_{period}-fat_control",
        "suffix": "_single_sample_control.tsv.gz",
        "p_col_regex": r"_p_value_tTest_control$",
    },
    {
        "name": "control_Wilcoxon",
        "dir_template": "single_sample_{period}-fat_control",
        "suffix": "_single_sample_control.tsv.gz",
        "p_col_regex": r"_p_value_Wilcoxon_control$",
    },
    {
        "name": "fat_tTest",
        "dir_template": "single_sample_{period}-fat_control",
        "suffix": "_single_sample_fat.tsv.gz",
        "p_col_regex": r"_p_value_tTest_fat$",
    },
    {
        "name": "fat_Wilcoxon",
        "dir_template": "single_sample_{period}-fat_control",
        "suffix": "_single_sample_fat.tsv.gz",
        "p_col_regex": r"_p_value_Wilcoxon_fat$",
    },
    {
        "name": "lmm",
        "dir_template": "lmm_{period}-fat_control",
        "suffix": "_lmm.tsv.gz",
        "p_col_regex": r"_p_value_LMM$",
    },
    {
        "name": "lmm_across_time_control",
        "dir_template": "lmm_across_time_{period}-fat_control",
        "suffix": "_lmm_across_time_control.tsv.gz",
        "p_col_regex": r"_p_value_LMM$",
    },
    {
        "name": "lmm_across_time_fat",
        "dir_template": "lmm_across_time_{period}-fat_control",
        "suffix": "_lmm_across_time_fat.tsv.gz",
        "p_col_regex": r"_p_value_LMM$",
    },
]

# Columns required for the detailed output file.
METADATA_COLS = ["contig", "gene_id", "position"]


def extract_data_from_files(
    file_paths: Dict[str, Path], p_col_regex: str
) -> pd.DataFrame:
    """
    Reads data files efficiently, extracting row-wise minimum p-values and metadata.
    This version first reads only the header to find columns, then reads only
    the necessary data to conserve memory.

    Args:
        file_paths: A dictionary mapping mag_ids to their full file paths (as Path objects).
        p_col_regex: The regex to identify p-value columns.

    Returns:
        A DataFrame containing all data, with 'mag_id' and 'min_p_value' columns.
    """
    all_data = []
    for mag_id, path in file_paths.items():
        # 1. Read only header to identify p-value columns dynamically
        header = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns
        p_value_cols = [c for c in header if re.search(p_col_regex, c)]

        if not p_value_cols:
            logging.warning(
                f"No p-value columns found in {path} for regex '{p_col_regex}'"
            )
            continue

        # 2. Read only the columns we absolutely need
        use_cols = METADATA_COLS + p_value_cols
        df = pd.read_csv(path, sep="\t", compression="gzip", usecols=use_cols)

        # 3. Create the output frame with metadata and mag_id
        df_out = df[METADATA_COLS].copy()
        df_out["mag_id"] = mag_id

        # 4. Calculate row-wise minimum p-value efficiently
        p_values_numeric = df[p_value_cols].apply(pd.to_numeric, errors="coerce")
        df_out["min_p_value"] = p_values_numeric.min(axis=1, skipna=True)

        all_data.append(df_out)

    if not all_data:
        logging.warning(
            "No data extracted from files. Ensure the provided paths and regex are correct."
        )
        return pd.DataFrame()

    return pd.concat(all_data, ignore_index=True)


def summarize_period(root_dir: Path, period: str, mag_ids: Set[str], out_path: Path):
    """
    Orchestrates the analysis for a single period (e.g., 'pre_end').
    Gathers data, performs BH correction, and writes summary and detail files.
    """
    logging.info(f"Starting analysis for period: {period}")

    all_significant_rows = []

    # Initialize a nested dictionary to hold summary results for each mag_id
    summary_results = {base: {} for base in mag_ids}

    # This is needed for the 'lmm_pre' test which has a different dir name structure
    # suffix_period = period.split("_", 1)[1]

    for test_config in TEST_CONFIGS:
        test_name = test_config["name"]

        # Build the subdirectory path from the template
        dir_template = test_config["dir_template"]
        subdir = dir_template.format(period=period)
        dirpath = root_dir / subdir

        if not dirpath.is_dir():
            logging.warning(f"Test directory not found, skipping: {dirpath}")
            continue

        logging.info(f"Processing test '{test_name}' for period '{period}'...")

        # Find all relevant files for the current test that match our master list of mag_ids
        suffix = test_config["suffix"]
        file_paths = {
            path.name[: -len(suffix)]: path
            for path in sorted(dirpath.glob(f"*{suffix}"))
            if path.name[: -len(suffix)] in mag_ids
        }

        if not file_paths:
            logging.warning(
                f"No files matching mag_ids found for test '{test_name}' in {dirpath}"
            )
            continue

        logging.info(f"Found {len(file_paths)} files for test '{test_name}'.")

        # 1. Read all data and extract p-values (one read pass per file)
        combined_df = extract_data_from_files(file_paths, test_config["p_col_regex"])
        if combined_df.empty:
            continue

        # 2. Perform multiple testing correction on all valid p-values at once
        valid_p_values_df = combined_df.dropna(subset=["min_p_value"])
        if valid_p_values_df.empty:
            logging.warning(f"No valid p-values found for test '{test_name}'.")
            continue

        logging.info(
            f"Correcting {len(valid_p_values_df):,} p-values for test '{test_name}'..."
        )

        _, qvals, _, _ = multipletests(
            valid_p_values_df["min_p_value"], method="fdr_bh"
        )

        # Add q-values back to the DataFrame
        corrected_df = valid_p_values_df.copy()
        corrected_df["q_value"] = qvals

        # 3. Identify and store all significant rows for the detailed report
        significant_rows = corrected_df[corrected_df["q_value"] < 0.05].copy()
        if not significant_rows.empty:
            significant_rows["period"] = period
            significant_rows["test"] = test_name
            all_significant_rows.append(significant_rows)

        # 4. Generate summary statistics (min q-value and count of significant rows) using groupby
        summary = (
            corrected_df.groupby("mag_id")
            .agg(
                min_q=("q_value", "min"), n_sig=("q_value", lambda x: (x < 0.05).sum())
            )
            .reindex(sorted(mag_ids))
        )  # Ensure all mag_ids are present in the summary

        # 5. Store summary results in our main results dictionary
        for mag_id, row in summary.iterrows():
            summary_results[mag_id][f"{test_name}_min_q"] = row["min_q"]
            summary_results[mag_id][f"{test_name}_n_sig"] = row["n_sig"]

    detail_path = out_path.with_name(f"{out_path.stem}_significant_rows.tsv")
    # Write detailed significant rows file
    if all_significant_rows:
        detailed_df = pd.concat(all_significant_rows, ignore_index=True)
        detailed_df.fillna(np.nan, inplace=True)
        # Reorder columns for the output file
        cols = ["period", "test", "mag_id"] + METADATA_COLS + ["q_value"]
        detailed_df[cols].to_csv(
            detail_path, sep="\t", index=False, float_format="%.4g"
        )
        logging.info(f"Wrote {len(detailed_df)} significant rows to {detail_path}")
    else:
        # Create an empty file with headers if no significant rows were found
        logging.info(
            f"No significant rows found for period '{period}'. Writing empty file to {detail_path}"
        )
        with open(detail_path, "w") as f:
            f.write(
                "\t".join(["period", "test", "mag_id"] + METADATA_COLS + ["q_value"])
                + "\n"
            )

    # Build and write the wide summary DataFrame
    summary_rows = []
    for mag_id in sorted(mag_ids):
        row = {"mag_id": mag_id, **summary_results[mag_id]}
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    # Fill any missing values with NaN for clarity in the final report
    summary_df.fillna(np.nan, inplace=True)
    summary_df.to_csv(out_path, sep="\t", index=False)
    logging.info(f"Wrote {period} summary for {len(summary_df)} mag_ids to {out_path}")


def main():
    """Main function to parse arguments and run the analysis."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    p = argparse.ArgumentParser(
        description="Efficiently compute BH-corrected p-values and summarize results for a given pair of timepoints.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--root_dir",
        required=True,
        type=Path,
        help="Path to the 'significance_tests/' folder.",
    )
    p.add_argument(
        "--time1", required=True, help="First timepoint for comparison (e.g., 'pre')."
    )
    p.add_argument(
        "--time2", required=True, help="Second timepoint for comparison (e.g., 'end')."
    )
    p.add_argument(
        "--out_dir", required=True, type=Path, help="Path for the output directory."
    )
    args = p.parse_args()

    # Create the output directory if it doesn't exist
    args.out_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Starting analysis.")

    # Determine master set of mag_ids from the test marked as master in the config.
    # This ensures that we only analyze a consistent set of samples across all tests,
    # based on the samples present in the 'pre_end' comparison.
    try:
        master_test_config = next(
            c for c in TEST_CONFIGS if c.get("is_master_for_mag_ids")
        )
    except StopIteration:
        logging.error(
            "No test in TEST_CONFIGS is marked with 'is_master_for_mag_ids: True'. Cannot determine mag_ids. Exiting."
        )
        return

    # The directory for mag_ids is fixed to pre_end, as per the script's original purpose.
    pe_dir = args.root_dir / master_test_config["dir_template"].format(period="pre_end")
    master_suffix = master_test_config["suffix"]

    try:
        if not pe_dir.is_dir():
            raise FileNotFoundError
        mag_ids = {
            path.name[: -len(master_suffix)]
            for path in pe_dir.glob(f"*{master_suffix}")
        }
        if not mag_ids:
            logging.error(
                f"No mag_ids found in master directory {pe_dir} using suffix '{master_suffix}'. Exiting."
            )
            return
        logging.info(
            f"Found {len(mag_ids)} master mag_ids to analyze, based on the 'pre_end' comparison."
        )
    except FileNotFoundError:
        logging.error(f"Master directory for mag_ids not found: {pe_dir}. Exiting.")
        return

    # Construct the period string and the output file path
    period = f"{args.time1}_{args.time2}"
    output_file = args.out_dir / f"{period}_summary.tsv"

    logging.info(f"Running analysis for period: {period}")
    logging.info(f"Output will be saved to: {output_file}")

    # Run the analysis for the single specified period
    summarize_period(args.root_dir, period, mag_ids, output_file)

    logging.info("Analysis complete.")


if __name__ == "__main__":
    main()
