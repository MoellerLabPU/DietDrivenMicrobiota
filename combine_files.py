#!/usr/bin/env python3
# pylint: disable=logging-fstring-interpolation

import argparse
import glob
import logging
import os
import sys
from pathlib import Path

import pandas as pd


def process_single_sample_type(
    single_dir,
    found_mag_ids,
    output_dir,
    file_suffix,
    output_filename,
    sample_type_label,
):
    """
    Finds, filters, and concatenates single-sample files of a specific type.

    Args:
        single_dir (str): The directory path containing the single-sample files.
        found_mag_ids (set): A set of MAG IDs to filter by.
        output_dir (str): The directory where the concatenated output file will be saved.
        file_suffix (str): The suffix of the files to process (e.g., '_single_sample_control.tsv.gz').
        output_filename (str): The name for the output concatenated TSV file.
        sample_type_label (str): A label for logging purposes (e.g., 'control').
    """
    logging.info(f"Processing {sample_type_label} files from: {single_dir}")
    files_to_process = glob.glob(os.path.join(single_dir, f"*{file_suffix}"))
    list_of_dfs = []
    processed_mag_ids_for_type = set()

    for file_path in files_to_process:
        filename = os.path.basename(file_path)
        mag_id = filename.removesuffix(file_suffix)
        # Check for duplicate MAG IDs for this specific sample type
        if mag_id in processed_mag_ids_for_type:
            logging.error(
                f"Duplicate MAG ID '{mag_id}' found for {sample_type_label} samples. Please ensure each MAG ID is unique. Offending file: {file_path}"
            )
            sys.exit(1)

        # Only process this file if its MAG ID was in the paired-sample set
        if mag_id in found_mag_ids:
            processed_mag_ids_for_type.add(mag_id)
            df = pd.read_csv(file_path, sep="\t", compression="gzip")
            df["MAG ID"] = mag_id
            list_of_dfs.append(df)

    # Concatenate and save the files
    if list_of_dfs:
        concatenated_df = pd.concat(list_of_dfs, ignore_index=True)
        output_path = os.path.join(output_dir, output_filename)
        concatenated_df.to_csv(output_path, sep="\t", index=False)
        logging.info(
            f"Successfully created concatenated {sample_type_label}-sample file at: {output_path}"
        )
    else:
        logging.warning(f"No matching {sample_type_label} files found for the MAG IDs.")


def process_files(paired_dir, single_dir, output_dir):
    """
    Processes and concatenates MAG sequencing data from paired and single-sample files.

    This function first finds all two-sample paired files, concatenates them,
    and extracts a set of unique MAG IDs. It then processes single-sample
    'control' and 'fat' files, filtering them to include only the MAG IDs
    found in the paired-sample files.

    Args:
        paired_dir (str): The directory path containing the two-sample paired files
                          (e.g., *_two_sample_paired.tsv.gz).
        single_dir (str): The directory path containing the single-sample files
                          (e.g., *_single_sample_control.tsv.gz and
                          *_single_sample_fat.tsv.gz).
        output_dir (str): The directory where the concatenated output files
                          will be saved.
    """
    # --- Step 1: Create output directory if it doesn't exist ---
    if not os.path.exists(output_dir):
        logging.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

    # --- Step 2: Process the two-sample paired files ---
    logging.info(f"Processing paired-sample files from: {paired_dir}")
    paired_files = glob.glob(os.path.join(paired_dir, "*_two_sample_paired.tsv.gz"))

    if not paired_files:
        logging.error(f"No '*_two_sample_paired.tsv.gz' files found in '{paired_dir}'.")
        # Exit if the primary files are not found, as MAG IDs are derived from them.
        sys.exit(1)

    list_of_paired_dfs = []
    found_mag_ids = set()

    for file_path in paired_files:
        # Extract the base filename to derive the MAG ID
        filename = os.path.basename(file_path)
        # The MAG ID is the filename minus the specific suffix
        mag_id = filename.removesuffix("_two_sample_paired.tsv.gz")

        # Check for duplicate MAG IDs
        if mag_id in found_mag_ids:
            logging.error(
                f"Duplicate MAG ID '{mag_id}' found in paired-sample files. Please ensure each MAG ID is unique. Offending file: {file_path}"
            )
            sys.exit(1)

        # Store the MAG ID for later use
        found_mag_ids.add(mag_id)

        # Read the gzipped tsv file into a pandas DataFrame
        df = pd.read_csv(file_path, sep="\t", compression="gzip")

        # Add the new 'MAG ID' column
        df["MAG ID"] = mag_id

        list_of_paired_dfs.append(df)

    # Concatenate all paired-sample dataframes into one
    if list_of_paired_dfs:
        concatenated_paired_df = pd.concat(list_of_paired_dfs, ignore_index=True)

        # Save the final concatenated dataframe
        output_path_paired = os.path.join(output_dir, "concatenated_paired_samples.tsv")
        concatenated_paired_df.to_csv(output_path_paired, sep="\t", index=False)
        logging.info(
            f"Successfully created concatenated paired-sample file at: {output_path_paired}"
        )
        logging.info(f"Found {len(found_mag_ids)} unique MAG IDs.")
    else:
        logging.error("No paired-sample files were successfully processed. Exiting.")
        return

    # --- Steps 3 & 4: Process single-sample files using the refactored function ---

    # Process the 'control' files
    process_single_sample_type(
        single_dir=single_dir,
        found_mag_ids=found_mag_ids,
        output_dir=output_dir,
        file_suffix="_single_sample_control.tsv.gz",
        output_filename="concatenated_single_sample_control.tsv",
        sample_type_label="control",
    )

    # Process the 'fat' files
    process_single_sample_type(
        single_dir=single_dir,
        found_mag_ids=found_mag_ids,
        output_dir=output_dir,
        file_suffix="_single_sample_fat.tsv.gz",
        output_filename="concatenated_single_sample_fat.tsv",
        sample_type_label="fat",
    )


def main():
    """
    Main function to parse arguments and run the file processing.
    """
    # --- Setup Logging ---
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        stream=sys.stdout,
    )

    # --- Setup Argument Parser for Command-Line usage ---
    parser = argparse.ArgumentParser(
        description="Concatenate MAG sequencing data files based on MAG IDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--paired_dir",
        type=Path,
        required=True,
        help="Directory containing the two-sample paired files (*_two_sample_paired.tsv.gz).",
    )

    parser.add_argument(
        "--single_dir",
        type=Path,
        required=True,
        help="Directory containing the single-sample files (*_single_sample_control.tsv.gz and *_single_sample_fat.tsv.gz).",
    )

    parser.add_argument(
        "--output_dir",
        type=Path,
        required=True,
        help="Directory where the output concatenated files will be saved.",
    )

    args = parser.parse_args()

    # Run the main function with the provided arguments
    process_files(args.paired_dir, args.single_dir, args.output_dir)


if __name__ == "__main__":
    main()
