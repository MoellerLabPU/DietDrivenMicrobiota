#!/usr/bin/env python3
"""
Extract a unified table of MAG, position, contig from per-MAG result files.

Given an input directory and a test type, this script finds all matching
TSV.GZ files, infers the MAG ID from the filename, and concatenates rows
from the columns [contig, position] into a single table with columns
[MAG, position, contig].

Example filenames:
  - SLG779_DASTool_bins_SLG779_bin.32_two_sample_paired.tsv.gz
  - SLG887_DASTool_bins_SLG887_bin.70_two_sample_paired.tsv.gz
  - SLG779_DASTool_bins_SLG779_bin.32_single_sample_fat.tsv.gz
  - SLG779_DASTool_bins_SLG779_bin.32_single_sample_control.tsv.gz

Recognized test types:
  - two_sample_paired
  - single_sample (matches single_sample_*)

Usage:
  python extract_mag_positions.py \
    --input-dir /path/to/dir \
    --test-type two_sample_paired \
    --out combined_mag_positions.tsv.gz

Notes:
  - Only the columns 'contig' and 'position' are read from input files.
  - Rows from multiple files for the same MAG are combined.
  - Duplicate rows across files are dropped by default (configurable).
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

logging.basicConfig(
    format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

# High-level: this script discovers per-MAG result files, extracts the
# contig/position info, validates uniqueness, and writes combined tables.


def discover_files(input_dir: Path, test_type: str) -> List[Path]:
    """Discover .tsv.gz files in input_dir matching test_type.

    For test_type == 'two_sample_paired': match '*_two_sample_paired.tsv.gz'
    For test_type == 'single_sample':     match '*_single_sample*.tsv.gz'
    """
    # Ensure the provided path exists and is a directory
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input dir not found or not a directory: {input_dir}")

    pattern = {
        "two_sample_paired": "*_two_sample_paired.tsv.gz",
        "single_sample": "*_single_sample*.tsv.gz",
    }.get(test_type)

    if pattern is None:
        raise ValueError(
            "Unsupported --test-type. Use one of: two_sample_paired, single_sample"
        )

    # Return a sorted list of file paths matching the test-type pattern
    return sorted(input_dir.glob(pattern))


def infer_mag_id(file_path: Path, test_type: str) -> Optional[str]:
    """Infer MAG ID from the filename using the specified test_type.

    Examples:
      SLG779_DASTool_bins_SLG779_bin.32_two_sample_paired.tsv.gz -> SLG779_DASTool_bins_SLG779_bin.32
      SLG779_DASTool_bins_SLG779_bin.32_single_sample_fat.tsv.gz -> SLG779_DASTool_bins_SLG779_bin.32
    """
    # Work with the filename only (no directory components)
    name = file_path.name
    # Build regex that captures MAG part before the test type and optional suffix, ending with .tsv.gz
    # For single_sample, allow extra label, e.g., _single_sample_fat
    if test_type == "two_sample_paired":
        regex = rf"^(?P<mag>.+?)_two_sample_paired\.tsv\.gz$"
    elif test_type == "single_sample":
        regex = rf"^(?P<mag>.+?)_single_sample(?:_[^.]+)?\.tsv\.gz$"
    else:
        return None

    # Use regex to capture the MAG identifier preceding the test-type suffix
    mag = re.match(regex, name)
    if mag:
        return mag.group("mag")
    return None


def extract_single_sample_group(file_path: Path) -> str:
    """Extract the group label for single-sample files.

    Examples:
      SLG779_..._single_sample_fat.tsv.gz -> fat
      SLG779_..._single_sample_control.tsv.gz -> control
      SLG779_..._single_sample.tsv.gz -> raises ValueError (no group)
    """
    # Extract the group name following '_single_sample_' in the filename.
    # If the filename doesn't follow the expected pattern, raise an error.
    name = file_path.name
    match = re.match(r"^.+?_single_sample(?:_(?P<group>[^.]+))?\.tsv\.gz$", name)
    if not match:
        raise ValueError(f"Filename does not match single_sample pattern: {name}")
    grp = match.group("group")
    # Require a non-empty group (user requested strict behavior)
    if not grp:
        raise ValueError(f"No group found in filename: {name}")
    return grp


def build_combined_table(files: List[Path], test_type: str) -> pd.DataFrame:
    """Build combined table and validate uniqueness of MAG/position/contig."""
    rows: List[pd.DataFrame] = []
    # Read each file and collect only the needed columns into a list of DataFrames
    for file in files:
        mag = infer_mag_id(file, test_type)
        if not mag:
            # Fail fast if filename format isn't as expected
            raise ValueError(f"Cannot infer MAG ID from filename: {file.name}")

        # Read only the contig and position columns to minimize memory usage
        df = pd.read_csv(
            file,
            sep="\t",
            usecols=["contig", "position"],
            dtype={"contig": "string", "position": "Int64"},
        )
        if df.empty:
            logging.warning(f"File {file.name} has no data rows, skipping")
            continue

        # Attach the MAG identifier as a column and keep the canonical column order
        df = df.assign(MAG=mag)[["MAG", "position", "contig"]]
        rows.append(df)

    if not rows:
        return pd.DataFrame(columns=["MAG", "position", "contig"], dtype="string")

    # Concatenate per-file frames into one combined DataFrame
    combined = pd.concat(rows, ignore_index=True)

    # Check for duplicates across files. We consider a duplicate any repeated
    # (MAG, position, contig) tuple because each such combination must be unique.
    duplicates = combined.duplicated(subset=["MAG", "position", "contig"])
    if duplicates.any():
        dup_rows = combined[duplicates].sort_values(["MAG", "contig", "position"])
        # Provide a helpful error with a small sample of the offending rows
        raise ValueError(
            f"Found {int(duplicates.sum())} duplicate MAG/position/contig entries across files:\n"
            f"{dup_rows.head(20).to_string(index=False)}"
        )

    # Sort the final table to make outputs stable and easier to inspect
    combined = combined.sort_values(["MAG", "contig", "position"], ignore_index=True)
    return combined


def write_output(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    compression = "infer"
    # If extension ends with .gz, pandas will gzip. Otherwise plain TSV.
    df.to_csv(out_path, sep="\t", index=False, compression=compression)


# Small helper: write_output is intentionally simple; it will create parent
# directories if needed and rely on pandas to infer compression from filename.


def main():
    parser = argparse.ArgumentParser(
        description="Combine MAG contig/position across files into one table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing per-MAG .tsv.gz result files",
    )
    parser.add_argument(
        "--test-type",
        required=True,
        choices=["two_sample_paired", "single_sample"],
        help="Test type to filter files by",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("./out"),
        help=(
            "Output directory. For two_sample_paired a single file is written. "
            "For single_sample separate files are created per group inside this directory."
        ),
    )
    args = parser.parse_args()

    # Find all files matching the requested test type in the input directory
    files = discover_files(args.input_dir, args.test_type)
    if not files:
        raise ValueError(
            f"No files found in {args.input_dir} for test type '{args.test_type}'"
        )
    logging.info(f"Discovered {len(files)} file(s).")

    # Partition files by group. For two_sample_paired, there's a single implicit group (None).
    # For single_sample, we extract the explicit group name from each filename.
    groups: Dict[str, List[Path]] = {}
    if args.test_type == "two_sample_paired":
        # All files belong to a single unnamed group
        groups[""] = files
    elif args.test_type == "single_sample":
        # single_sample: partition files by their group label
        for file in files:
            grp = extract_single_sample_group(file)
            groups.setdefault(grp, []).append(file)

    if not groups:
        raise ValueError(f"No groups discovered for test type '{args.test_type}'")

    # Ensure output directory exists before writing files
    args.out.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    # Process each group: build combined table and write to disk
    for grp in sorted(groups.keys()):
        grp_files = groups[grp]
        group_label = f"'{grp}'" if grp else "all files"
        logging.info(
            f"Group {group_label}: {len(grp_files)} file(s). Building table..."
        )

        combined = build_combined_table(grp_files, args.test_type)

        if combined.empty:
            logging.warning(
                f"No rows identified for group {group_label}. All files are empty."
            )
        # If this is a single-sample run, add a 'group' column to the table
        if args.test_type == "single_sample":
            # attach the group name (same value for all rows) so downstream
            # consumers can easily filter/identify the group
            combined = combined.assign(group=grp)[
                ["MAG", "position", "contig", "group"]
            ]

        # Filename: <test_type>_<group>_mag_positions.tsv.gz
        # For two_sample_paired (empty group), omit the group suffix
        filename = (
            f"{args.test_type}_{grp}_mag_positions.tsv.gz"
            if grp
            else f"{args.test_type}_mag_positions.tsv.gz"
        )
        out_path = args.out / filename
        write_output(combined, out_path)
        logging.info(f"Wrote {len(combined):,} rows to {out_path}")
        total_rows += len(combined)


if __name__ == "__main__":
    main()
