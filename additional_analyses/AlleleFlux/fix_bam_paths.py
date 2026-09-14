#!/usr/bin/env python3
"""
Rewrite the `bam_path` column of the diet-manipulation metadata sheet.

The original metadata (`metadata_md_bam.tsv`) points at
`/scratch/gpfs/sg4230/popgentoolkit/bt2/`, which no longer exists. The BAMs now
live under the AMOELLER project copy. The old AlleleFlux config worked around
this with `input.bam_dir`; the current schema has no such key and validates
every `bam_path` at DAG construction, so the sheet itself must be corrected.

Only the `bam_path` column changes; every other column is copied verbatim.

Usage:
    python fix_bam_paths.py                      # defaults below
    python fix_bam_paths.py --bam_dir DIR --output OUT.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

DEFAULT_METADATA = "/scratch/gpfs/AMOELLER/diet_manip/copy_sg4230_scratch/sg4230/popgentoolkit/metadata_md_bam.tsv"
DEFAULT_BAM_DIR = "/scratch/gpfs/AMOELLER/diet_manip/copy_sg4230_scratch/sg4230/popgentoolkit/bt2"
DEFAULT_OUTPUT = str(Path(__file__).parent / "metadata_md_bam_fixed.tsv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=DEFAULT_METADATA, help="Input metadata TSV")
    parser.add_argument("--bam_dir", default=DEFAULT_BAM_DIR, help="Directory holding the BAM files")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Output metadata TSV")
    parser.add_argument(
        "--extension",
        default=".sorted.bam",
        help="BAM filename suffix appended to sample_id (default: .sorted.bam)",
    )
    parser.add_argument(
        "--require_index",
        action="store_true",
        help="Also require a .bai next to each BAM",
    )
    args = parser.parse_args()

    bam_dir = Path(args.bam_dir)
    if not bam_dir.is_dir():
        sys.exit(f"ERROR: bam_dir is not a directory: {bam_dir}")

    df = pd.read_csv(args.metadata, sep="\t")
    if "sample_id" not in df.columns:
        sys.exit("ERROR: metadata is missing the required 'sample_id' column")

    df["bam_path"] = df["sample_id"].map(lambda s: str(bam_dir / f"{s}{args.extension}"))

    missing = [p for p in df["bam_path"] if not Path(p).exists()]
    if missing:
        sys.exit(
            f"ERROR: {len(missing)}/{len(df)} BAM files not found, e.g. {missing[:3]}\n"
            "Check --bam_dir and --extension."
        )

    if args.require_index:
        no_idx = [p for p in df["bam_path"] if not Path(p + ".bai").exists()]
        if no_idx:
            sys.exit(f"ERROR: {len(no_idx)}/{len(df)} BAMs have no .bai index, e.g. {no_idx[:3]}")

    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.output}")
    print(f"All {len(df)} bam_path entries verified to exist under {bam_dir}")


if __name__ == "__main__":
    main()
