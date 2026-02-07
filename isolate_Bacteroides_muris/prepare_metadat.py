#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Convert Moeller Isolate Strain Library TSV -> AlleleFlux metadata TSV."
    )
    ap.add_argument(
        "--input", required=True, help="Input TSV exported from the isolate sheet."
    )
    ap.add_argument(
        "--output", required=True, help="Output metadata TSV for AlleleFlux."
    )
    ap.add_argument("--bam-dir", required=True, help="Directory containing BAMs.")
    args = ap.parse_args()

    logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.INFO)

    in_path = Path(args.input)
    out_path = Path(args.output)
    bam_dir = Path(args.bam_dir)

    df = pd.read_csv(in_path, sep="\t", dtype=str)
    logging.info(f"Loaded {len(df)} rows from {in_path}")

    # Keep only rows that actually have a genome ID (everything else can't become a sample)
    df = df[df["User.Genome"].str.strip() != ""].copy()
    logging.info(f"Rows with non-empty User.Genome: {len(df)}")

    df["sample_id"] = (
        df["User.Genome"].str.strip().str.replace(r"\.fna$", "", regex=True)
    )
    df.rename(
        columns={
            "Host.Diet": "group",
            "Host.Collection.Timepoint": "time",
            "Host.Litter.Replicate": "replicate",
            "Host.ID": "subjectID",
        },
        inplace=True,
    )

    # Construct bam_path
    df["bam_path"] = df["sample_id"].map(
        lambda s: str((bam_dir / f"{s}.sorted.bam").resolve())
    )

    meta = df[
        ["sample_id", "bam_path", "group", "time", "replicate", "subjectID"]
    ].copy()

    # Drop any broken rows (should be rare)
    before = len(meta)
    meta = meta[
        (meta["sample_id"].str.strip() != "")
        & (meta["group"].str.strip() != "")
        & (meta["time"].str.strip() != "")
    ].copy()
    logging.info(f"Dropped {before - len(meta)} rows missing sample/group/time")

    # Check BAM existence and remove rows with missing BAMs
    before_bams = len(meta)
    bam_exists = meta["bam_path"].map(lambda p: Path(p).exists())
    missing_bam_rows = meta[~bam_exists].copy()
    meta = meta[bam_exists].copy()
    bam_dropped = before_bams - len(meta)
    if bam_dropped > 0:
        logging.warning(f"Dropped {bam_dropped} rows with missing BAMs:")
        for sample_id in missing_bam_rows["sample_id"].tolist():
            logging.warning(f"  {sample_id}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out_path, sep="\t", index=False)
    logging.info(f"Wrote AlleleFlux metadata: {out_path} ({len(meta)} rows)")


if __name__ == "__main__":
    main()
