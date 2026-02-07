#!/usr/bin/env python3
"""
subset_bacteroides_muris.py

Pipeline:
1) Read the Moeller Lab Bacterial Isolate Library table (tsv).
2) Subset rows where Isolate.Classification matches "Bacteroides muris".
3) Save the subset table.
4) For each row, derive an assembly name from User.Genome by removing trailing ".fna".
5) Find paired-end reads for each assembly in configured search directories.
6) Copy reads into the destination directory.
7) Verify that all reads were copied successfully (exist + same byte size).
8) Verify that the subset's User.Genome list matches an expected list.

Exit codes:
  0  success
  2  subset does not match expected list
  3  missing reads for at least one assembly
  4  copy verification failed for at least one file
"""
from __future__ import annotations

import argparse
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

log = logging.getLogger(__name__)

CLASS_COL = "Isolate.Classification"
GENOME_COL = "User.Genome"

EXPECTED_USER_GENOMES = [
    "1102_B1_10471829_HT5CLAFX5.fna",
    "1102_B2_10471829_HT5CLAFX5.fna",
    "1102_C2_10471829_HT5CLAFX5.fna",
    "1102_F8_10471829_HT5CLAFX5.fna",
    "1104_A2_10471829_HT5CLAFX5.fna",
    "1104_A5_10471829_HT5CLAFX5.fna",
    "1104_B4_10471829_HT5CLAFX5.fna",
    "1104_B6_10471829_HT5CLAFX5.fna",
    "1104_C3_10471829_HT5CLAFX5.fna",
    "1104_C6_10471829_HT5CLAFX5.fna",
    "1104_D3_10471829_HT5CLAFX5.fna",
    "1104_D4_10471829_HT5CLAFX5.fna",
    "1104_E3_10471829_HT5CLAFX5.fna",
    "1104_F3_10471829_HT5CLAFX5.fna",
    "1104_H4_10471829_HT5CLAFX5.fna",
    "1173_A8_10471829_HT5CLAFX5.fna",
    "1173_B10_10471829_HT5CLAFX5.fna",
    "1173_B9_10471829_HT5CLAFX5.fna",
    "1173_C8_10471829_HT5CLAFX5.fna",
    "1173_E9_10471829_HT5CLAFX5.fna",
    "1173_F8_10471829_HT5CLAFX5.fna",
    "1173_H7_10471829_HT5CLAFX5.fna",
    "1173_H8_10471829_HT5CLAFX5.fna",
    "14229_33255_206112_HN2G7AFX5_Plate_5_B8_GTGCCATA_TACCGGAT.fna",
    "191_G4_10471829_HT5CLAFX5.fna",
    "530_B10_10471829_HT5CLAFX5.fna",
    "530_E10_10471829_HT5CLAFX5.fna",
    "530_E9_10471829_HT5CLAFX5.fna",
    "SLG1107_A7_10473203_HT53MAFX5.fna",
    "SLG1107_D5_10473203_HT53MAFX5.fna",
    "SLG1107_H3_10473203_HT53MAFX5.fna",
    "SLG1107_H5_10473203_HT53MAFX5.fna",
    "SLG1205_D5_10473203_HT53MAFX5.fna",
    "SLG1205_H5_10473203_HT53MAFX5.fna",
    "SLG1207_D6_10473203_HT53MAFX5.fna",
    "SLG1207_F6_10473203_HT53MAFX5.fna",
    "SLG1215_B10_10473203_HT53MAFX5.fna",
    "SLG441_A4_10473203_HT53MAFX5.fna",
    "SLG616_E1_10473203_HT53MAFX5.fna",
    "SLG616_F1_10473203_HT53MAFX5.fna",
    "SLG888_A4_10473203_HT53MAFX5.fna",
    "SLG888_E3_10473203_HT53MAFX5.fna",
    "SLG942_C9_10473203_HT53MAFX5.fna",
]

READ_PATTERNS = [
    "{asm}_{end}.fastq.gz",
    "{asm}_{end}*.fastq.gz",
    "{asm}*{end}*.fastq.gz",
    "{asm}_{end}*.fq.gz",
    "{asm}*{end}*.fq.gz",
]


def subset_bmuris(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Subset dataframe to rows matching target classification."""
    s = df[CLASS_COL].astype(str).str.strip().str.casefold()
    target_cf = target.strip().casefold()
    mask = s == target_cf
    if not mask.any():
        mask = s.str.contains(target_cf, na=False)
    return df.loc[mask].copy()


def strip_fna(user_genome: str) -> str:
    """Remove .fna or .fna.gz suffix from genome name."""
    ug = str(user_genome).strip()
    for suffix in (".fna.gz", ".fna"):
        if ug.endswith(suffix):
            return ug[: -len(suffix)]
    return ug


def ensure_fna(user_genome: str) -> str:
    """Ensure genome name ends with .fna."""
    ug = str(user_genome).strip()
    if ug.endswith(".fna.gz"):
        ug = ug[:-3]
    if not ug or ug.lower() == "nan" or ug.endswith(".fna"):
        return ug
    return ug + ".fna"


@dataclass
class ReadPair:
    """Paired-end read files for an assembly."""

    assembly: str
    r1: Path
    r2: Path


def find_reads_for_assembly(
    assembly: str, search_dirs: Sequence[Path]
) -> Optional[ReadPair]:
    """Find R1/R2 reads whose filenames start with assembly name."""
    for dir in search_dirs:
        r1_hits = _glob_reads(assembly, "R1", dir)
        r2_hits = _glob_reads(assembly, "R2", dir)
        if len(r1_hits) == 1 and len(r2_hits) == 1:
            return ReadPair(assembly, r1_hits[0], r2_hits[0])
        if r1_hits and r2_hits:
            log.warning(
                f"Multiple read files found for {assembly} in {dir}, "
                f"R1 candidates: {[p.name for p in r1_hits]}, "
                f"R2 candidates: {[p.name for p in r2_hits]}, "
            )
            for h1 in r1_hits:
                for h2 in r2_hits:
                    if h1.parent == h2.parent:
                        log.warning(f"matched pair: R1={h1.name}, R2={h2.name}")
                        return ReadPair(assembly, h1, h2)
    return None


def _glob_reads(assembly: str, end: str, base: Path) -> list[Path]:
    """Glob for read files matching assembly and end (R1/R2)."""
    hits: list[Path] = []
    for pattern in READ_PATTERNS:
        hits.extend(base.glob(pattern.format(asm=assembly, end=end)))
    return sorted(set(hits))


def copy_and_verify(src: Path, dst: Path) -> tuple[bool, str]:
    """Copy file from src to dst and verify by size."""
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists():
        if dst.stat().st_size == src.stat().st_size:
            return True, "exists_ok_size"
        dst.unlink()

    shutil.copy2(src, dst)

    if not dst.exists():
        return False, "copy_failed_missing_dst"
    if dst.stat().st_size != src.stat().st_size:
        return (
            False,
            f"size_mismatch src={src.stat().st_size} dst={dst.stat().st_size}",
        )
    return True, "copied_ok"


def verify_subset_against_expected(
    sub: pd.DataFrame,
    expected_list: list[str],
) -> tuple[bool, list[str], list[str]]:
    """Check if subset matches expected genome list."""
    expected_set = {ensure_fna(x) for x in expected_list if x and x.lower() != "nan"}
    subset_set = {
        ensure_fna(x)
        for x in sub[GENOME_COL].astype(str)
        if x and str(x).lower() != "nan"
    } - {""}
    missing = sorted(expected_set - subset_set)
    extra = sorted(subset_set - expected_set)
    return (not missing and not extra), missing, extra


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Subset Bacteroides muris isolates and copy reads.",
    )
    ap.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to isolate library table (.tsv).",
    )
    ap.add_argument(
        "--target",
        default="Bacteroides muris",
        help="Target classification string.",
    )
    ap.add_argument(
        "--output-table",
        required=True,
        type=Path,
        help="Where to write the subset table (.tsv).",
    )
    ap.add_argument(
        "--dest-reads-dir",
        type=Path,
        default=Path("/workdir1/sidd/for_sam/isolate_Bacteroides_muris/reads"),
    )
    ap.add_argument(
        "--batch1-dir",
        type=Path,
        default=Path(
            "/local/workdir/lab_data/Isolate_Sequencing/"
            "10471829_BRC_HackFlex_DietManipulation-Batch1-Isolates"
        ),
    )
    ap.add_argument(
        "--batch2-dir",
        type=Path,
        default=Path(
            "/local/workdir/lab_data/Isolate_Sequencing/"
            "10473203_BRC_HackFlex_DietManipulation-Batch2-Isolates"
        ),
    )
    return ap


def main() -> int:
    """Main entry point."""
    logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.INFO)
    args = build_parser().parse_args()

    # Read & subset
    df = pd.read_csv(args.input, sep="\t")
    log.info(f"Read {len(df)} rows from {args.input}")

    sub = subset_bmuris(df, args.target)
    log.info(f"Subset: {len(sub)} rows matching '{args.target}'")

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(args.output_table, sep="\t", index=False)
    log.info(f"Wrote subset table -> {args.output_table}")

    # Derive assemblies
    assemblies = sorted(
        {
            strip_fna(ug)
            for ug in sub[GENOME_COL].astype(str)
            if ug.strip() and ug.strip().lower() != "nan"
        }
    )
    log.info(f"Unique assemblies to fetch reads for: {len(assemblies)}")

    # Locate read pairs
    search_dirs = [args.batch2_dir, args.batch1_dir]
    read_pairs: dict[str, ReadPair] = {}
    missing_reads: list[str] = []

    for assembly in assemblies:
        rp = find_reads_for_assembly(assembly, search_dirs)
        if rp:
            read_pairs[assembly] = rp
        else:
            missing_reads.append(assembly)

    if missing_reads:
        log.warning(f"Missing reads for {len(missing_reads)} assemblies:")
        for assembly in missing_reads:
            log.warning(f"  - {assembly}")
        # return 3

    # Copy + verify
    any_copy_fail = False

    for assembly, rp in read_pairs.items():
        dst_r1 = args.dest_reads_dir / rp.r1.name
        dst_r2 = args.dest_reads_dir / rp.r2.name

        ok1, msg1 = copy_and_verify(rp.r1, dst_r1)
        ok2, msg2 = copy_and_verify(rp.r2, dst_r2)

        if not (ok1 and ok2):
            log.error(f"Copy failed for {assembly}: R1={msg1}, R2={msg2}")
            any_copy_fail = True

    if any_copy_fail:
        log.error("One or more read files failed copy/verification.")
        return 4

    log.info(f"All reads copied and verified (pairs: {len(read_pairs)})")

    # Expected-list check
    ok, missing, extra = verify_subset_against_expected(sub, EXPECTED_USER_GENOMES)
    if ok:
        log.info("Subset User.Genome list matches expected list")
        return 0

    log.warning("Subset User.Genome list does NOT match expected list.")
    if missing:
        log.warning(f"Missing from subset ({len(missing)}):")
        for x in missing:
            log.warning(f"  - {x}")
    if extra:
        log.warning(f"Extra in subset ({len(extra)}):")
        for x in extra:
            log.warning(f"  - {x}")

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
