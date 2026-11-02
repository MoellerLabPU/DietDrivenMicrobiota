#!/usr/bin/env python3
"""Pivot one MAG's AlleleFlux per-mouse PRE/END allele frequencies into a wide table.

Input: the per-subject changes table AlleleFlux writes before averaging per
replicate (the averaged version is what the preprocessing step reads):
    <run_root>/allele_analysis/allele_analysis_<comparison>/<MAG>_allele_frequency_changes.tsv.gz
    long format, one row per site x mouse, columns:
        subjectID gene_id contig position replicate group
        total_coverage_pre total_coverage_end total_coverage_combined
        A/T/G/C_frequency_pre  A/T/G/C_frequency_end  A/T/G/C_frequency_diff
    A row exists only where the mouse had coverage at BOTH timepoints.

Output (two files in --out_dir):
    <MAG>_allele_frequency_changes_wide.tsv
        rows    = every (contig, position, gene_id) in the input, nothing dropped
        columns = <group>_rep<replicate>_mouse<subjectID>_<value> for every
                  value column above, copied as-is; blank where that mouse has
                  no row at the site. Mice of the same replicate are adjacent
                  (rep1 control mice, rep1 fat mice, rep2 control mice, ...).
    <MAG>_sample_metadata.tsv
        sample_id subjectID time group replicate -- copied from the run's
        per-MAG roster <run_root>/inputMetadata/inputMetadata_<comparison>/<MAG>_metadata.tsv

Usage:
    python wide_allele_frequency_table.py --run_root <...>/longitudinal \
        --mag <MAG> --out_dir <dir> [--comparison pre_end-fat_control]
"""
import argparse
from pathlib import Path

import pandas as pd

SITE = ["contig", "position", "gene_id"]
KEY = ["group", "replicate", "subjectID"]


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--run_root", required=True, help="AlleleFlux run root, .../longitudinal")
    p.add_argument("--comparison", default="pre_end-fat_control")
    p.add_argument("--mag", required=True)
    p.add_argument("--out_dir", required=True)
    a = p.parse_args()

    root = Path(a.run_root)
    src = (root / "allele_analysis" / f"allele_analysis_{a.comparison}"
           / f"{a.mag}_allele_frequency_changes.tsv.gz")
    roster = root / "inputMetadata" / f"inputMetadata_{a.comparison}" / f"{a.mag}_metadata.tsv"
    for f in (src, roster):
        if not f.exists():
            raise FileNotFoundError(f)

    # Key columns as category: the table has tens of millions of rows.
    df = pd.read_csv(src, sep="\t",
                     dtype={"contig": "category", "gene_id": "category", "subjectID": "category",
                            "replicate": "category", "group": "category"})
    # gene_id carries a trailing space and is empty (read as NaN) at intergenic sites.
    df["gene_id"] = df["gene_id"].astype("object").fillna("").str.strip().astype("category")
    values = [c for c in df.columns if c not in SITE + KEY]

    # (site, group, replicate, subjectID) is unique in the input, so this is a pure reshape.
    wide = df.pivot(index=SITE, columns=KEY, values=values)
    # Column order: replicate -> group -> mouse -> value.
    keys = sorted(df[KEY].drop_duplicates().itertuples(index=False),
                  key=lambda k: (int(k.replicate), k.group, int(k.subjectID)))
    cols = [(v, k.group, k.replicate, k.subjectID) for k in keys for v in values]
    wide = wide.reindex(columns=pd.MultiIndex.from_tuples(cols))
    wide.columns = [f"{g}_rep{r}_mouse{s}_{v}" for v, g, r, s in wide.columns]
    # Blanks make pandas upcast the coverage counts to float; keep them integer.
    for c in [c for c in wide.columns if "total_coverage" in c]:
        wide[c] = wide[c].astype("Int64")
    wide = wide.sort_index().reset_index()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out / f"{a.mag}_allele_frequency_changes_wide.tsv", sep="\t", index=False)
    meta = pd.read_csv(roster, sep="\t", dtype=str)[["sample_id", "subjectID", "time", "group", "replicate"]]
    meta.to_csv(out / f"{a.mag}_sample_metadata.tsv", sep="\t", index=False)
    print(f"wrote {out / f'{a.mag}_allele_frequency_changes_wide.tsv'}: "
          f"{wide.shape[0]:,} sites x {len(keys)} mice x {len(values)} values")


if __name__ == "__main__":
    main()
