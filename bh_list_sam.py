#!/usr/bin/env python3
"""
Compute BH‐corrected minima and counts for each test‐type across all files in a given period,
but only for basenames present in the 'pre_end' paired_tTest set. Produces two wide summary
TSVs (pre_end and pre_post) AND two detailed TSVs listing every significant row
(basename + contig + gene_id + position) for each period.
"""

import argparse
import csv
import os
import re

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def row_minima(path, col_regex):
    df = pd.read_csv(path, sep="\t", compression="gzip")
    cols = [c for c in df.columns if re.search(col_regex, c)]
    return df[cols].apply(pd.to_numeric, errors="coerce").min(axis=1, skipna=True)


def summarize_period(root_dir, period, basenames, out_path):
    # Prepare a detailed‐rows output alongside the summary; derive filename:
    detail_path = out_path.replace(".tsv", "_significant_rows.tsv")
    detail_f = open(detail_path, "w", newline="")
    writer = csv.writer(detail_f, delimiter="\t")
    writer.writerow(["period", "test", "basename", "contig", "gene_id", "position"])

    # Determine suffix for lmm_pre directories
    suffix_period = period.split("_", 1)[1]

    tests = [
        (
            "paired_tTest",
            f"two_sample_paired_{period}-fat_control",
            "_two_sample_paired.tsv.gz",
            r"_p_value_tTest",
        ),
        (
            "paired_Wilcoxon",
            f"two_sample_paired_{period}-fat_control",
            "_two_sample_paired.tsv.gz",
            r"_p_value_Wilcoxon",
        ),
        (
            "control_tTest",
            f"single_sample_{period}-fat_control",
            "_single_sample_control.tsv.gz",
            r"_p_value_tTest_control",
        ),
        (
            "control_Wilcoxon",
            f"single_sample_{period}-fat_control",
            "_single_sample_control.tsv.gz",
            r"_p_value_Wilcoxon_control",
        ),
        (
            "fat_tTest",
            f"single_sample_{period}-fat_control",
            "_single_sample_fat.tsv.gz",
            r"_p_value_tTest_fat",
        ),
        (
            "fat_Wilcoxon",
            f"single_sample_{period}-fat_control",
            "_single_sample_fat.tsv.gz",
            r"_p_value_Wilcoxon_fat",
        ),
        (
            "lmm_pre",
            f"lmm_pre_{suffix_period}-fat_control",
            "_lmm.tsv.gz",
            r"_p_value_LMM",
        ),
        (
            "lmm_across_time_control",
            f"lmm_across_time_{period}-fat_control",
            "_lmm_across_time_control.tsv.gz",
            r"_p_value_LMM",
        ),
        (
            "lmm_across_time_fat",
            f"lmm_across_time_{period}-fat_control",
            "_lmm_across_time_fat.tsv.gz",
            r"_p_value_LMM",
        ),
    ]

    # store per-test -> {basename: (min_q, n_sig)}
    all_results = {name: {} for name, *_ in tests}

    for name, subdir, suffix, col_regex in tests:
        dirpath = os.path.join(root_dir, subdir)
        if not os.path.isdir(dirpath):
            continue

        segments, lengths, bases, file_paths = [], [], [], []
        for fn in sorted(os.listdir(dirpath)):
            if not fn.endswith(suffix):
                continue
            base = fn[: -len(suffix)]
            if base not in basenames:
                continue
            path = os.path.join(dirpath, fn)
            mins = row_minima(path, col_regex)
            segments.append(mins)
            lengths.append(len(mins))
            bases.append(base)
            file_paths.append(path)

        if not segments:
            continue

        concat = pd.concat(segments, ignore_index=True)
        mask = ~concat.isna()
        vals = concat[mask].to_numpy()
        _, qvals, _, _ = multipletests(vals, method="fdr_bh")

        full_q = np.full(len(concat), np.nan)
        full_q[mask.to_numpy()] = qvals

        idx = 0
        for base, ln, path in zip(bases, lengths, file_paths):
            seg = full_q[idx : idx + ln]
            idx += ln

            # record detailed rows for q < 0.05
            if ln > 0:
                df = pd.read_csv(path, sep="\t", compression="gzip")
                for i, q in enumerate(seg):
                    if q < 0.05:
                        writer.writerow(
                            [
                                period,
                                name,
                                base,
                                df.at[i, "contig"],
                                df.at[i, "gene_id"],
                                df.at[i, "position"],
                            ]
                        )

            # record summary
            if ln > 0:
                all_results[name][base] = (
                    float(np.nanmin(seg)),
                    int(np.nansum(seg < 0.05)),
                )
            else:
                all_results[name][base] = (np.nan, np.nan)

    detail_f.close()  # finish writing the detailed rows file

    # Build and write the wide summary DataFrame
    rows = []
    for base in sorted(basenames):
        row = {"basename": base}
        for name in all_results:
            min_q, n_sig = all_results[name].get(base, (np.nan, np.nan))
            row[f"{name}_min_q"] = min_q
            row[f"{name}_n_sig"] = n_sig
        rows.append(row)

    df = pd.DataFrame(rows)
    df.fillna("NA", inplace=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {period} summary to {out_path}")
    print(f"Wrote {period} significant‐rows to {detail_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--root_dir", required=True, help="Path to significance_tests/ folder"
    )
    p.add_argument("--pre_end_out", required=True, help="TSV for pre_end summary")
    p.add_argument("--pre_post_out", required=True, help="TSV for pre_post summary")
    args = p.parse_args()

    pe_dir = os.path.join(args.root_dir, "two_sample_paired_pre_end-fat_control")
    basenames = {
        fn[: -len("_two_sample_paired.tsv.gz")]
        for fn in os.listdir(pe_dir)
        if fn.endswith("_two_sample_paired.tsv.gz")
    }

    summarize_period(args.root_dir, "pre_end", basenames, args.pre_end_out)
    summarize_period(args.root_dir, "pre_post", basenames, args.pre_post_out)


if __name__ == "__main__":
    main()
