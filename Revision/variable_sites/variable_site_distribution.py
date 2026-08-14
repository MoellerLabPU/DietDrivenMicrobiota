#!/usr/bin/env python3
"""Variable-site counts and spacing for the 62 divergence-tested Fig 1 SGBs, END vs PRE.

A **variable site is a site that was TESTED** -- one that appears in p_value_summary
for a given test family. The between-group and within-group tests preprocess
separately, so their tested sets differ (21.5% of HF-tested sites are not
divergence-tested). Q1/Q3/Q4 are therefore reported under four definitions:

    div    tested by two_sample_paired_tTest        (the divergence test)
    hf     tested by single_sample_tTest / fat      (HF parallelism)
    lf     tested by single_sample_tTest / control  (LF parallelism)
    union  tested by any of the above

Each significance question is judged against its OWN tested set, so a hit can never
fall outside the denominator it is expressed against: Q5/Q6 vs div, Q7/Q8 vs hf,
Q9/Q10 vs lf.

Answers, per contig and per SGB:
  1  contigs containing >=1 variable site       (x4 definitions)
  2  contigs in the SGB reference
  3  mean distance (nt) between variable sites  (x4 definitions)
  4  mean distance between randomly drawn sites (x4 definitions)
  5/6   contigs with, and spacing of, divergence-significant sites
  7/8   same for HF parallelism
  9/10  same for LF parallelism

Distances are computed WITHIN a contig only -- a gap spanning two contigs is
undefined -- then pooled across the SGB's contigs.

Outputs four things: one table per MAG (a row per contig), a combined contig-level
table, an SGB-level table, and a summary. Read-only w.r.t. the AlleleFlux run.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

def setup_logging(level=logging.INFO) -> None:
    """Configure logging once, in main(). Same format as the AlleleFlux scripts."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


logger = logging.getLogger(__name__)

# Pinned in ONE place on purpose. Several AlleleFlux run trees hold files with
# identical names and different contents -- scores_two_sample_paired-pre_end-...
# has 72 / 67 / 62 / 22 rows in AlleleFlux_final / _mapq2 / _mapq20 / revision.
# Published Fig 1 is mapq20 (62 SGBs). Never glob across AlleleFlux_*, and note
# /scratch/.../sidd/diet_manip/AlleleFlux_mapq20 is a DIFFERENT tree from this one.
RUN = Path("/scratch/gpfs/AMOELLER/diet_manip/AlleleFlux_mapq20/longitudinal")
# Reference geometry. The fasta path written in figures/Fig1/alleleflux_config_mapq20.yml
# is broken (missing a sg4230/ component); these two are the working copies.
FAI = Path("/scratch/gpfs/AMOELLER/diet_manip/redo_representative_megamag.fa.fai")
MAPPING = Path("/scratch/gpfs/AMOELLER/diet_manip/megamag_mapping.tsv")
# "END vs PRE" in the ask. pre_post also exists in this run and is excluded.
COMPARISON = "pre_end-fat_control"

# The SGB set. Verified: the MAGs tested by two_sample_paired_tTest are EXACTLY the
# 62 rows of this table, so "Fig 1 SGBs" and "MAGs tested in the divergence test"
# name the same set. HF is tested in 72 MAGs and LF in 68; both are cut to these 62.
FIG1_SCORES = (
    RUN / "scores/processed/combined/MAG"
    / "scores_two_sample_paired-pre_end-fat_control-MAGs.tsv"
)

# Divergence = between-group. Four between-group tests were run; the paired t-test
# is the project convention AND the only one with a clean permutation null (57/62
# MAGs significant in real data, 0/69 under all three label permutations, whereas
# CMH and LMM give MORE hits when permuted). Wilcoxon cannot reach q<0.05 at n=8.
DIVERGENCE_TEST = "two_sample_paired_tTest"
# Parallelism = within-group: does one diet group move consistently over time?
PARALLELISM_TEST = "single_sample_tTest"
# HF/LF -> fat/control is an ASSUMPTION. The Fig 1 code only ever says Fat/Control;
# the sole HF/LF mapping in the repo is Revision/.../per_litter_boxplots.ipynb.
HF_GROUP = "fat"
LF_GROUP = "control"
# BH-corrected. Note the correction is applied genome-wide across all MAGs pooled,
# not per SGB, so "significant for this SGB" means "significant in a global ranking".
Q_THRESHOLD = 0.05

# The three tested sets, in report order; "union" is derived, not read from disk.
TEST_SETS = ("div", "hf", "lf")
DEFINITIONS = TEST_SETS + ("union",)
# Only these three carry a significance question (Q5/6, Q7/8, Q9/10).
SIG_SETS = TEST_SETS

# The p_value_summary files are 386 MB and 356 MB of uncompressed TSV. Reading one
# whole with default dtypes balloons to several GB, so stream it.
CHUNK = 2_000_000

# Only the species label is carried through. `domain` is constant across all 62
# SGBs and class/order are near-duplicates of phylum (7/7/9 distinct values), so the
# other ranks add columns without adding information. Add "phylum" here if a figure
# ever needs to group or colour by clade.
TAXONOMY_COLS = ["species"]

# Returned for a contig a definition never tested. Typed int64 on purpose: a bare
# np.array([]) is float64 and would quietly change the dtype of position arrays.
EMPTY = np.empty(0, dtype=np.int64)


def load_reference_geometry() -> tuple[pd.Series, pd.Series]:
    """contig -> length, and mag_id -> set of its reference contigs.

    The .fai is a plain 2-column table (name, length), so pandas reads it directly.
    Biopython/pysam would also parse it, but this analysis never touches sequence --
    only contig lengths -- so neither would buy anything but a dependency.
    """
    fai = pd.read_csv(
        FAI, sep="\t", header=None, usecols=[0, 1], names=["contig", "length"]
    )
    contig_len = fai.set_index("contig")["length"]
    mapping = pd.read_csv(MAPPING, sep="\t")
    # A hard fail, not a warning: if the mapping and the reference disagree then every
    # Q2 count and every random-gap expectation is silently wrong. The join is exactly
    # 1:1 for this pair of files, and the summed lengths reproduce the genome_size
    # column AlleleFlux writes in its own QC tables for all 160 MAGs.
    missing = ~mapping["contig_id"].isin(contig_len.index)
    if missing.any():
        logger.error(
            f"{int(missing.sum())} contigs in the mapping have no length in the .fai "
            "— the mapping and reference do not correspond."
        )
        raise SystemExit(1)
    return contig_len, mapping.groupby("mag_id")["contig_id"].apply(set)


def load_fig1_mags() -> tuple[list[str], pd.DataFrame]:
    """The Fig-1 SGB list, plus whatever GTDB taxonomy the scores table carries."""
    df = pd.read_csv(FIG1_SCORES, sep="\t")
    mags = df["MAG_ID"].dropna().unique().tolist()
    # Inferred, not stated by the figure code: scores_figs.Rmd gates every panel on a
    # variable `two_samp_pre_end` that is read 19x and never assigned anywhere in the
    # repo. This 62-MAG list matches paired_test_eligible==TRUE in the eligibility
    # table exactly, which is the strongest available evidence for what Fig 1 showed.
    have = [c for c in TAXONOMY_COLS if c in df.columns]
    taxonomy = df[["MAG_ID"] + have].drop_duplicates("MAG_ID")
    logger.info(f"Fig-1 SGB set: {len(mags)} MAGs from {FIG1_SCORES.name}")
    return mags, taxonomy


def load_tested_sites(
    path: Path, test_type: str, group: str | None, mags: set[str]
) -> pd.DataFrame:
    """Every site this test EVALUATED for the given MAGs, with its q_value attached.

    Returns tested sites, not significant ones. The q_value rides along so the
    significant subset is a filter on the result rather than a second pass over
    hundreds of megabytes.

    The file holds one row per (site x test_type [x group_analyzed]); pinning BOTH is
    what makes a row a site, and uniqueness is then asserted rather than imposed. Do
    NOT "fix" a duplicate count with drop_duplicates(): forgetting group_analyzed
    collides the two diet groups on the same position (721,676 rows for 517,704 sites
    over the unthresholded single_sample_tTest rows), and deduping would silently
    merge fat and control into one plausible-looking but wrong set.
    """
    # usecols + explicit dtypes keep this to ~100 MB instead of several GB; the files
    # also carry gene_id and source_file, which we never need.
    label = f"{test_type}{'/' + group if group else ''}"
    usecols = ["mag_id", "contig", "position", "test_type", "q_value"]
    if group is not None:
        usecols.append("group_analyzed")

    kept = []
    total = 0
    # Stream in chunks: these files are 386 MB and 356 MB, and reading one whole with
    # default dtypes balloons to several GB.
    for chunk in pd.read_csv(
        path, sep="\t", usecols=usecols, chunksize=CHUNK,
        dtype={"mag_id": "string", "contig": "string", "position": "int64",
               "test_type": "string", "q_value": "float64"},
    ):
        total += len(chunk)
        sel = (chunk["test_type"] == test_type) & chunk["mag_id"].isin(mags)
        if group is not None:
            sel &= chunk["group_analyzed"] == group
        if sel.any():
            kept.append(chunk.loc[sel, ["mag_id", "contig", "position", "q_value"]])

    if not kept:
        logger.error(f"{path.name} [{test_type}]: no tested sites at all — wrong file?")
        raise SystemExit(1)

    out = pd.concat(kept, ignore_index=True)
    # A site is (mag_id, contig, position). position is 0-based and position 0 really
    # occurs in this data, so never treat it as missing or falsy.
    key = ["mag_id", "contig", "position"]
    dupes = len(out) - len(out.drop_duplicates(key))
    if dupes:
        logger.error(
            f"{path.name} [{label}]: {dupes:,} duplicate rows on {key} after filtering. "
            "The filters no longer pin one row per site — check the file's schema "
            "before trusting any count."
        )
        raise SystemExit(1)
    # Categoricals: 62 MAGs and ~8k contigs across ~1M rows, so a large saving.
    out["mag_id"] = out["mag_id"].astype("category")
    out["contig"] = out["contig"].astype("category")
    n_sig = int((out["q_value"] < Q_THRESHOLD).sum())
    logger.info(
        f"{path.name} [{label}]: {len(out):,} tested sites, {n_sig:,} significant "
        f"(q<{Q_THRESHOLD:g}), {out['mag_id'].nunique()} MAGs "
        f"(scanned {total:,} rows)"
    )
    return out


def gap_stats(positions: np.ndarray) -> dict:
    """Spacing within ONE contig. n sites give n-1 gaps, so <2 sites give none."""
    if len(positions) < 2:
        # NaN rather than 0: "no gap could be measured" is not "the gap is zero".
        return {"mean_gap": np.nan, "median_gap": np.nan, "n_gaps": 0}
    # sort first: positions arrive grouped by contig, not in coordinate order. The
    # gaps telescope to (max - min), so the mean is span/(n-1) whatever the pattern.
    g = np.diff(np.sort(positions))
    return {
        "mean_gap": float(g.mean()),
        "median_gap": float(np.median(g)),
        "n_gaps": int(g.size),
    }


def expected_gap(contig_length: float, n: int) -> float:
    """Uniform-null spacing on one contig: (L+1)/(n+1).

    For n sites placed uniformly among L positions the expected distance between
    neighbours is (L+1)/(n+1). E[first site] simplifies from C(L+1,n+1)/C(L,n) to
    (L+1)/(n+1); mirroring gives E[last], so E[span] = (L+1)(n-1)/(n+1), and dividing
    by the n-1 gaps cancels the (n-1). Checked against exhaustive enumeration.

    Deliberately NOT L/n: sites never reach the contig ends, so the span they cover
    is under L and the gaps are tighter than the density implies. The two agree to
    0.002% at n ~ 15k but differ ~25% on a short contig.
    """
    # The `contig_length == contig_length` test is a NaN check: a contig missing from
    # the .fai has no length, so no expectation can be formed.
    if n < 2 or not (contig_length == contig_length):
        return np.nan
    return (contig_length + 1) / (n + 1)


def _by_contig(df: pd.DataFrame) -> dict[str, np.ndarray]:
    """contig -> positions. Built once per MAG; filtering per contig would be O(n*m)."""
    if not len(df):
        return {}
    return {c: g["position"].to_numpy()
            for c, g in df.groupby("contig", observed=True)}


def contig_rows(mag: str, tested: dict[str, pd.DataFrame], contig_len: pd.Series):
    """One row per contig for this MAG, covering all four definitions."""
    # Two lookup tables per definition: all tested positions, and the significant
    # subset. Both built once here rather than re-filtered inside the contig loop.
    pos = {d: _by_contig(tested[d]) for d in DEFINITIONS}
    sig = {d: _by_contig(tested[d][tested[d]["q_value"] < Q_THRESHOLD])
           for d in SIG_SETS}

    # Every contig touched by ANY definition gets a row, so a per-MAG table is a
    # complete picture rather than only the divergence-tested contigs.
    contigs = sorted(set().union(*(set(pos[d]) for d in DEFINITIONS)))

    rows = []
    for c in contigs:
        L = float(contig_len.get(c, np.nan))
        # -1 flags a contig with no length in the .fai. It should never happen (the
        # mapping/fai join is checked at startup) but a sentinel beats a silent NaN.
        row = {"MAG_ID": mag, "contig": c,
               "contig_len": int(L) if L == L else -1}
        # One block per definition, so every div_* / hf_* / lf_* / union_* column
        # sits together in the output rather than being split across the table.
        for d in DEFINITIONS:
            # Tested sites -> Q1/Q3/Q4. A contig may be tested by one definition and
            # not another, hence .get() with a typed empty array rather than [c].
            p = pos[d].get(c, EMPTY)
            g = gap_stats(p)
            row[f"{d}_n_sites"] = int(len(p))          # Q1
            row[f"{d}_mean_gap"] = g["mean_gap"]       # Q3
            row[f"{d}_median_gap"] = g["median_gap"]   # Q3
            row[f"{d}_n_gaps"] = g["n_gaps"]           # Q3 support count
            row[f"{d}_expected_gap"] = expected_gap(L, len(p))  # Q4

            # Significant subset -> Q5-Q10. Skipped for "union", which is a union of
            # tested sets and has no test, hence no significance question, of its own.
            if d not in SIG_SETS:
                continue
            p = sig[d].get(c, EMPTY)
            g = gap_stats(p)
            row[f"{d}_n_sig"] = int(len(p))
            row[f"{d}_has_sig"] = bool(len(p))         # Q5/Q7/Q9 count contigs by this
            row[f"{d}_sig_mean_gap"] = g["mean_gap"]   # Q6/Q8/Q10
            row[f"{d}_sig_median_gap"] = g["median_gap"]
            row[f"{d}_sig_n_gaps"] = g["n_gaps"]
        rows.append(row)
    return pd.DataFrame(rows)


def _pool(means: pd.Series, counts: pd.Series) -> float:
    """Gap-count-weighted mean, reconstructed from per-contig means.

    Pooling weights each contig by how many gaps it contributes, so a dense contig
    dominates a sparse one; averaging per-contig means instead would weight them
    equally, which is a different statistic and not a rounding difference. Deriving
    it from the contig table guarantees the two levels can never disagree.
    """
    # Rows with no gaps carry a NaN mean and would poison np.average, so drop them
    # rather than zero-filling: with a weight of 0 they contribute nothing either
    # way, and "count > 0 implies the mean exists" is an invariant of gap_stats.
    ok = counts > 0
    return float(np.average(means[ok], weights=counts[ok])) if ok.any() else np.nan


def sgb_row(mag: str, cr: pd.DataFrame, n_ref: int, ref_len: int) -> dict:
    """Aggregate one MAG's contig table into a single SGB-level row."""
    row = {"MAG_ID": mag, "n_ref_contigs": n_ref, "ref_genome_len": ref_len}
    # Q1-Q4 rolled up. Every figure is rebuilt from the contig table rather than
    # recomputed from the sites, so the two levels cannot drift apart.
    for d in DEFINITIONS:
        n_sites = int(cr[f"{d}_n_sites"].sum())
        row[f"{d}_n_contigs_with_sites"] = int((cr[f"{d}_n_sites"] > 0).sum())
        row[f"{d}_n_sites"] = n_sites
        row[f"{d}_mean_gap"] = _pool(cr[f"{d}_mean_gap"], cr[f"{d}_n_gaps"])
        row[f"{d}_expected_gap"] = _pool(cr[f"{d}_expected_gap"], cr[f"{d}_n_gaps"])
        row[f"{d}_n_gaps"] = int(cr[f"{d}_n_gaps"].sum())
        row[f"{d}_pct_genome_tested"] = 100 * n_sites / ref_len if ref_len else np.nan
    for d in SIG_SETS:
        # Denominator is this test's OWN tested contigs, so a significant site can
        # never fall outside the set its fraction is expressed against.
        denom = int((cr[f"{d}_n_sites"] > 0).sum())
        nsig = int(cr[f"{d}_has_sig"].sum())
        row[f"{d}_n_contigs_sig"] = nsig
        row[f"{d}_pct_contigs_sig"] = 100 * nsig / denom if denom else np.nan
        row[f"{d}_n_sig_sites"] = int(cr[f"{d}_n_sig"].sum())
        row[f"{d}_sig_mean_gap"] = _pool(cr[f"{d}_sig_mean_gap"], cr[f"{d}_sig_n_gaps"])
        row[f"{d}_sig_n_gaps"] = int(cr[f"{d}_sig_n_gaps"].sum())
    return row


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Per-contig and per-SGB variable-site counts and spacing for the 62 "
            "divergence-tested Fig 1 SGBs, END vs PRE. A variable site is a site "
            "that was TESTED; see METHODS.md for the four definitions. Inputs are "
            "pinned as module constants, not flags: several AlleleFlux run trees "
            "hold identically named files with different contents, so pointing this "
            "at another run should be a deliberate edit."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Directory to write contigs/<MAG>.tsv, contig_level_all.tsv, "
             "sgb_level.tsv and summary.tsv into.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N SGBs (smoke test).",
    )
    args = parser.parse_args()

    setup_logging()
    (args.outdir / "contigs").mkdir(parents=True, exist_ok=True)

    contig_len, ref_contigs = load_reference_geometry()
    mags, taxonomy = load_fig1_mags()
    if args.limit:
        mags = mags[: args.limit]
    magset = set(mags)

    # Three passes over ~740 MB, done ONCE up front rather than per MAG. Tested sites
    # arrive with their q_value, so significance costs no extra pass. Note this no
    # longer touches allele_analysis at all -- "tested" lives in p_value_summary.
    logger.info("Loading tested sites (3 streaming passes)...")
    pv = RUN / "p_value_summary" / COMPARISON
    period = COMPARISON.split("-", 1)[0]
    paired = pv / f"p_value_summary_two_sample_paired_{period}.tsv"
    single = pv / f"p_value_summary_single_sample_{period}.tsv"
    # One flat table per definition. hf and lf read the SAME file and differ only in
    # group_analyzed.
    tested_tables = {
        "div": load_tested_sites(paired, DIVERGENCE_TEST, None, magset),
        "hf": load_tested_sites(single, PARALLELISM_TEST, HF_GROUP, magset),
        "lf": load_tested_sites(single, PARALLELISM_TEST, LF_GROUP, magset),
    }
    # Re-key each table by MAG so the loop below is a dict lookup. Filtering a
    # ~400k-row table inside the loop instead would mean 62 full scans per definition.
    # observed=True because mag_id is categorical: without it pandas emits an empty
    # group for every category value, including MAGs absent from this frame.
    by_mag = {
        d: dict(tuple(df.groupby("mag_id", observed=True)))
        for d, df in tested_tables.items()
    }

    # Stand-in for a MAG a given test never evaluated, so the loop needs no special
    # case; contig_rows() then yields zero rows for that definition.
    empty = pd.DataFrame({"mag_id": [], "contig": [], "position": [], "q_value": []})
    contig_tables, sgb_rows = [], []
    for i, mag in enumerate(mags, 1):
        tested = {d: by_mag[d].get(mag, empty) for d in TEST_SETS}
        # "union" = tested by anything. A site tested by two families appears in both
        # frames, so dedupe on (contig, position); keep the smallest q_value so the
        # union's own significance count stays meaningful.
        u = pd.concat([tested[d] for d in TEST_SETS], ignore_index=True)
        tested["union"] = (
            u.sort_values("q_value").drop_duplicates(["contig", "position"])
            if len(u) else u
        )

        cr = contig_rows(mag, tested, contig_len)
        if cr.empty:
            logger.warning(f"[{i}/{len(mags)}] {mag}: no tested sites — skipped")
            continue
        cr.to_csv(args.outdir / "contigs" / f"{mag}.tsv", sep="\t", index=False,
                  float_format="%.2f")
        contig_tables.append(cr)

        # Q2 comes from the reference alone, so it counts every contig the SGB has,
        # including ones no test ever reached.
        ref = ref_contigs.get(mag, set())
        sgb_rows.append(
            sgb_row(mag, cr, len(ref), int(sum(contig_len.get(c, 0) for c in ref)))
        )
        r = sgb_rows[-1]
        logger.info(
            f"[{i}/{len(mags)}] {mag}: {len(cr)} contigs | "
            f"div {r['div_n_sites']:,} sites on {r['div_n_contigs_with_sites']} contigs "
            f"| union {r['union_n_sites']:,}"
        )

    # ---- combined contig-level table ----
    allc = pd.concat(contig_tables, ignore_index=True)
    allc.to_csv(args.outdir / "contig_level_all.tsv", sep="\t", index=False,
                float_format="%.2f")

    # ---- SGB level ----
    sgb = pd.DataFrame(sgb_rows).merge(taxonomy, on="MAG_ID", how="left")
    # Identify-the-row columns first, metrics after.
    lead = ["MAG_ID"] + [c for c in TAXONOMY_COLS if c in sgb.columns]
    sgb = sgb[lead + [c for c in sgb.columns if c not in lead]]
    sgb.to_csv(args.outdir / "sgb_level.tsv", sep="\t", index=False, float_format="%.2f")

    # ---- summary: one row per definition, aggregated over all 62 SGBs ----
    # Pooled over SGBs the same way the SGB level pools over contigs: weighted by
    # gap count, never a mean of means.
    srows = []
    for d in DEFINITIONS:
        r = {
            "definition": d,
            "n_SGBs": int((sgb[f"{d}_n_sites"] > 0).sum()),
            "n_ref_contigs": int(sgb["n_ref_contigs"].sum()),
            "n_contigs_with_sites": int(sgb[f"{d}_n_contigs_with_sites"].sum()),
            "n_sites": int(sgb[f"{d}_n_sites"].sum()),
            "mean_gap": _pool(sgb[f"{d}_mean_gap"], sgb[f"{d}_n_gaps"]),
            "expected_gap": _pool(sgb[f"{d}_expected_gap"], sgb[f"{d}_n_gaps"]),
            "n_gaps": int(sgb[f"{d}_n_gaps"].sum()),
        }
        # "union" has no significance question of its own, so these stay absent
        # rather than being filled with a meaningless value.
        if d in SIG_SETS:
            r["n_contigs_sig"] = int(sgb[f"{d}_n_contigs_sig"].sum())
            r["pct_contigs_sig"] = (
                100 * r["n_contigs_sig"] / r["n_contigs_with_sites"]
                if r["n_contigs_with_sites"] else np.nan
            )
            r["n_sig_sites"] = int(sgb[f"{d}_n_sig_sites"].sum())
            r["sig_mean_gap"] = _pool(sgb[f"{d}_sig_mean_gap"], sgb[f"{d}_sig_n_gaps"])
            r["sig_n_gaps"] = int(sgb[f"{d}_sig_n_gaps"].sum())
        srows.append(r)
    summary = pd.DataFrame(srows)
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False, float_format="%.2f")

    logger.info(f"Wrote {len(contig_tables)} per-MAG tables to {args.outdir / 'contigs'}")
    logger.info(f"Wrote contig_level_all.tsv  {allc.shape[0]:,} rows x {allc.shape[1]} cols")
    logger.info(f"Wrote sgb_level.tsv         {sgb.shape[0]} rows x {sgb.shape[1]} cols")
    logger.info(f"Wrote summary.tsv           {summary.shape[0]} rows x {summary.shape[1]} cols")
    return 0


if __name__ == "__main__":
    sys.exit(main())
