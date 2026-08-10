#!/usr/bin/env python3
"""
Per-cell relative abundance for the diet-manipulation AlleleFlux run
====================================================================

Builds ``rel_abundance_by_cell.tsv``: one row per significant-sites heatmap *cell*, carrying that
MAG's relative abundance at PRE and END and its change between them, aggregated over exactly the
samples the corresponding statistical test used.

See DESIGN.md for the full rationale. The three things that matter most:

1. **The unit of observation is the REPLICATE (cage), not the mouse.** AlleleFlux collapses every
   mouse sharing a ``replicate`` into a single unweighted mean before any test runs
   (``allele_freq.py:638-657``, ``get_mean_change``). ``single_sample`` therefore runs a one-sample
   t-test on 8 replicate values, and ``two_sample_paired`` pairs the two arms *by replicate*
   (``two_sample_paired.py:80-86``). Averaging over mice instead of cages would describe a different
   population than the p-value does -- cage 2 holds one mouse per arm while the rest hold two, so the
   two averages genuinely differ.

2. **Relative abundance is used as published.** Table S3 is already a per-sample composition summing
   to 100% over all 160 genomes; it is not renormalized to the subset AlleleFlux tested.

3. **Change is paired within each mouse.** ``delta = RA(end) - RA(pre)`` per subjectID, then averaged
   into cages, then over cages.

Aggregation, in three steps:

    step 1  per mouse       delta(mouse) = RA(end) - RA(pre)
    step 2  per (replicate, arm)  d(r) = unweighted mean over that arm's mice in cage r
    step 3  per cell        aggregate the replicate values the cell's test used

A *cell* is one row of ``significant_sites_mag_cell_stats_long.tsv``, keyed by
(comparison, test_type, group_analyzed, mag_id). One MAG appears in three of them and gets three
different abundance changes:

    single_sample_tTest / fat       -> the 8 fat-arm replicate values
    single_sample_tTest / control   -> the 8 control-arm replicate values
    two_sample_paired_tTest / ""    -> all 16 pooled, plus the within-replicate fat-minus-control
                                       contrast, which is what that test is actually sensitive to

Scope is PRE vs END only (``pre_end-fat_control``) and the two t-test families.

Error policy: every structural expectation is asserted, and a violation raises. Nothing is silently
dropped, repaired, or averaged away -- see ``_check`` and the assertions through ``main``.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

COMPARISON = "pre_end-fat_control"
T1, TN = "pre", "end"
# cell_stats `test_type` -> the `test_family` label carried in the output (and used to join).
TEST_FAMILY = {
    "single_sample_tTest": "single_sample",
    "two_sample_paired_tTest": "two_sample",
}
N_REPLICATES = (
    8  # cages; also == min_sample_num for this run, which is what makes step 3 exact
)
FLOAT_TOL = 1e-9

OUT_COLS = [
    "comparison",
    "period",
    "group_pair",
    "test_family",
    "group",  # single_sample: the analyzed arm ; two_sample: "" (pooled across both arms)
    "mag_id",
    "ra_pre",
    "ra_end",
    "ra_mean",
    "delta_ra",
    "abs_mean_delta_ra",  # |mean of the replicate deltas| -- size of the net shift
    "mean_abs_delta_ra",  # mean of |replicate delta|      -- size of the typical cage's shift
    "n_replicates",
    "n_mice",
    # two_sample rows only; blank on single_sample rows.
    "delta_ra_fat",
    "delta_ra_control",
    "delta_ra_contrast",
    "abs_mean_delta_ra_contrast",
    "mean_abs_delta_ra_contrast",
]


class DataIntegrityError(RuntimeError):
    """A structural expectation about the inputs was violated. Never caught -- always fatal."""


def _check(condition: bool, message: str) -> None:
    """Raise ``DataIntegrityError(message)`` unless ``condition`` holds.

    Used instead of bare ``assert`` so the checks survive ``python -O`` and carry a message aimed at
    whoever has to fix the inputs, not at whoever wrote this file.
    """
    # A plain function call, deliberately NOT `assert`: `python -O` strips assert statements from
    # the bytecode, which would silently disable every integrity check in this script.
    if not condition:
        # Raise (never log-and-continue): DataIntegrityError is uncaught by design, so a violated
        # input assumption kills the run instead of producing a plausible-looking wrong table.
        raise DataIntegrityError(message)


def load_table_s3(path: Path) -> pd.DataFrame:
    """Read Table S3 into a long frame of ``(mag_id, sample_id, ra)``.

    The file carries a one-line title above the real header, so row 0 is skipped and row 1 becomes
    the column names; the leading (unnamed) column holds the genome IDs. Trailing all-empty columns
    produced by the spreadsheet export are dropped by name.

    The 100%-per-sample check is the load-bearing one: it confirms the file really is the published
    composition (a share of the whole community) rather than raw counts, which is the assumption the
    whole "used as published, not renormalized" decision rests on.
    """
    # skiprows=1 jumps the one-line title so the real header becomes the columns; index_col=0 takes
    # the leading unnamed column as the genome IDs.
    wide = pd.read_csv(path, sep="\t", skiprows=1, index_col=0)
    # Spreadsheet exports pad trailing empty columns, which pandas auto-names "Unnamed: N".
    wide = wide.loc[:, [c for c in wide.columns if not str(c).startswith("Unnamed")]]
    # mag_id and sample_id are merge keys downstream -- a stray space would turn a match into a
    # silent miss, so strip both axes now.
    wide.index = wide.index.astype(str).str.strip()
    wide.columns = [str(c).strip() for c in wide.columns]

    _check(
        wide.notna().all().all(),
        f"{path} has missing values; Table S3 must be fully populated (absent genomes are 0, not NA).",
    )
    # No separate finiteness check: NaNs are excluded above, and an inf would blow the 100% column-sum
    # check below.
    col_sums = wide.sum()
    # 1e-6 on a 0-100 scale: loose enough for export rounding, far too tight for raw counts to pass.
    off = col_sums[(col_sums - 100.0).abs() > 1e-6]
    _check(
        off.empty,
        f"{path}: every sample column must sum to 100% (a composition). Offenders: "
        f"{off.head().to_dict()}",
    )
    logger.info(
        f"Table S3: {wide.shape[0]} genomes x {wide.shape[1]} samples, all columns sum to 100%."
    )
    # Long form (one row per genome x sample) so the design attaches with a single sample_id merge.
    return (
        wide.rename_axis("mag_id")
        .reset_index()
        .melt(id_vars="mag_id", var_name="sample_id", value_name="ra")
    )


def load_metadata(path: Path) -> pd.DataFrame:
    """Read the sample sheet, restrict to PRE/END, and assert the design is the balanced one.

    Everything downstream assumes a complete, balanced panel: every mouse present at both timepoints,
    every cage present in both arms, and ``N_REPLICATES`` cages. These are exactly the conditions
    under which the flat pool over (replicate, timepoint) equals ``(ra_pre + ra_end) / 2`` and under
    which the paired ``delta`` is defined for every mouse, so they are checked once, here, rather
    than being re-derived at each aggregation step.
    """
    # dtype=str keeps every ID literal: numeric-looking replicate/subject values would otherwise be
    # sniffed as int64 and silently stop matching their string twins in merges and group keys.
    md = pd.read_csv(path, sep="\t", dtype=str)
    need = ["sample_id", "subjectID", "replicate", "group", "time"]
    missing_cols = [c for c in need if c not in md.columns]
    _check(not missing_cols, f"{path} missing column(s): {missing_cols}")
    md = md[need].copy()
    for c in need:
        md[c] = md[
            c
        ].str.strip()  # every kept column is a key; stray whitespace would fork groups

    md = md[md["time"].isin([T1, TN])].copy()
    _check(not md.empty, f"{path} has no rows with time in {{{T1}, {TN}}}.")
    _check(
        md["sample_id"].is_unique,
        f"{path}: sample_id must be unique within the {T1}/{TN} subset.",
    )

    reps = sorted(md["replicate"].unique())
    _check(
        len(reps) == N_REPLICATES,
        f"expected {N_REPLICATES} replicates, found {len(reps)}: {reps}. The metadata-only shortcut "
        f"assumes the study's cage count; see DESIGN.md section 2.",
    )

    # Each mouse belongs to exactly one cage and one arm -- otherwise "the arm's mice in cage r" is
    # not well defined and step 2 would double-count.
    per_subj = md.groupby("subjectID")[["replicate", "group"]].nunique()
    bad = per_subj[(per_subj["replicate"] > 1) | (per_subj["group"] > 1)]
    _check(
        bad.empty, f"subjectID(s) spanning >1 replicate or group: {bad.index.tolist()}"
    )

    # Complete pairing: the paired delta needs both timepoints for every mouse.
    times = md.groupby("subjectID")["time"].agg(set)
    missing = times[times != {T1, TN}]
    _check(
        missing.empty,
        f"subjectID(s) missing a timepoint (need both {T1} and {TN}): {missing.to_dict()}",
    )

    # Both arms in every cage: required for the paired two-sample contrast to exist per replicate.
    arms = md.groupby("replicate")["group"].agg(set)
    lopsided = arms[arms.map(len) != 2]
    _check(
        lopsided.empty, f"replicate(s) not present in both arms: {lopsided.to_dict()}"
    )

    logger.info(
        f"Metadata: {md['subjectID'].nunique()} mice, {len(reps)} replicates, "
        f"{md['group'].nunique()} arms, {len(md)} {T1}/{TN} samples."
    )
    return md


def load_cells(cell_stats_path: Path) -> pd.DataFrame:
    """The tested cells: one row per (test_family, group, mag_id) we must produce abundance for.

    Cell stats is derived from ``p_value_summary``, which only ever contains MAGs that produced
    p-values -- so this *is* the tested set, and iterating over it is what keeps untested MAGs out of
    the output. Restricted here to ``COMPARISON`` and the two t-test families.
    """
    cs = pd.read_csv(cell_stats_path, sep="\t", dtype={"group_analyzed": str})
    # two_sample rows leave group_analyzed blank (read back as NaN); "" is the working sentinel for
    # "pooled across both arms" and doubles as a group/join key, so it must be a real string.
    cs["group_analyzed"] = cs["group_analyzed"].fillna("").astype(str).str.strip()
    cs = cs[
        (cs["comparison"] == COMPARISON) & (cs["test_type"].isin(TEST_FAMILY))
    ].copy()
    _check(
        not cs.empty,
        f"{cell_stats_path} has no rows for comparison={COMPARISON!r} and test types "
        f"{sorted(TEST_FAMILY)}.",
    )
    cs["test_family"] = cs["test_type"].map(TEST_FAMILY)
    # two_sample cells carry no arm; single_sample cells carry the arm they analyzed.
    cs["group"] = np.where(cs["test_family"] == "two_sample", "", cs["group_analyzed"])
    # One output row per unique (family, arm, MAG), however many source rows described the cell.
    cells = (
        cs[["test_family", "group", "mag_id"]].drop_duplicates().reset_index(drop=True)
    )
    logger.info(
        "Cells: "
        + ", ".join(
            f"{fam}/{grp or '(pooled)'}={n}"
            for (fam, grp), n in cells.groupby(["test_family", "group"]).size().items()
        )
    )
    return cells


def per_replicate_values(ra_long: pd.DataFrame, md: pd.DataFrame) -> pd.DataFrame:
    """Steps 1 and 2: mouse-level paired delta, then the unweighted collapse into (replicate, arm).

    Returns one row per ``(mag_id, replicate, group)`` with ``pre``, ``end``, ``delta`` and
    ``n_mice`` -- the frame every cell is then aggregated from.

    Step 2 is the unweighted mean that mirrors ``get_mean_change``: a cage with one mouse contributes
    exactly as much as a cage with two. That is the whole point, and it is why ``n_mice`` is carried
    (for reporting) but never used as a weight.
    """
    # Inner join drops Table S3 samples outside the PRE/END design (other timepoints/cohorts);
    # validate re-asserts, at the point of use, that no sample_id maps to two metadata rows.
    df = ra_long.merge(md, on="sample_id", how="inner", validate="many_to_one")
    _check(
        not df.empty,
        "no Table S3 sample matched the metadata after the PRE/END filter.",
    )

    # ---- step 1: per mouse, END - PRE ----
    # aggfunc="first" is a selector, not an aggregation: the design has one sample per
    # (mouse, timepoint), so each pivot cell sees at most one value. The check below catches
    # *absent* cells, not duplicates -- this line is where that design assumption is spent.
    wide = df.pivot_table(
        index=["mag_id", "subjectID", "replicate", "group"],
        columns="time",
        values="ra",
        aggfunc="first",
    ).reset_index()
    _check(
        {T1, TN}.issubset(wide.columns) and wide[[T1, TN]].notna().all().all(),
        f"every (mag, mouse) must have both a {T1} and an {TN} abundance after the pivot.",
    )
    wide["delta"] = wide[TN] - wide[T1]

    # ---- step 2: unweighted mean over the mice of that arm within each cage ----
    rep = wide.groupby(["mag_id", "replicate", "group"], as_index=False).agg(
        pre=(T1, "mean"),
        end=(TN, "mean"),
        delta=("delta", "mean"),
        n_mice=("subjectID", "nunique"),
    )
    # No delta-vs-(end - pre) check here: the pivot assertion above guarantees every mouse carries
    # both timepoints, which makes the two algebraically identical within each (mag, replicate,
    # group). The end-to-end version of this identity is still asserted once, in build_cells.
    logger.info(
        f"Replicate-level values: {len(rep):,} rows over {rep['mag_id'].nunique()} MAGs "
        f"({rep['n_mice'].min()}-{rep['n_mice'].max()} mice per cage-arm)."
    )
    return rep


def _aggregate(sub: pd.DataFrame) -> pd.Series:
    """Step 3 for one cell: collapse its replicate values into the reported columns.

    ``sub`` holds the replicate rows the cell's test used -- one arm's 8 for a single-sample cell,
    both arms' 16 for a two-sample cell.

    ``abs_mean_delta_ra`` and ``mean_abs_delta_ra`` differ whenever cages move in opposite
    directions: the first averages first and so lets them cancel, the second takes the absolute value
    first and so does not. They coincide only when every cage moves the same way.
    """
    ra_pre, ra_end = sub["pre"].mean(), sub["end"].mean()
    delta = sub["delta"].mean()
    return pd.Series(
        {
            "ra_pre": ra_pre,
            "ra_end": ra_end,
            # Mean of the two timepoint means. Equals a flat pool over the (replicate, timepoint)
            # values because both timepoints carry the same replicates -- verified in verify_flat_pool.
            "ra_mean": (ra_pre + ra_end) / 2.0,
            "delta_ra": delta,
            "abs_mean_delta_ra": abs(delta),
            "mean_abs_delta_ra": sub["delta"].abs().mean(),
            "n_replicates": sub["replicate"].nunique(),
            # Sum of per-(cage, arm) mouse counts -- reporting only, never a weight (see step 2).
            "n_mice": int(sub["n_mice"].sum()),
        }
    )


def build_cells(cells: pd.DataFrame, rep: pd.DataFrame) -> pd.DataFrame:
    """Step 3 across every cell, plus the two-sample within-replicate contrast columns.

    Single-sample cells aggregate their own arm's replicate rows. Two-sample cells aggregate both
    arms pooled, and additionally carry the difference-in-differences: within each replicate,
    ``c(r) = d_fat(r) - d_control(r)``, averaged over replicates. That contrast -- not the pooled
    delta -- is the abundance quantity shaped like what ``two_sample_paired`` tests, because that test
    compares the arms *inside* a replicate (``two_sample_paired.py:80-86``). A MAG that rose equally
    in both arms has a large pooled delta and a contrast near zero, and is exactly the MAG the
    divergence test should not flag.
    """
    arms = sorted(rep["group"].unique())
    _check(arms == ["control", "fat"], f"expected arms ['control','fat'], found {arms}")

    out = []
    for (fam, grp), grp_cells in cells.groupby(["test_family", "group"]):
        wanted = set(grp_cells["mag_id"])
        # A single-sample cell sees only its own arm; a two-sample cell sees both.
        sub = rep[rep["mag_id"].isin(wanted)]
        if fam == "single_sample":
            sub = sub[sub["group"] == grp]
        # include_groups=False keeps mag_id out of the frame handed to _aggregate (and opts out of
        # the pandas-2.2-deprecated behavior of passing group keys along).
        agg = (
            sub.groupby("mag_id").apply(_aggregate, include_groups=False).reset_index()
        )

        missing = wanted - set(agg["mag_id"])
        _check(
            not missing,
            f"{fam}/{grp or '(pooled)'}: no abundance for MAG(s) {sorted(missing)[:5]}",
        )

        expect_reps = N_REPLICATES
        _check(
            (agg["n_replicates"] == expect_reps).all(),
            f"{fam}/{grp or '(pooled)'}: every cell must aggregate {expect_reps} replicates; found "
            f"{sorted(agg['n_replicates'].unique())}.",
        )
        agg["test_family"], agg["group"] = fam, grp
        out.append(agg)

    cell_df = pd.concat(out, ignore_index=True)

    # ---- two-sample contrast: form the difference INSIDE each replicate, then average ----
    # One row per (mag, replicate) with a `fat` and a `control` column. aggfunc="first" is a pure
    # selector here: rep is unique on (mag_id, replicate, group) by construction of its groupby.
    d = rep.pivot_table(
        index=["mag_id", "replicate"], columns="group", values="delta", aggfunc="first"
    )
    _check(
        d.notna().all().all(),
        "every (mag, replicate) needs a delta in BOTH arms to form the paired contrast; "
        "some are missing.",
    )
    d["c"] = d["fat"] - d["control"]
    contrast = d.groupby("mag_id").agg(
        delta_ra_fat=("fat", "mean"),
        delta_ra_control=("control", "mean"),
        delta_ra_contrast=("c", "mean"),
        # mean of |within-cage contrast|: how far apart the arms moved in a typical cage,
        # regardless of which direction. Large here with a small delta_ra_contrast is the
        # cages-disagree case -- high variance, which is what weakens a paired t-test.
        mean_abs_delta_ra_contrast=("c", lambda v: v.abs().mean()),
    )
    contrast["abs_mean_delta_ra_contrast"] = contrast["delta_ra_contrast"].abs()
    # No contrast-vs-(fat - control) check: the d.notna assertion above pins both arms to the same
    # replicate set, under which averaging commutes with the subtraction exactly.

    is_two = cell_df["test_family"] == "two_sample"
    cell_df = cell_df.merge(contrast.reset_index(), on="mag_id", how="left")
    # The contrast is a property of the two-sample test only; blank it elsewhere.
    for c in [
        "delta_ra_fat",
        "delta_ra_control",
        "delta_ra_contrast",
        "abs_mean_delta_ra_contrast",
        "mean_abs_delta_ra_contrast",
    ]:
        cell_df.loc[~is_two, c] = np.nan

    # ---- the identity that catches most aggregation mistakes ----
    drift = (cell_df["delta_ra"] - (cell_df["ra_end"] - cell_df["ra_pre"])).abs().max()
    _check(drift < FLOAT_TOL, f"delta_ra != ra_end - ra_pre (max |diff| = {drift:.3e})")

    cell_df["comparison"] = COMPARISON
    # "pre_end-fat_control" -> period="pre_end", group_pair="fat_control" (split at the first "-").
    cell_df["period"], _, cell_df["group_pair"] = COMPARISON.partition("-")
    # reindex pins the published column order and drops any helper columns not in OUT_COLS.
    return (
        cell_df.reindex(columns=OUT_COLS)
        .sort_values(["test_family", "group", "mag_id"])
        .reset_index(drop=True)
    )


def verify_flat_pool(
    cell_df: pd.DataFrame, cells: pd.DataFrame, rep: pd.DataFrame
) -> None:
    """Confirm ``ra_mean`` equals a flat pool over the (replicate, timepoint) values.

    DRIDO's ``both`` variation flat-pools every observation rather than nesting means, precisely
    because its panel is ragged. Here the two definitions coincide -- but only while every replicate
    appears at both timepoints. Checking it makes the equivalence a verified property rather than an
    assumption, and it fails under exactly the condition that also invalidates the paired delta.
    """
    worst = 0.0
    for (fam, grp), grp_cells in cells.groupby(["test_family", "group"]):
        sub = rep[rep["mag_id"].isin(set(grp_cells["mag_id"]))]
        if fam == "single_sample":
            sub = sub[sub["group"] == grp]
        # Flat pool: stack the pre and end columns into one vector of observations, then one mean.
        flat = (
            sub.melt(id_vars="mag_id", value_vars=["pre", "end"], value_name="v")
            .groupby("mag_id")["v"]
            .mean()
        )
        got = cell_df[
            (cell_df["test_family"] == fam) & (cell_df["group"] == grp)
        ].set_index("mag_id")["ra_mean"]
        # reindex aligns the recomputed pool to the cell rows so the subtraction is by mag_id,
        # not by row position.
        worst = max(worst, float((got - flat.reindex(got.index)).abs().max()))
    _check(
        worst < FLOAT_TOL,
        f"ra_mean != flat pool over (replicate, timepoint) values (max |diff| = {worst:.3e}); "
        f"some replicate is missing a timepoint.",
    )
    logger.info(
        f"ra_mean == flat pool over (replicate, timepoint): max |diff| = {worst:.3e}."
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Per-cell relative abundance and PRE->END change for the diet AlleleFlux run.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    base = Path("/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026")
    # ArgumentDefaultsHelpFormatter only appends "(default: ...)" to arguments that HAVE a help
    # string -- without one, -h shows nothing, so every argument here must carry help=.
    ap.add_argument(
        "--table-s3",
        type=Path,
        default=base / "relative_abundance/TableS3_goldman.tsv",
        help="published per-sample relative-abundance matrix (Table S3)",
    )
    ap.add_argument(
        "--metadata",
        type=Path,
        default=base / "metadata_md_bam_fixed.tsv",
        help="sample sheet with sample_id, subjectID, replicate, group, time",
    )
    ap.add_argument(
        "--cell-stats",
        type=Path,
        default=base
        / "relative_abundance/significant_sites_summary/significant_sites_mag_cell_stats_long.tsv",
        help="significant-sites cell stats defining the tested (test, group, MAG) cells",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).parent / "tables",
        help="directory for rel_abundance_by_cell.tsv",
    )
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    # A missing file is a usage error, not a data-integrity violation -- plain exit, no _check.
    for p in (args.table_s3, args.metadata, args.cell_stats):
        if not p.exists():
            logger.error(f"Input does not exist: {p}")
            sys.exit(1)

    ra_long = load_table_s3(args.table_s3)
    md = load_metadata(args.metadata)

    # Table S3 and the metadata must describe the same samples: an unmatched sample on either side
    # means the composition and the design are from different runs.
    s3_samples, md_samples = set(ra_long["sample_id"]), set(md["sample_id"])
    _check(
        md_samples <= s3_samples,
        f"metadata sample(s) absent from Table S3: {sorted(md_samples - s3_samples)[:5]}",
    )

    cells = load_cells(args.cell_stats)

    # Every tested MAG needs a composition row, or its abundance would silently be absent.
    orphans = set(cells["mag_id"]) - set(ra_long["mag_id"])
    _check(not orphans, f"tested MAG(s) with no Table S3 row: {sorted(orphans)[:5]}")

    # Restrict Table S3 to the tested MAGs up front: untested genomes never reach the output, and
    # the per-mouse pivot scales with rows.
    rep = per_replicate_values(
        ra_long[ra_long["mag_id"].isin(set(cells["mag_id"]))], md
    )
    cell_df = build_cells(cells, rep)
    verify_flat_pool(cell_df, cells, rep)

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / "rel_abundance_by_cell.tsv"
    # %.10g keeps ~10 significant digits -- diff-stable output without float representation noise.
    cell_df.to_csv(out, sep="\t", index=False, float_format="%.10g")
    logger.info(
        f"Wrote {out} ({len(cell_df)} rows, {cell_df['mag_id'].nunique()} MAGs)."
    )
    logger.info(
        "\n"
        + cell_df.groupby(["test_family", "group"])
        .agg(
            n_MAGs=("mag_id", "nunique"),
            ra_mean_median=("ra_mean", "median"),
            delta_ra_median=("delta_ra", "median"),
            n_mice=("n_mice", "first"),
        )
        .to_string()
    )


if __name__ == "__main__":
    main()
