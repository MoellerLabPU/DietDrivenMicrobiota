#!/usr/bin/env python3
"""Figures for the variable-site tables.

Every figure is a function taking DataFrames and returning a matplotlib Figure, so
they can be called individually from a REPL as well as written to disk by the CLI.
Nothing is computed here that the tables do not already contain -- these are views
of the analysis output, not a second analysis.

    python plots.py --indir <stats dir> --outdir <figs dir>
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, LogNorm


def setup_logging(level=logging.INFO) -> None:
    """Configure logging once, in main(). Same format as the analysis script."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


logger = logging.getLogger(__name__)

# Okabe-Ito colourblind-safe triple for the three real tests; validated to pass the
# lightness band, chroma floor, all-pairs CVD separation and contrast checks against
# both a light and a dark surface. "union" is deliberately NEUTRAL rather than a
# fourth hue: it is an aggregate of the other three, not a peer category, and adding
# a fourth hue pushed the worst CVD pair to deltaE 7.6 (below the 8 floor).
COLORS = {"div": "#D55E00", "hf": "#009E73", "lf": "#0072B2", "union": "#666666"}
# Display names. The TABLE columns keep the short keys (div_/hf_/lf_/union_) --
# these are for axes, titles and legends only.
LABELS = {
    "div": "divergence (high fat vs low fat)",
    "hf": "high-fat parallelism",
    "lf": "low-fat parallelism",
    "union": "union of tested sets",
}
# Compact form for panel titles and crowded legends.
SHORT = {"div": "divergence", "hf": "high fat", "lf": "low fat", "union": "union"}
DEFINITIONS = ("div", "hf", "lf", "union")
SIG_SETS = ("div", "hf", "lf")

# Where variable_site_distribution.py writes its tables, and where the figures go
# beside them. Defaults rather than required flags because there is one canonical
# location; pass --indir/--outdir to point at a different run.
INDIR = Path("/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/variable_site_stats")
OUTDIR = INDIR / "figs"

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "legend.frameon": False,
})


def load(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The three tables. contigs/<MAG>.tsv is not read: it is a split of contig_level_all."""
    sgb = pd.read_csv(indir / "sgb_level.tsv", sep="\t")
    contig = pd.read_csv(indir / "contig_level_all.tsv", sep="\t")
    summary = pd.read_csv(indir / "summary.tsv", sep="\t")
    logger.info(
        f"loaded {len(sgb)} SGBs, {len(contig):,} contigs, "
        f"{len(summary)} summary rows from {indir}"
    )
    return sgb, contig, summary


def _identity_line(ax, lo: float, hi: float) -> None:
    """y = x. Everything below it is 'observed smaller than expected'."""
    ax.plot([lo, hi], [lo, hi], ls="--", lw=1, color="0.45", zorder=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)


# --------------------------------------------------------------------- SGB level

def fig_significance_ranked(sgb: pd.DataFrame):
    """Per-SGB share of tested contigs carrying a significant site (Q5/Q7/Q9).

    Each test is a fraction of ITS OWN tested contigs, which is why the three series
    are comparable at all despite testing different numbers of contigs.
    """
    d = sgb.sort_values("div_pct_contigs_sig", ascending=False).reset_index(drop=True)
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(11, 4.2))

    # Vertical connectors first so the markers sit on top of them.
    ax.vlines(x, d.div_pct_contigs_sig, d.hf_pct_contigs_sig,
              color="0.8", lw=0.8, zorder=1)
    for k in ("div", "hf"):
        ax.scatter(x, d[f"{k}_pct_contigs_sig"], s=26, color=COLORS[k],
                   label=LABELS[k], zorder=3, edgecolor="white", linewidth=0.5)
    # Low fat is drawn as a flat line at zero rather than as a third scatter series,
    # and the annotation states that as fact -- so assert it instead of trusting it.
    # A different comparison, a looser threshold or per-SGB FDR could all make low fat
    # non-zero, and the figure would go on drawing a line at 0 and captioning it.
    if not (d.lf_pct_contigs_sig == 0).all():
        raise ValueError(
            "low-fat parallelism is no longer zero everywhere "
            f"(max {d.lf_pct_contigs_sig.max():.2f}%) -- this figure hardcodes it at 0. "
            "Plot it as a third series before using this again."
        )
    ax.axhline(0, color=COLORS["lf"], lw=2, zorder=2)
    # Bottom-left is the only reliably empty corner: the series descends rightward.
    ax.annotate(
        f"low-fat parallelism: 0% for all {len(d)} SGBs\n"
        "(no low-fat site reaches q < 0.05)",
        xy=(0.4, 3.5), fontsize=8, color=COLORS["lf"], va="bottom", ha="left",
    )

    for k in ("div", "hf"):
        ax.axhline(d[f"{k}_pct_contigs_sig"].median(), color=COLORS[k],
                   ls=":", lw=1, alpha=0.7)

    ax.set_xlabel(f"{len(d)} SGBs, ranked by divergence")
    ax.set_ylabel("% of that test's tested contigs\nwith a significant site")
    ax.set_title("Significance reaches a minority of contigs, and never in the low-fat group",
                 loc="left", fontsize=10)
    ax.set_xticks([])
    ax.set_ylim(-3, 105)
    ax.legend(loc="upper right")
    fig.text(0.005, -0.02,
             f"dotted lines = medians ({d.div_pct_contigs_sig.median():.0f}% divergence, "
             f"{d.hf_pct_contigs_sig.median():.0f}% high fat)", fontsize=7, color="0.4")
    return fig


def fig_gap_obs_vs_exp(sgb: pd.DataFrame):
    """Observed spacing against the uniform null, per SGB, per definition (Q3 vs Q4)."""
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), sharex=True, sharey=True)
    vals = pd.concat([sgb[[f"{k}_mean_gap", f"{k}_expected_gap"]].stack()
                      for k in DEFINITIONS]).replace(0, np.nan).dropna()
    lo, hi = vals.min() * 0.6, vals.max() * 1.6

    for ax, k in zip(axes, DEFINITIONS):
        o, e = sgb[f"{k}_mean_gap"], sgb[f"{k}_expected_gap"]
        ok = o.notna() & e.notna()
        ax.scatter(e[ok], o[ok], s=20, color=COLORS[k], alpha=0.8,
                   edgecolor="white", linewidth=0.4)
        _identity_line(ax, lo, hi)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ratio = (o[ok] / e[ok]).median()
        ax.set_title(f"{SHORT[k]}  (n={int(ok.sum())}, median ratio {ratio:.2f})",
                     fontsize=9)
        ax.set_xlabel("expected gap, uniform null (nt, log)")
    axes[0].set_ylabel("observed mean gap (nt, log)")
    fig.suptitle("Tested sites sit closer together than random placement predicts", fontsize=10, x=0.007, ha="left")
    return fig


def fig_div_vs_hf(sgb: pd.DataFrame):
    """Do divergence and high-fat parallelism single out the same SGBs?"""
    fig, ax = plt.subplots(figsize=(5.2, 5))
    x, y = sgb.hf_pct_contigs_sig, sgb.div_pct_contigs_sig
    ok = x.notna() & y.notna()
    ax.scatter(x[ok], y[ok], s=34, color=COLORS["div"], alpha=0.75,
               edgecolor="white", linewidth=0.5)
    _identity_line(ax, -3, 105)
    ax.set_xlabel("% high-fat-tested contigs with a high-fat-significant site")
    ax.set_ylabel("% divergence-tested contigs\nwith a divergence-significant site")
    ax.set_title(f"r = {x[ok].corr(y[ok]):+.2f}  (n={int(ok.sum())} SGBs)",
                 loc="left", fontsize=10)
    ax.annotate("above the line:\nmore divergence than high fat", (4, 92),
                fontsize=7, color="0.45")
    return fig


def fig_contig_nesting(sgb: pd.DataFrame):
    """reference -> tested -> significant, one line per SGB (Q2, Q1, Q5).

    Divergence only: three definitions overlaid would be unreadable at 62 lines.
    """
    fig, ax = plt.subplots(figsize=(5.6, 5))
    cols = ["n_ref_contigs", "div_n_contigs_with_sites", "div_n_contigs_sig"]
    for _, r in sgb.iterrows():
        # 0 cannot be drawn on a log axis; nudge to 0.5 and mark the axis accordingly.
        ax.plot([0, 1, 2], [max(r[c], 0.5) for c in cols],
                color=COLORS["div"], alpha=0.28, lw=0.9, marker="o", ms=2.6)
    med = [sgb[c].median() for c in cols]
    ax.plot([0, 1, 2], med, color="black", lw=2.4, marker="o", ms=6, zorder=5,
            label="median SGB")
    for i, v in enumerate(med):
        ax.annotate(f"{v:,.0f}", (i, v), xytext=(6, 6), textcoords="offset points",
                    fontsize=8, fontweight="bold")
    ax.set_yscale("log")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["in the\nreference", "with a\ntested site",
                        "with a significant\nsite"])
    ax.set_ylabel("contigs per SGB (log)")
    ax.set_title("Most contigs are never tested, and most tested contigs are not significant",
                 loc="left", fontsize=9.5)
    # The reference stage is definition-independent; the other two are divergence.
    # Upper right: the lines descend left-to-right, so that corner is always free.
    # Lower centre collides with the legend.
    ax.annotate("tested / significant = divergence\n(two_sample_paired_tTest, q < 0.05)",
                xy=(0.98, 0.97), xycoords="axes fraction", ha="right", va="top",
                fontsize=7.5, color=COLORS["div"])
    ax.legend(loc="lower left")
    fig.text(0.005, -0.03, "zeros drawn at 0.5 so they are visible on a log axis",
             fontsize=7, color="0.4")
    return fig


# ------------------------------------------------------------------ contig level

def fig_len_vs_sites(contig: pd.DataFrame):
    """Does the amount of testing track contig length?"""
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharex=True, sharey=True)
    for ax, k in zip(axes, SIG_SETS):
        sub = contig[contig[f"{k}_n_sites"] > 0]
        ax.scatter(sub.contig_len, sub[f"{k}_n_sites"], s=5, alpha=0.25,
                   color=COLORS[k], edgecolor="none")
        r = np.corrcoef(np.log10(sub.contig_len), np.log10(sub[f"{k}_n_sites"]))[0, 1]
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{SHORT[k]}  (n={len(sub):,} contigs, log-log r={r:+.2f})",
                     fontsize=9)
        ax.set_xlabel("contig length (bp, log)")
    axes[0].set_ylabel("tested sites on the contig (log)")
    fig.suptitle("Tested sites per contig against contig length",
                 fontsize=10, x=0.007, ha="left")
    return fig


def fig_gap_ecdf(contig: pd.DataFrame):
    """Distribution of per-contig mean gaps, observed against the uniform null."""
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, k in zip(axes, SIG_SETS):
        for col, ls, lab in ((f"{k}_mean_gap", "-", "observed"),
                             (f"{k}_expected_gap", "--", "uniform null")):
            # ax.ecdf sorts internally and draws a proper step, collapsing tied
            # values into one step rather than a vertical stack of points.
            v = contig.loc[contig[f"{k}_n_gaps"] > 0, col].dropna()
            v = v[v > 0]
            if not len(v):
                continue
            ax.ecdf(v, ls=ls, lw=1.6, color=COLORS[k], label=lab,
                    alpha=1 if ls == "-" else 0.55)
        ax.set_xscale("log")
        ax.set_title(f"{SHORT[k]}  "
                     f"(n={int((contig[f'{k}_n_gaps'] > 0).sum()):,} contigs)",
                     fontsize=9)
        ax.set_xlabel("mean gap on a contig (nt, log)")
        ax.legend(loc="lower right", fontsize=8)
    axes[0].set_ylabel("cumulative fraction of contigs")
    fig.suptitle("Observed spacing is shifted left of the null at every quantile",
                 fontsize=10, x=0.007, ha="left")
    return fig


# --------------------------------------------------------------------- per-MAG
# One page per SGB. The two spacing views from the SGB-level figures, but with
# CONTIGS as the unit and one panel per definition, beside that SGB's numbers.


def _mag_gap_scatter(ax, cr: pd.DataFrame, k: str, lo: float, hi: float) -> None:
    """Observed vs uniform-null mean gap, one point per contig (cf. fig_gap_obs_vs_exp)."""
    o, e = cr[f"{k}_mean_gap"], cr[f"{k}_expected_gap"]
    ok = o.notna() & e.notna() & (o > 0) & (e > 0)
    if not ok.any():
        ax.text(0.5, 0.5, "no contig has\n>=2 tested sites", ha="center", va="center",
                fontsize=7.5, color="0.5", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(SHORT[k], fontsize=9, color=COLORS[k])
        return
    ax.scatter(e[ok], o[ok], s=16, color=COLORS[k], alpha=0.75,
               edgecolor="white", linewidth=0.3)
    ax.plot([lo, hi], [lo, hi], ls="--", lw=1, color="0.45", zorder=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ratio = (o[ok] / e[ok]).median()
    ax.set_title(f"{SHORT[k]}  (n={int(ok.sum())}, median {ratio:.2f}x)",
                 fontsize=8.5, color=COLORS[k])
    ax.set_xlabel("expected gap (nt, log)", fontsize=8)


def _mag_gap_ecdf(ax, cr: pd.DataFrame, k: str) -> None:
    """Cumulative distribution of per-contig mean gaps, observed vs null (cf. fig_gap_ecdf)."""
    has = cr[f"{k}_n_gaps"] > 0
    if not has.any():
        ax.text(0.5, 0.5, "no contig has\n>=2 tested sites", ha="center", va="center",
                fontsize=7.5, color="0.5", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    for col, ls, lab in ((f"{k}_mean_gap", "-", "observed"),
                         (f"{k}_expected_gap", "--", "uniform null")):
        v = cr.loc[has, col].dropna()
        v = v[v > 0]
        if not len(v):
            continue
        # ax.ecdf draws a step, so with only a few contigs the staircase stays
        # visible -- honest about how little the curve rests on.
        ax.ecdf(v, lw=1.6, ls=ls, color=COLORS[k],
                alpha=1 if ls == "-" else 0.55, label=lab)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("mean gap on a contig (nt, log)", fontsize=8)
    ax.set_title(f"{SHORT[k]}  (n={int(has.sum())} contigs)", fontsize=8.5,
                 color=COLORS[k])
    ax.legend(fontsize=6.5, loc="lower right")


def mag_card(mag: str, cr: pd.DataFrame, row: pd.Series):
    """One page per SGB: spacing vs null and the gap ECDF, split by definition."""
    fig = plt.figure(figsize=(16.5, 7.4))
    gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 1, 0.9],
                          hspace=0.45, wspace=0.32)

    # A shared range across all four scatter panels, so they can be read against
    # each other rather than each silently rescaling to its own data.
    vals = pd.concat([cr[[f"{k}_mean_gap", f"{k}_expected_gap"]].stack()
                      for k in DEFINITIONS]).dropna()
    vals = vals[vals > 0]
    lo, hi = (vals.min() * 0.5, vals.max() * 2) if len(vals) else (1, 10)

    for col, k in enumerate(DEFINITIONS):
        ax = fig.add_subplot(gs[0, col])
        _mag_gap_scatter(ax, cr, k, lo, hi)
        if col == 0:
            ax.set_ylabel("observed mean gap (nt, log)", fontsize=8)
        ax = fig.add_subplot(gs[1, col])
        _mag_gap_ecdf(ax, cr, k)
        if col == 0:
            ax.set_ylabel("cumulative fraction\nof contigs", fontsize=8)

    fig.text(0.005, 0.955, "spacing vs uniform null, one point per contig — "
             "below y = x means clustered", fontsize=8.5, color="0.35")
    fig.text(0.005, 0.475, "distribution of per-contig mean gaps, observed vs null",
             fontsize=8.5, color="0.35")

    tx = fig.add_subplot(gs[:, 4])
    tx.axis("off")
    lines = [f"contigs in this table:   {len(cr):,}",
             f"reference contigs:       {int(row.n_ref_contigs):,}",
             f"reference length:        {int(row.ref_genome_len):,} bp", ""]
    for k in DEFINITIONS:
        lines.append(f"[{SHORT[k]}]")
        lines.append(f"  contigs tested:        {int(row[f'{k}_n_contigs_with_sites']):,}")
        lines.append(f"  tested sites:          {int(row[f'{k}_n_sites']):,}")
        mg, eg = row[f"{k}_mean_gap"], row[f"{k}_expected_gap"]
        lines.append(f"  mean gap:              {mg:,.0f} nt" if pd.notna(mg)
                     else "  mean gap:              n/a")
        lines.append(f"  expected gap:          {eg:,.0f} nt" if pd.notna(eg)
                     else "  expected gap:          n/a")
        if pd.notna(mg) and pd.notna(eg) and eg:
            lines.append(f"  observed / expected:   {mg / eg:.2f}")
        if k in SIG_SETS:
            lines.append(f"  contigs significant:   {int(row[f'{k}_n_contigs_sig']):,} "
                         f"({row[f'{k}_pct_contigs_sig']:.1f}%)")
            lines.append(f"  significant sites:     {int(row[f'{k}_n_sig_sites']):,}")
        lines.append("")
    tx.text(0, 1, "\n".join(lines), va="top", ha="left", family="monospace",
            fontsize=7.2, transform=tx.transAxes)

    fig.suptitle(mag, x=0.005, ha="left", fontsize=11, y=0.995)
    return fig


def write_cards(sgb: pd.DataFrame, contig: pd.DataFrame, dest: Path) -> int:
    by_mag = dict(tuple(contig.groupby("MAG_ID", observed=True)))
    n = 0
    with PdfPages(dest) as pdf:
        for _, row in sgb.sort_values("div_n_sites", ascending=False).iterrows():
            cr = by_mag.get(row.MAG_ID)
            if cr is None or cr.empty:
                logger.warning(f"{row.MAG_ID}: no contig rows — page skipped")
                continue
            fig = mag_card(row.MAG_ID, cr, row)
            pdf.savefig(fig)
            plt.close(fig)
            n += 1
    return n


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Figures for the variable-site tables. Reads only the analysis "
                    "outputs; computes no new statistics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--indir", type=Path, default=INDIR,
                        help="Directory holding sgb_level.tsv, contig_level_all.tsv "
                             "and summary.tsv.")
    parser.add_argument("--outdir", type=Path, default=OUTDIR,
                        help="Directory to write the PNGs and the per-MAG PDF into.")
    args = parser.parse_args()

    # Non-interactive backend for file output, set HERE and not at import time: the
    # notebook imports this module, and forcing Agg on import would silently disable
    # its inline rendering. Switching after the pyplot import is supported.
    matplotlib.use("Agg")

    setup_logging()
    args.outdir.mkdir(parents=True, exist_ok=True)
    sgb, contig, _ = load(args.indir)

    figures = {
        "sgb_significance_ranked": lambda: fig_significance_ranked(sgb),
        "sgb_gap_obs_vs_exp": lambda: fig_gap_obs_vs_exp(sgb),
        "sgb_div_vs_hf": lambda: fig_div_vs_hf(sgb),
        "sgb_contig_nesting": lambda: fig_contig_nesting(sgb),
        "contig_len_vs_sites": lambda: fig_len_vs_sites(contig),
        "contig_gap_ecdf": lambda: fig_gap_ecdf(contig),
    }
    for name, build in figures.items():
        fig = build()
        dest = args.outdir / f"{name}.png"
        fig.savefig(dest)
        plt.close(fig)
        logger.info(f"wrote {dest.name}")

    dest = args.outdir / "per_mag_cards.pdf"
    logger.info(f"wrote {dest.name} ({write_cards(sgb, contig, dest)} pages)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
