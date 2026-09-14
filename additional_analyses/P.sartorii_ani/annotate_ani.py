#!/usr/bin/env python3
"""Attach mouse metadata to skani's output and classify every genome pair.

skani hands back rows of (genome_A, genome_B, ANI) with filenames like
SLG888_A4_10473203_HT53MAFX5.fna.gz. That is unanswerable as it stands -- you
cannot see which mouse, which timepoint, or which diet. This turns it into a
table you can reason about.

The question: did the same P. sartorii strain persist in a mouse from
Pre-Treatment to End, in both Fat and Control?

The trick is that a high within-mouse ANI proves nothing on its own. If two
UNRELATED isolates also score 99.9 on these assemblies, then 99.9 is simply
what this species scores against itself here and the persistence claim is
empty. So every pair is labelled, and the between-mouse pairs are kept as the
null to compare against.

Usage: 04_annotate.py <metadata.tsv> <skani_pairs.tsv> <pairs_annotated.tsv>
"""

import logging
import os
import sys

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(message)s")

meta_in, pairs_in, out_path = sys.argv[1:4]

TIMEPOINTS = ["Pre-Treatment", "End"]

# ---------------------------------------------------------------------------
# 1. Metadata, indexed by the filename skani will report back to us.
# ---------------------------------------------------------------------------
meta = pd.read_csv(meta_in, sep="\t")
meta = meta[meta["Host.Collection.Timepoint"].isin(TIMEPOINTS)]
meta = meta.set_index("Assembly.File")

# ---------------------------------------------------------------------------
# 2. skani output. It echoes back the full paths it was given, so strip them
#    to bare filenames so they match the metadata index.
# ---------------------------------------------------------------------------
# float_precision="round_trip": pandas' DEFAULT parser is off by one ULP on
# long float strings -- it reads 99.55754999999999 as 99.55755, a different
# double -- which would silently alter ANI values on the way through.
df = pd.read_csv(pairs_in, sep="\t", float_precision="round_trip")
df["a"] = df["Ref_file"].map(os.path.basename)
df["b"] = df["Query_file"].map(os.path.basename)

# Short names. The two percentile columns only exist if skani ran with --ci,
# so they are picked up conditionally rather than assumed.
rename = {"ANI": "ani",
          "Align_fraction_ref": "af_a",
          "Align_fraction_query": "af_b",
          "ANI_5_percentile": "ani_lo",
          "ANI_95_percentile": "ani_hi"}
df = df[["a", "b"] + [c for c in rename if c in df.columns]].rename(columns=rename)

# ---------------------------------------------------------------------------
# 3. Attach mouse, diet and timepoint to BOTH sides of each pair. A pair only
#    means something once you know both ends of it.
# ---------------------------------------------------------------------------
fields = {"Isolate.ID": "iso", "Host.ID": "host",
          "Host.Diet": "diet", "Host.Collection.Timepoint": "tp"}
for side in ("a", "b"):
    for col, short in fields.items():
        df[f"{short}_{side}"] = df[side].map(meta[col])

# A genome outside the Pre/End set maps to NaN; drop those pairs rather than
# letting them fall through into a comparison class.
df = df.dropna(subset=["host_a", "host_b"])

# ---------------------------------------------------------------------------
# 4. Classify. This is the analytical step -- everything else is plumbing.
# ---------------------------------------------------------------------------
def classify(r):
    """Which of the four kinds of comparison is this pair?"""
    # Does the pair span the two timepoints, or sit inside one?
    crosses = {r["tp_a"], r["tp_b"]} == {"Pre-Treatment", "End"}

    if r["host_a"] == r["host_b"]:
        # Same mouse. Across time = persistence, the question being asked.
        # Same timepoint = diversity co-existing in one gut at one moment,
        # which sets the floor for how similar two "different" isolates are.
        return "within_host_across_time" if crosses else "within_host_same_time"

    # Different mice: the null. Same diet vs different diet is split out in
    # case diet itself structures which strains a mouse carries.
    if r["diet_a"] == r["diet_b"]:
        return "between_host_same_diet"
    return "between_host_diff_diet"


df["comparison"] = df.apply(classify, axis=1)

# Diet of the pair, for the Fat vs Control contrast. "mixed" for pairs that
# span diets, which by definition are between-mouse.
df["diet"] = df.apply(
    lambda r: r["diet_a"] if r["diet_a"] == r["diet_b"] else "mixed", axis=1)

# A pair is only as trustworthy as its worse-covered genome, so keep the lower
# of the two aligned fractions. This is the column to filter on: a 99.9% ANI
# computed over 20% of a genome is not the same evidence as one over 95%.
df["af_min"] = df[["af_a", "af_b"]].min(axis=1)

df = df.sort_values("ani", ascending=False)
df.to_csv(out_path, sep="\t", index=False)
logging.info(f"{len(df)} pairs -> {out_path}\n")

# ---------------------------------------------------------------------------
# 5. Summarise.
# ---------------------------------------------------------------------------
logging.info("=== ANI by comparison class ===")
logging.info(df.groupby("comparison")
        .agg(n=("ani", "size"), ani_min=("ani", "min"),
             ani_med=("ani", "median"), ani_max=("ani", "max"),
             af_min=("af_min", "min")).to_string())

w = df[df["comparison"] == "within_host_across_time"]
logging.info("\n=== within-host Pre vs End, per mouse (the key contrast) ===")
if w.empty:
    logging.info("none -- no mouse has isolates at both timepoints")
else:
    logging.info(w.groupby(["diet", "host_a"])
           .agg(n_pairs=("ani", "size"), ani_min=("ani", "min"),
                ani_med=("ani", "median"), af_min=("af_min", "min"))
           .to_string())

    # The comparison that actually carries the argument.
    bet = df[df["comparison"].str.startswith("between_host")]
    logging.info("\n=== separation from the null ===")
    logging.info(f"within-host  Pre x End : median {w['ani'].median()}  "
          f"min {w['ani'].min()}  (n={len(w)})")
    logging.info(f"between-host           : median {bet['ani'].median()}  "
          f"max {bet['ani'].max()}  (n={len(bet)})")
    logging.info("\nIf the within-host minimum sits above the between-host maximum")
    logging.info("the two distributions do not overlap, which is a much stronger")
    logging.info("statement than clearing a fixed 99.9% threshold -- 99.9 is at the")
    logging.info("resolution limit of skani, so a threshold alone cannot carry it.")
    logging.info("Pairs are pseudoreplicates within a few mice, so any significance")
    logging.info("test must permute at the MOUSE level, not the pair level.")
