#!/usr/bin/env python3
"""Sites significant for BOTH divergence and high-fat parallelism, per species.

For a reviewer response. Counts, for each species in SPECIES below:
  1. sites with q < 0.05 in the divergence test AND q < 0.05 in HF parallelism
  2. contigs carrying at least one such site

A site qualifies only if the SAME (contig, position) clears q < 0.05 in both tests,
which is an inner join on position -- not something the per-contig counts in
sgb_level.tsv can give, since those store totals and discard the pairing. A contig
with 65 divergence-significant and 84 HF-significant sites could have anywhere from
0 to 65 in common.

Writes one TSV and logs the per-species averages.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd


def setup_logging(level=logging.INFO) -> None:
    """Configure logging once, in main(). Same format as the analysis script."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


logger = logging.getLogger(__name__)

# Same pinned run as variable_site_distribution.py. Several AlleleFlux trees hold
# identically named files with different contents, so this is not a flag.
RUN = Path("/scratch/gpfs/AMOELLER/diet_manip/AlleleFlux_mapq20/longitudinal")
COMPARISON = "pre_end-fat_control"
OUTDIR = Path("/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026")

# Carries both MAG_ID and the GTDB species labels, so it is the species -> MAG map.
FIG1_SCORES = (
    RUN
    / "scores/processed/combined/MAG"
    / "scores_two_sample_paired-pre_end-fat_control-MAGs.tsv"
)

DIVERGENCE_TEST = "two_sample_paired_tTest"  # between-group: HF vs LF
PARALLELISM_TEST = "single_sample_tTest"  # within-group
HF_GROUP = "fat"
LF_GROUP = "control"
# BH-corrected. NOTE the correction is pooled genome-wide across all MAGs, not per
# species, so "q < 0.05 for this species" means significant in a global ranking.
Q_THRESHOLD = 0.05

CHUNK = 2_000_000

# The species to report on, supplied by the reviewer request: those that did NOT
# show strain replacement in >50% of HF replicates. The leading number is their
# numbering, carried through so the output can be pasted back alongside theirs.
SPECIES = [
    ("110_1", "Odoribacter sp910578105"),
    ("111_1", "Odoribacter sp910589025"),
    ("113_1", "Odoribacter sp910578985"),
    ("114_1", "Pullibacteroides sp910586395"),
    ("115_1", "Pullibacteroides sp910578275"),
    ("119_1", "UBA7173 sp001689485"),
    ("127_1", "UBA3263 sp001689615"),
    ("128_1", "Phocaeicola sartorii"),
    ("129_1", "Bacteroides sp910586915"),
    ("12_1", "Ligilactobacillus murinus"),
    ("131_1", "Bacteroides sp910578895"),
    ("136_1", "Alistipes sp002428825"),
    ("138_1", "Cryptobacteroides sp009774765"),
    ("139_1", "Cryptobacteroides sp002298075"),
    ("13_1", "Lactococcus lactis"),
    ("148_1", "MGBC100320 sp910588305"),
    ("155_1", "MGBC113161 sp910579245"),
    ("3_1", "Helicobacter_C typhlonius"),
    ("4_1", "Mucispirillum sp910586745"),
    ("126_1", "CAG-873 sp009775535"),
    ("122_1", "Muribaculum intestinale"),
    ("130_1", "Bacteroides sp002491635"),
    ("135_1", "Alistipes sp009774895"),
    ("149_1", "Coproplasma sp013316055"),
    ("133_1", "CAG-485 sp002361155"),
    ("68_1", "1XD42-69 sp011959925"),
    ("123_1", "CAG-873 sp002490635"),
    ("52_1", "Kineothrix sp910588855"),
    ("124_1", "CAG-873 sp009775265"),
]


def load_species_map() -> pd.DataFrame:
    """Attach a MAG_ID to each requested species, failing loudly on any mismatch."""
    scores = pd.read_csv(FIG1_SCORES, sep="\t")
    # GTDB labels arrive with an "s__" rank prefix in some tables and not others.
    scores["sp"] = (
        scores["species"].astype(str).str.replace("s__", "", regex=False).str.strip()
    )

    want = pd.DataFrame(SPECIES, columns=["species_num", "species"])
    merged = want.merge(
        scores[["sp", "MAG_ID"]], left_on="species", right_on="sp", how="left"
    ).drop(columns="sp")

    # Both directions matter: a species with no MAG would silently vanish, and one
    # matching several MAGs would silently double-count its sites.
    missing = merged[merged.MAG_ID.isna()]
    if len(missing):
        logger.error(
            f"{len(missing)} species have no MAG in {FIG1_SCORES.name}: "
            f"{missing.species.tolist()}"
        )
        raise SystemExit(1)
    dupes = merged[merged.duplicated("species_num", keep=False)]
    if len(dupes):
        logger.error(
            f"{dupes.species.nunique()} species map to more than one MAG: "
            f"{sorted(set(dupes.species))}"
        )
        raise SystemExit(1)

    logger.info(f"matched all {len(merged)} requested species to a MAG")
    return merged


def load_significant(
    path: Path, test_type: str, group: str | None, mags: set[str]
) -> tuple[pd.DataFrame, pd.Series]:
    """Significant sites for one test plus per-MAG tested counts, for our MAGs.

    Pinning test_type AND group_analyzed is what makes a row a site; uniqueness is
    then asserted rather than imposed. Deduping instead would silently merge the two
    diet groups if the group filter were ever dropped.

    Returns (significant_sites, tested_per_mag). "Tested" = every row surviving the
    test/group/MAG pin, BEFORE the q cut -- the denominator each significance count
    should be read against. The test families test different site universes, so
    each analysis carries its own denominator; there is no shared one.
    """
    usecols = ["mag_id", "contig", "position", "test_type", "q_value"]
    if group is not None:
        usecols.append("group_analyzed")

    label = f"{test_type}{'/' + group if group else ''}"
    tested_chunks = []  # all pinned rows (keys only); sig is a q-cut view of these
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=CHUNK):
        sel = (chunk["test_type"] == test_type) & chunk["mag_id"].isin(mags)
        if group is not None:
            sel &= chunk["group_analyzed"] == group
        if sel.any():
            # Carry q_value along so significance can be cut after the dedupe check.
            tested_chunks.append(
                chunk.loc[sel, ["mag_id", "contig", "position", "q_value"]]
            )

    tested = (
        pd.concat(tested_chunks, ignore_index=True)
        if tested_chunks
        else pd.DataFrame(columns=["mag_id", "contig", "position", "q_value"])
    )
    key = ["mag_id", "contig", "position"]
    # Uniqueness over ALL tested rows, not just significant ones -- a strictly
    # stronger check, and required for the tested row count to equal a site count.
    dupes = len(tested) - len(tested.drop_duplicates(key))
    if dupes:
        logger.error(
            f"{path.name} [{label}]: {dupes:,} duplicate rows on {key}; "
            "the filters no longer pin one row per site"
        )
        raise SystemExit(1)

    tested_per_mag = tested.groupby("mag_id").size()
    out = tested.loc[tested["q_value"] < Q_THRESHOLD, key].reset_index(drop=True)
    logger.info(
        f"{path.name} [{label}]: {len(out):,} significant of {len(tested):,} tested "
        f"sites across {tested.mag_id.nunique()} of the requested MAGs"
    )
    return out, tested_per_mag


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Per species, count sites significant for BOTH divergence and high-fat "
            "parallelism, and the contigs carrying them. Inputs are pinned as module "
            "constants; the species list is SPECIES in this file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=OUTDIR,
        help="Directory to write divergence_and_hf_sites.tsv into.",
    )
    args = parser.parse_args()

    setup_logging()
    args.outdir.mkdir(parents=True, exist_ok=True)

    species = load_species_map()
    mags = set(species.MAG_ID)

    pv = RUN / "p_value_summary" / COMPARISON
    period = COMPARISON.split("-", 1)[0]
    div, div_tested = load_significant(
        pv / f"p_value_summary_two_sample_paired_{period}.tsv",
        DIVERGENCE_TEST,
        None,
        mags,
    )
    hf, hf_tested = load_significant(
        pv / f"p_value_summary_single_sample_{period}.tsv",
        PARALLELISM_TEST,
        HF_GROUP,
        mags,
    )
    # LF parallelism is loaded ONLY for its tested denominator; its significant
    # sites play no part in the divergence-AND-HF join below.
    _, lf_tested = load_significant(
        pv / f"p_value_summary_single_sample_{period}.tsv",
        PARALLELISM_TEST,
        LF_GROUP,
        mags,
    )

    # The whole point: an inner join on the SITE, so a position only survives if it
    # cleared q < 0.05 in both tests. Merging on all three key columns means a site
    # on one contig can never pair with the same coordinate on another.
    both = div.merge(hf, on=["mag_id", "contig", "position"], how="inner")
    logger.info(
        f"{len(both):,} sites significant in BOTH tests "
        f"(of {len(div):,} divergence and {len(hf):,} HF)"
    )

    per_mag = (
        both.groupby("mag_id")
        .agg(
            n_sites_both=("position", "size"),
            n_contigs_both=("contig", "nunique"),
        )
        .reset_index()
    )

    # left join, then fill: a species with no qualifying site must appear as 0, not
    # drop out of the table.
    out = species.merge(per_mag, left_on="MAG_ID", right_on="mag_id", how="left").drop(
        columns="mag_id"
    )
    # Per-analysis tested denominators (index = mag_id). A MAG absent from a
    # stratum genuinely had 0 sites tested there, so the fillna below is correct
    # for these columns too. Do NOT compare counts across columns as if they
    # shared a universe -- each test family tests its own site set.
    for col, counts in [
        ("n_sites_tested_div", div_tested),
        ("n_sites_tested_hf", hf_tested),
        ("n_sites_tested_lf", lf_tested),
    ]:
        out[col] = out["MAG_ID"].map(counts)
    count_cols = [
        "n_sites_both",
        "n_contigs_both",
        "n_sites_tested_div",
        "n_sites_tested_hf",
        "n_sites_tested_lf",
    ]
    out[count_cols] = out[count_cols].fillna(0).astype(int)

    dest = args.outdir / "divergence_and_hf_sites.tsv"
    out.to_csv(dest, sep="\t", index=False)

    # The numbers that fill the XXs in the manuscript sentence.
    logger.info(f"wrote {dest}  ({len(out)} species)")
    logger.info(
        f"  sites per species   mean {out.n_sites_both.mean():.1f}   "
        f"median {out.n_sites_both.median():.0f}   "
        f"range {out.n_sites_both.min()}-{out.n_sites_both.max()}"
    )
    logger.info(
        f"  contigs per species mean {out.n_contigs_both.mean():.1f}   "
        f"median {out.n_contigs_both.median():.0f}   "
        f"range {out.n_contigs_both.min()}-{out.n_contigs_both.max()}"
    )
    logger.info(
        f"  species with no qualifying site: {int((out.n_sites_both == 0).sum())}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
