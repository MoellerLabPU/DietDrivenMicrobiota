#!/usr/bin/env python
"""Drop mice 534/538/539/540/541 from the Hackflex metadata sheet.

Only the Hackflex sheet needs filtering — these mice have zero rows in the
TruSeq/longitudinal metadata (verified 2026-08-19), so the anchor step is
unaffected. Run inside the alleleflux conda env.
"""

import os

import pandas as pd

SRC = "/scratch/gpfs/AMOELLER/diet_manip/copy_sg4230_scratch/sg4230/popgentoolkit/hackflex_metadata_final.txt"
# Inside the workflow's output_dir (user preference, 2026-08-20). NOTE: if you
# wipe the output tree for a clean rerun, re-run this script first — the
# filtered sheet is an input and goes down with the tree.
OUT = "/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/AlleleFlux_revision/plotting_SLG443_bin96_excl5mice/hackflex_metadata_final.excl534_538-541.tsv"
EXCLUDE = ["534", "538", "539", "540", "541"]  # str: subjectID dtype forced below

os.makedirs(os.path.dirname(OUT), exist_ok=True)
df = pd.read_csv(SRC, sep="\t", dtype={"subjectID": str})
# every excluded ID must exist in the sheet, or a typo silently under-excludes
missing = set(EXCLUDE) - set(df["subjectID"])
assert not missing, f"IDs not found in {SRC}: {sorted(missing)}"

kept = df[~df["subjectID"].isin(EXCLUDE)]
kept.to_csv(OUT, sep="\t", index=False)
print(f"{len(df)} -> {len(kept)} rows; wrote {OUT}")
