#!/usr/bin/env python3
"""Copy the Phocaeicola sartorii assemblies out of the shared folder and write
a metadata table for labelling the ANI results."""

import logging
import os
import shutil
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(message)s")

XLSX = Path("/local1/workdir1/sidd/for_sam/Revision/Psartorii_ANI/Moeller_Lab_Bacterial_Isolate_Library.xlsx")
SRC = Path("/local1/workdir1/Shared_Folder/IsolateDataForNate/assemblies")
DEST = Path("/local1/workdir1/sidd/for_sam/Revision/Psartorii_ANI/assemblies")
META = DEST.parent / "sartorii_metadata.tsv"
MISSING_META = DEST.parent / "sartorii_missing_metadata.tsv"

DEST.mkdir(exist_ok=True)

# The "Isolate Strain Library" sheet is a superset of both "Isolate
# Classifications" sheets, so it is the single source of truth here.
df = pd.read_excel(XLSX, sheet_name="Isolate Strain Library")
hits = df[df["Isolate.Classification"].str.contains("sartorii", na=False)
          & df["Isolate.ID"].str.startswith("SLG", na=False)].copy()

# The spreadsheet records names ending in ".fna"; the shared folder stores
# them gzipped.
on_disk = set(os.listdir(SRC))
hits["Assembly.File"] = [n + ".gz" if n + ".gz" in on_disk else ""
                         for n in hits["User.Genome"]]

COLS = ["Isolate.ID", "User.Genome", "Assembly.File", "Host.Diet",
        "Host.Collection.Timepoint", "Host.Species", "Host.ID",
        "Isolate.Classification"]

hits[COLS].to_csv(META, sep="\t", index=False, na_rep="NA")

# Same columns again for any isolates with no assembly on disk, so they can be
# chased down or reported as excluded. Only written if there are any.
absent = hits[hits["Assembly.File"] == ""]
if len(absent):
    absent[COLS].to_csv(MISSING_META, sep="\t", index=False, na_rep="NA")

found = [f for f in hits["Assembly.File"] if f]
for name in found:
    shutil.copy2(SRC / name, DEST / name)

logging.info(f"{len(hits)} sartorii isolates; copied {len(found)} to {DEST}")
logging.info(f"wrote {META}")
if len(absent):
    logging.info(f"{len(absent)} have no assembly in {SRC}")
    logging.info(f"wrote {MISSING_META}")
