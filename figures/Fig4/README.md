# Fig 4 — Phase Variation Analysis

Figure 4 presents the phase variation analysis of *Bacteroides muris* isolates using [PhaseFinder](https://github.com/XiaofangJ/PhaseFinder).

## Code

- **`Snakefile_miceisolate_phasevariation`** — Snakemake workflow that runs PhaseFinder to:
  1. Locate invertible DNA regions in the reference genome
  2. Create synthetic inversions for read mapping
  3. Calculate inversion ratios from paired-end isolate reads

- **`filter_out_file.py`** — Filters PhaseFinder ratio output files by Pe_ratio, Pe_R/Pe_F, and Span_R/Span_F thresholds.

- **`combine.py`** — Combines filtered PhaseFinder output files across samples, extracting MAG and SLG identifiers from filenames.

## Usage

```bash
# Run PhaseFinder pipeline
snakemake -s Snakefile_miceisolate_phasevariation --cores <N>

# Filter results
python filter_out_file.py /path/to/ratio_files

# Combine filtered results
python combine.py /path/to/filtered_files
```
