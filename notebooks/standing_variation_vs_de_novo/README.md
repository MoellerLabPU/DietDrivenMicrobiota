# Standing variation or de novo mutation?

For the sites AlleleFlux called significant under the high-fat diet, was the allele
that rose already present in the mice at PRE (standing variation), or did it appear
afterwards (de novo candidate)? Three steps, each feeding the next. All read the
published `AlleleFlux_mapq20` run.

## 1. `div_and_hf_sites_in_both.py`

For each species that did **not** show strain replacement in more than half of the
HF replicates, counts the sites with q < 0.05 in **both** the divergence test
(`two_sample_paired_tTest`) and HF parallelism (`single_sample_tTest`, `fat`), and
the contigs carrying at least one. A site qualifies only if the same (contig,
position) clears the threshold in both tests, an inner join on position rather
than a comparison of per-contig totals. Writes one TSV; `--help` lists the paths.

## 2. `pre_allele_presence.ipynb`

Two species from step 1 carry only a handful of such sites: `SLG1171_DASTool_bins_9`
(*Odoribacter* sp910578105, 2 sites) and `SLG421_DASTool_bins_82`
(*Cryptobacteroides* sp009774765, 3 sites). The notebook re-reads the raw
`profiles/` for those sites, calls each mouse's major allele, builds a consensus
per diet × timepoint, and asks whether the HF END allele was already present in
the HF mice at PRE, in how many mice and at what depth. Section 4 repeats the
question for the divergent allele rather than the major allele; section 5 draws
per-cage allele-frequency trajectories as an eye test.

## 3. `baseline_presence/`

The same question for **every** significant site, via the packaged
`alleleflux-baseline-presence` command. `run_baseline_presence.sbatch` is the
SLURM submission as run (`pre_end-fat_control`, `two_sample_paired_tTest`,
q < 0.05, minimum 5 reads and 5% frequency). Its README documents the flow, the
filters in the order they apply, and both output tables.
