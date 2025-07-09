#!/usr/bin/env python3
"""
Analyzes the functional impact (synonymous vs. non-synonymous) of significant
single-nucleotide variants (SNVs) where a major allele switch has occurred.

This script takes a list of significant SNV sites and, for each site, checks
for changes in the major allele between two timepoints across different replicates.
For each detected switch, it identifies the affected codon, translates it to an
amino acid, and determines if the change was synonymous or non-synonymous.

This updated version incorporates the calculation of dN/dS ratios by first
calculating the number of potential synonymous (S) and non-synonymous (N)
sites for each gene using a cached implementation of the Nei-Gojobori method.
The summary statistics are now stratified by the 'group' column to allow for
direct comparison between experimental groups.

WORKFLOW OVERVIEW:
1. Load significant sites and filter for paired t-tests
2. For each MAG:
   - Load allele frequency data and ORF sequences
   - Calculate potential S/N sites for all genes
   - Identify timepoint pairs for comparison
   - Parallelize analysis across sites
3. For each site:
   - Check for major allele switches between timepoints
   - Analyze mutation effects (synonymous vs non-synonymous)
   - Handle forward/reverse strand orientations
4. Generate summaries at position, gene, and MAG levels
5. Calculate dN/dS ratios for evolutionary analysis
"""

import argparse
import gzip
import logging
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq, Seq
from tqdm import tqdm

# --- Module-level Constants for Efficiency and Clarity ---
# A constant map for complementing DNA bases. Used for reverse-strand genes.
# When analyzing reverse-strand genes, we need to complement the alleles from
# the frequency file to match the reverse-complemented sequence from Prodigal
COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _calculate_codon_sites(
    codon: str, table: CodonTable.CodonTable
) -> tuple[float, float]:
    """
    Helper function to calculate S and N sites for a single codon.

    For a given codon, this function iterates through all 9 possible single-base
    changes and classifies each as synonymous or non-synonymous. It returns
    the fractional count of potential S and N sites for that codon.

    The Nei-Gojobori method calculates potential sites by considering all possible
    single nucleotide changes at each position and determining what fraction would
    be synonymous vs non-synonymous. Since there are 3 possible changes at each
    of the 3 positions (9 total), each position contributes 1/3 to the final count.

    Args:
        codon: A 3-base DNA string (e.g., "ATG").
        table: A Biopython CodonTable object for translation.

    Returns:
        A tuple containing the fractional number of synonymous sites and
        non-synonymous sites for the codon.
    """
    # Initialize counters for this specific codon.
    s_codon = 0.0  # Count of synonymous changes
    n_codon = 0.0  # Count of non-synonymous changes

    # Explicitly handle stop codons and invalid codons for clarity.
    if codon in table.stop_codons:
        aa_orig = "*"  # Use '*' to represent a stop codon.
    else:
        # Use .get() to handle codons that might not be in the forward table.
        aa_orig = table.forward_table.get(codon)

    # If the codon is invalid (e.g., contains 'N' or '-'), we can't analyze it.
    # We treat all 3 positions as non-synonymous sites by convention.
    if aa_orig is None:
        raise ValueError(f"Invalid codon: {codon}. Cannot analyze non-standard codons.")
        # return 0.0, 3.0

    # Iterate through each of the 3 positions in the codon (0, 1, 2).
    for i in range(3):
        original_base = codon[i]
        # Iterate through the 3 possible alternative bases at this position.
        # Note: We test all 4 bases but skip the original one
        for new_base in "ATGC":
            if new_base == original_base:
                continue  # Skip the original base, as it's not a change.

            # Create the new codon by substituting the base.
            new_codon_list = list(codon)
            new_codon_list[i] = new_base
            new_codon_str = "".join(new_codon_list)

            # Check if the new codon is a stop codon.
            if new_codon_str in table.stop_codons:
                aa_new = "*"
            else:
                aa_new = table.forward_table.get(new_codon_str)

            # If the new codon is invalid, count it as non-synonymous.
            if aa_new is None:
                raise ValueError(
                    f"New codon is invalid: {new_codon_str}. Cannot analyze non-standard codons."
                )
                # n_codon += 1
                # continue

            # Compare amino acids to classify the hypothetical mutation.
            if aa_new == aa_orig:
                s_codon += 1  # The change was synonymous (same amino acid).
            else:
                n_codon += 1  # The change was non-synonymous (different amino acid).

    # Each site contributes a fraction (s/3 or n/3) to the total.
    # This averages the potential outcomes over the 3 positions.
    # For example, if position 1 had 2 synonymous and 1 non-synonymous change,
    # position 2 had 1 synonymous and 2 non-synonymous, etc., we average them.
    return s_codon / 3.0, n_codon / 3.0


def _precompute_codon_sites_cache(table_id: int = 11) -> dict:
    """
    Pre-computes and caches the S and N site counts for all 64 possible codons.
    This avoids redundant calculations and significantly improves performance.

    Since there are only 64 possible codons (4^3), we can calculate the S/N sites
    for each one upfront and store them in a dictionary. This cache is used by
    calculate_potential_sites() to avoid recalculating the same codon multiple times.

    Args:
        table_id: The NCBI genetic code table ID (11 = bacterial/archaeal).

    Returns:
        A dictionary mapping each codon string to its (potential_S_sites, potential_N_sites) tuple.
        Example: {"ATG": (0.33, 2.67), "TTT": (0.67, 2.33), ...}
    """
    # Load the specified genetic code table.
    table = CodonTable.unambiguous_dna_by_id[table_id]
    bases = "ATGC"
    cache = {}

    # Generate all 64 possible codons programmatically (4^3 = 64).
    for a in bases:
        for b in bases:
            for c in bases:
                codon = a + b + c
                # Calculate S/N sites for each codon and store in the cache.
                # The key is the codon string, the value is the tuple of (S, N) sites.
                cache[codon] = _calculate_codon_sites(codon, table)
    return cache


# --- Pre-computed cache for all 64 codons, created once at script startup ---
_CODON_SITE_CACHE = _precompute_codon_sites_cache(table_id=11)


def calculate_potential_sites(gene_seq_str: str) -> dict:
    """
    Calculates potential synonymous (S) and non-synonymous (N) sites for a gene
    using a pre-computed cache for maximum efficiency.

    Args:
        gene_seq_str: The string of the protein-coding DNA sequence.

    Returns:
        A dictionary containing the total counts for 'N' and 'S' sites.
        Example: {"S": 45.67, "N": 123.33}
    """
    # Initialize running totals for the entire gene.
    potential_S_sites = 0.0
    potential_N_sites = 0.0

    # Iterate over the sequence in 3-base steps (codons).
    for i in range(0, len(gene_seq_str), 3):
        # Slice out the current 3-base codon.
        codon = gene_seq_str[i : i + 3].upper()

        # Skip partial codons at the end of a sequence.
        if len(codon) != 3:
            logging.warning(
                f"Gene sequence '{gene_seq_str}' of length {len(gene_seq_str):,} is not a multiple of 3. "
                f"Skipping partial codon at position {i}."
            )
            continue

        # Look up the pre-calculated S and N values from the cache. This is much
        # faster than re-calculating for every codon in every gene.
        # 1. First, get the result from the cache into a single variable.
        potential_codon_sites = _CODON_SITE_CACHE.get(codon)

        # 2. Now, check if the result is None (meaning the codon was not found).
        if potential_codon_sites is None:
            # If the codon is invalid, log a warning and skip it.
            # This handles non-standard bases (like 'Y', 'N') gracefully.
            logging.warning(
                f"Invalid codon '{codon}' found in gene '{gene_seq_str}' at position {i}. Skipping."
            )
            continue

        # 3. If the result is valid, unpack it into the two variables.
        s_codon, n_codon = potential_codon_sites

        # Add the fractional counts for this codon to the gene's total.
        potential_S_sites += s_codon
        potential_N_sites += n_codon

    # Return the final summed totals for the gene.
    return {"S": potential_S_sites, "N": potential_N_sites}


def init_worker(freq_df, gene_data, mag_id, before_timepoint_suffix, cli_args):
    """
    Initializer for each worker process in the pool.

    This function is called once per worker. It sets up global variables
    within the worker's own memory space. This is the most efficient way to
    provide large, read-only objects to worker processes, as it avoids
    repeatedly sending the data with each task.

    Args:
        freq_df (pd.DataFrame): The allele frequency data for the current MAG.
        gene_data (dict): The parsed ORF/gene data for the current MAG.
        mag_id (str): The ID of the MAG being processed.
        before_timepoint_suffix (str): The suffix for the 'before' timepoint.
        cli_args (argparse.Namespace): The parsed command-line arguments.
    """
    # These globals are specific to each worker process, not shared between them.
    global g_freq_df, g_gene_data, g_mag_id, g_before_timepoint_suffix, g_args
    g_freq_df = freq_df
    g_gene_data = gene_data
    g_mag_id = mag_id
    g_before_timepoint_suffix = before_timepoint_suffix
    g_args = cli_args


def setup_and_load_data(args):
    """
    Loads and filters the initial significant sites data.

    This function reads the significant sites file and applies filtering to ensure
    we only analyze sites that meet our criteria (e.g., paired t-tests, single gene assignments).

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Returns:
        A pandas DataFrame of sites to be processed, or None if no sites are found.
    """
    logging.info(f"Loading significant sites from {args.significant_sites}")
    sig_sites_df = pd.read_csv(args.significant_sites, sep="\t")

    # Filter for paired t-tests if the column exists
    if "test_type" in sig_sites_df.columns:
        logging.info("Filtering significant sites for test_type == 'paired_tTest'")
        sig_sites_df = sig_sites_df[sig_sites_df["test_type"] == "paired_tTest"].copy()
        if sig_sites_df.empty:
            logging.warning("No sites found with test_type 'paired_tTest'. Exiting.")
            return None
    else:
        logging.warning(
            "Column 'test_type' not found in significant_sites file. Proceeding without filtering."
        )

    logging.info(
        "Pre-filtering sites to include only those with a valid, single gene_id."
    )
    sites_to_process = sig_sites_df.dropna(subset=["gene_id"]).copy()
    # Exclude sites that map to multiple genes (comma-separated)
    sites_to_process = sites_to_process[
        ~sites_to_process["gene_id"].astype(str).str.contains(",")
    ]

    return sites_to_process


def process_single_mag(mag_id, sites_for_mag, args):
    """
    Handles all processing for a single MAG: loading files, calculating
    potential sites, and running the parallel analysis.

    Args:
        mag_id (str): The identifier for the MAG to process.
        sites_for_mag (pd.DataFrame): A DataFrame of significant sites for this MAG.
        args (argparse.Namespace): The parsed command-line arguments.

    Returns:
        A list of result dictionaries for all observed mutation events in the MAG.
    """
    # --- 1. Load MAG-specific data ---
    freq_path = (
        args.frequency_dir / f"{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
    )
    # Use glob to find the ORF file, whether it's gzipped or not.
    orf_path = next(args.orf_dir.glob(f"{mag_id}.fna*"), None)
    # Skip this MAG if essential files are missing.
    if not freq_path.exists() or not orf_path:
        logging.warning(f"Missing frequency or ORF file for {mag_id}. Skipping.")
        return []

    # Read only necessary columns to save memory.
    # We exclude frequency difference columns as they're not needed for this analysis.
    header = pd.read_csv(freq_path, sep="\t", compression="gzip", nrows=0).columns
    cols_to_use = [c for c in header if not c.endswith("_frequency_diff")]
    freq_df = pd.read_csv(freq_path, sep="\t", compression="gzip", usecols=cols_to_use)
    # Set multi-index for fast lookups by contig and position
    freq_df.set_index(["contig", "position"], inplace=True)

    # Parse the entire ORF file into a dictionary for fast access.
    gene_data = parse_orf_file(orf_path)
    if not gene_data:
        logging.warning(f"No ORF data loaded for {mag_id}. Skipping.")
        return []

    # --- 2. Pre-calculate N and S sites for all genes in the MAG ---
    # This step calculates the denominator for dN/dS calculations.
    # We do this once per gene rather than once per site to avoid redundancy.
    logging.info(
        f"Calculating potential synonymous and nonsynonymous sites for {len(gene_data):,} genes in {mag_id}..."
    )
    for gene_id, info in gene_data.items():
        seq_str = str(info["record"].seq)
        if len(seq_str) % 3 != 0:
            logging.warning(
                f"Gene {gene_id} in {mag_id} has length {len(seq_str):,} (not a multiple of 3)."
            )
        # Update the gene data dict with potential S and N sites.
        # The enhanced gene_data is later passed to the workers.
        potential_n_s_counts = calculate_potential_sites(seq_str)
        info["potential_S_sites"] = potential_n_s_counts["S"]
        info["potential_N_sites"] = potential_n_s_counts["N"]

    # --- 3. Determine timepoints and prepare for parallel analysis ---
    # Dynamically find the timepoint suffixes (e.g., 'pre', 'post') from the column names.
    timepoint_cols = [c for c in freq_df.columns if c.startswith("A_frequency_")]
    suffixes = {c.split("_")[-1] for c in timepoint_cols}

    # Validate that our focus timepoint exists in the data
    if args.focus_timepoint not in suffixes:
        logging.warning(
            f"Focus timepoint '{args.focus_timepoint}' not in {freq_path}. Skipping MAG."
        )
        return []

    # Find the other timepoint to use as the reference ("before").
    other_suffixes = list(suffixes - {args.focus_timepoint})
    if len(other_suffixes) != 1:
        logging.warning(
            f"Expected one other timepoint, found {other_suffixes}. Skipping MAG."
        )
        return []
    before_timepoint_suffix = other_suffixes[0]

    # Convert the dataframe rows into a list of tasks for the pool.
    # Each task represents one significant site to analyze.
    tasks = [row for _, row in sites_for_mag.iterrows()]

    # Package all the data needed by the workers into a tuple.
    init_args = (freq_df, gene_data, mag_id, before_timepoint_suffix, args)
    logging.info(
        f"\nStarting analysis for {mag_id} with {len(tasks):,} sites. "
        f"Before timepoint: '{before_timepoint_suffix}', After timepoint: '{args.focus_timepoint}'"
    )

    # --- 4. Run analysis in parallel ---
    mag_results = []
    # Don't use more threads than there are tasks to avoid overhead
    threads_to_use = min(args.threads, len(tasks))
    with Pool(
        processes=threads_to_use, initializer=init_worker, initargs=init_args
    ) as pool:
        # `imap_unordered` is memory-efficient. It applies the function to each task
        # and yields results as they complete, which is great for progress bars.
        for site_results_list in pool.imap_unordered(analyze_site_for_mag, tasks):
            if site_results_list:
                mag_results.extend(site_results_list)
    return mag_results


def parse_orf_file(orf_path: Path) -> dict:
    """
    Parses a Prodigal ORF file (.fna or .fna.gz) into a dictionary using Biopython.
    This function pre-loads all gene data into memory for a given MAG, allowing for
    very fast lookups later during the site-by-site analysis.

    Prodigal outputs FASTA files where each gene has a header like:
    >k141_1_24 # 24124 # 24882 # 1 # ID=1_24;partial=00;start_type=ATG;...

    The header contains: gene_id # start # end # strand # metadata
    - start/end are 1-based coordinates on the contig
    - strand is +1 (forward) or -1 (reverse)
    - The sequence is always provided in the 5'->3' coding direction

    Args:
        orf_path: The file path to the Prodigal-generated ORF file.

    Returns:
        A dictionary mapping each gene_id to its sequence record and metadata.
        Example:
        {
            'k141_1_24': {
                'record': SeqRecord(...), # The full Biopython object
                'contig': 'k141_1',       # Contig name (derived from gene_id)
                'start': 24124,           # 1-based start on contig
                'end': 24882,             # 1-based end on contig
                'strand': 1               # +1 for forward, -1 for reverse
            },
            ...
        }
    """
    gene_data = {}
    # If the file doesn't exist (e.g., Prodigal failed for a MAG), return an empty dict.
    if not orf_path.exists():
        logging.warning(f"ORF file not found: {orf_path}")
        return gene_data

    # Handle both gzipped (.gz) and regular text files automatically.
    _open = gzip.open if orf_path.suffix == ".gz" else open

    # Use Biopython's optimized FASTA parser.
    with _open(orf_path, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            # The gene_id is the unique identifier (e.g., 'contig_1_1').
            gene_id = record.id
            # The rest of the header contains metadata separated by '#'.
            # Example: >contig_1_1 # 123 # 456 # 1 # ...
            header_parts = [p.strip() for p in record.description.split("#")]
            # Reconstruct the contig ID by removing the ORF number suffix.
            # For gene_id 'k141_1_24', contig would be 'k141_1'
            contig = gene_id.rsplit("_", 1)[0]

            # Store all relevant information in the dictionary.
            gene_data[gene_id] = {
                "record": record,  # Full Biopython SeqRecord
                "contig": contig,  # Parent contig name
                "start": int(header_parts[1]),  # 1-based start position
                "end": int(header_parts[2]),  # 1-based end position
                "strand": int(header_parts[3]),  # Strand (+1 or -1)
            }

    return gene_data


def get_major_allele(row: pd.Series, timepoint_suffix: str) -> str:
    """
    Determines the major allele for a given timepoint from a row of frequency data.
    The major allele is the one with the highest frequency among A, T, G, C.

    This function looks for columns named like "A_frequency_pre", "T_frequency_post", etc.
    and identifies which allele has the highest frequency at the specified timepoint.

    Args:
        row: A pandas Series representing one row from the allele frequency file.
        timepoint_suffix: The string identifying the timepoint (e.g., 'pre', 'post').

    Returns:
        The major allele ('A', 'T', 'G', or 'C') or None if there is a tie,
        no alleles are found, or frequencies are not valid numbers.
    """
    alleles = ["A", "T", "G", "C"]
    # Create a dictionary of allele frequencies for the specified timepoint.
    # e.g., {'A': 0.9, 'T': 0.1, 'G': 0.0, 'C': 0.0} for suffix 'post'.
    freqs = {
        allele: row.get(f"{allele}_frequency_{timepoint_suffix}", np.nan)
        for allele in alleles
    }

    # Filter out any alleles with non-numeric or missing frequencies.
    # This handles cases where some allele frequencies might be NaN or missing.
    valid_freqs = {a: f for a, f in freqs.items() if pd.notna(f)}
    if not valid_freqs:
        logging.warning(
            f"No valid allele frequencies found for timepoint '{timepoint_suffix}' in row: {row.name}"
        )
        return None

    # Find the highest frequency value.
    max_freq = max(valid_freqs.values())
    # Find all alleles that have this maximum frequency.
    major_alleles = [allele for allele, freq in valid_freqs.items() if freq == max_freq]

    # Return the allele only if it's unique to avoid ambiguity from ties.
    if len(major_alleles) == 1:
        return major_alleles[0]

    # If there's a tie (e.g., A=0.5, T=0.5), there is no single major allele.
    return None


def analyze_mutation_effect(
    gene_info: dict, position: int, allele_before: str, allele_after: str
) -> dict:
    """
    Determines the functional effect of a mutation by handling coordinate systems
    and allele orientations for both forward and reverse strand genes.

    This is the core logic function. It translates the SNV's contig-based position
    into a gene-based index, builds the 'before' and 'after' codons, translates
    them, and classifies the mutation as synonymous (S) or non-synonymous (NS).

    COORDINATE SYSTEM EXPLANATION:
    - Input position: 0-based position on the contig (from frequency file)
    - Prodigal coordinates: 1-based start/end positions on the contig
    - Gene sequence: Always 5'->3' coding sequence, regardless of strand

    For FORWARD genes: position in gene = SNV_pos - (gene_start - 1)
    For REVERSE genes: position in gene = (gene_end - 1) - SNV_pos

    STRAND HANDLING:
    - Forward strand: Use alleles directly from frequency file
    - Reverse strand: Complement alleles to match reverse-complement sequence

    Args:
        gene_info: The dictionary for one gene from parse_orf_file().
        position: The 0-based position of the SNV on the contig.
        allele_before: The major allele at the first timepoint (e.g., 'A').
        allele_after: The major allele at the second timepoint (e.g., 'G').

    Returns:
        A dictionary with the analysis results, or an empty dictionary if the
        analysis cannot be completed (e.g., SNV in partial codon).
    """
    # --- Input Standardization ---
    # Ensure all inputs are uppercase to prevent case-sensitivity issues. This handles
    # potential lowercase alleles from the freq file so they are not skipped.
    allele_before = allele_before.upper()
    allele_after = allele_after.upper()

    # Unpack gene information for clarity.
    gene_record = gene_info["record"]  # Biopython SeqRecord object
    gene_seq = gene_info["record"].seq  # The 5'-3' coding sequence from Prodigal
    gene_start = gene_info["start"]  # 1-based start on forward contig
    gene_end = gene_info["end"]  # 1-based end on forward contig
    strand = gene_info["strand"]  # +1 or -1

    # --- CRITICAL LOGIC: Calculate SNV index within the gene sequence ---
    # This section handles the different coordinate systems for forward vs. reverse strands.
    # We subtract 1 from gene's position to make it 0-based index to match with the
    # index of the positions in the freq file.

    # EXAMPLE:
    #   - Contig is 10000 bp long.
    #   - SNV Position (0-based from freq file): 8025

    #   - FORWARD GENE EXAMPLE:
    #     - Prodigal Header: # 8000 # 8500 # 1
    #     - The gene sequence starts at contig position 8000.
    #     - Calculation: 8025 (SNV) - (8000 - 1) (gene start) = index 26
    #     - This is intuitive: the SNV is 26 bases into the gene.

    #   - REVERSE GENE EXAMPLE:
    #     - Prodigal Header: # 8000 # 8500 # -1
    #     - The provided gene sequence is the REVERSE COMPLEMENT. Its first base
    #       corresponds to contig position 8500.
    #     - Calculation: (8500 - 1) (gene end) - 8025 (SNV) = index 474
    #     - This is counter-intuitive but correct: the SNV is near the 'start'
    #       coordinate on the contig, which means it's at the *end* of the
    #       reverse-complemented sequence string.

    if strand == 1:
        # --- FORWARD STRAND ---
        # The position is a direct offset from the gene's start.
        pos_in_gene = position - (gene_start - 1)
        # Alleles from the frequency file are on the forward strand, so use them directly.
        effective_allele_before = allele_before
        effective_allele_after = allele_after

    elif strand == -1:
        # --- REVERSE STRAND ---
        # The position is an offset from the gene's end (which is the start of the rev-comp seq).
        pos_in_gene = (gene_end - 1) - position

        # The alleles ('A', 'T', 'G', 'C') were detected on the FORWARD strand.
        # We must complement them to match the reverse-complemented `gene_seq`.
        # Example: An 'A' on the forward strand is a 'T' on the reverse strand.
        effective_allele_before = COMPLEMENT_MAP.get(allele_before)
        effective_allele_after = COMPLEMENT_MAP.get(allele_after)

        # If the allele isn't one of the four bases (e.g., 'N'), we can't analyze it.
        if effective_allele_before is None or effective_allele_after is None:
            logging.warning(
                f"Allele '{allele_before}' or '{allele_after}' is not a valid base for gene {gene_record.id} "
                f"on reverse strand. Skipping."
            )
            return {}  # Return empty dict to skip this site.
    else:
        logging.warning(
            f"Unknown strand '{strand}' for gene {gene_record.id}. Skipping."
        )
        return {}

    # --- Codon Identification ---
    # Find the 0-based start of the 3-base codon that contains our mutation.
    # '//' returns the quotient (integer division) effectively rounding down.
    # Examples: 15//3 = 5, 16//3 = 5, 17//3 = 5 (all in codon 5)
    codon_start_in_gene = (pos_in_gene // 3) * 3
    # Get the codon from the reference ORF sequence.
    codon_context = gene_seq[codon_start_in_gene : codon_start_in_gene + 3].upper()

    # Edge case: If the SNV is in a partial codon at the very end of a gene, we can't analyze it.
    if len(codon_context) != 3:
        logging.warning(
            f"SNV at pos {position} in gene {gene_record.id} is in a partial codon. Skipping."
        )
        return {}

    # Find the position within the 3-base codon (0, 1, or 2).
    # '%' returns the remainder. Examples: 15%3 = 0, 16%3 = 1, 17%3 = 2
    pos_in_codon = pos_in_gene % 3

    # --- Construct 'Before' and 'After' Codons ---
    # Create mutable versions of the codon to substitute the alleles.
    mutable_codon_before = MutableSeq(codon_context)
    mutable_codon_before[pos_in_codon] = effective_allele_before
    codon_before = Seq(mutable_codon_before)

    mutable_codon_after = MutableSeq(codon_context)
    mutable_codon_after[pos_in_codon] = effective_allele_after
    codon_after = Seq(mutable_codon_after)

    # --- Translation and Classification ---
    # Prodigal's .fna provides the correct 5'-to-3' coding sequence for ALL genes.
    # Therefore, we can translate directly without further manipulation.
    # `table=11` is the standard bacterial/archaeal genetic code.
    # `cds=False` allows translation of sequences that might not start with start codons
    aa_before = codon_before.translate(table=11, cds=False)
    aa_after = codon_after.translate(table=11, cds=False)

    # Classify based on whether the amino acid changed.
    # S = Synonymous (same amino acid), NS = Non-synonymous (different amino acid)
    mutation_type = "S" if aa_before == aa_after else "NS"

    # Return a dictionary containing all analysis results.
    return {
        "pos_in_gene": pos_in_gene,
        "pos_in_codon (0-based index)": pos_in_codon,
        "codon_before": str(codon_before),
        "codon_after": str(codon_after),
        "aa_before": str(aa_before),
        "aa_after": str(aa_after),
        "mutation_type": mutation_type,
        "strand": strand,
    }


def analyze_site_for_mag(site_row: pd.Series) -> list:
    """
    Analyzes a single significant site for the currently loaded MAG.

    This function is called by a multiprocessing pool and relies on the global
    variables set by the `init_worker` function. For each significant site,
    it checks if there are major allele switches between timepoints across
    different replicates/subjects, and analyzes the functional impact.

    WORKFLOW:
    1. Extract site coordinates and gene info
    2. Look up frequency data for this site
    3. For each replicate/subject at this site:
       - Determine major alleles at both timepoints
       - Check if there's a clear allele switch
       - Analyze mutation effect (synonymous vs non-synonymous)
    4. Return list of results (one per replicate with allele switch)

    Args:
        site_row: A pandas Series representing one row from the significant sites file.

    Returns:
        A list of result dictionaries, one for each replicate with an allele switch.
    """
    # These global variables are guaranteed to exist in the worker process
    # because they were set by the `init_worker` function.
    global g_freq_df, g_gene_data, g_args, g_before_timepoint_suffix, g_mag_id

    site_results = []
    # Extract site information from the input row.
    contig, position, gene_id = (
        site_row["contig"],
        site_row["position"],
        site_row["gene_id"],
    )

    # Look up frequency data for this specific site (contig and position).
    # Using .loc[[(...)]] ensures the result is always a DataFrame.
    site_freq_data = g_freq_df.loc[[(contig, position)]]
    gene_id = str(gene_id).strip()
    # Get the pre-parsed gene information using the gene_id.
    gene_info = g_gene_data.get(gene_id)

    # Check if the gene_id from the significant sites file was found in the ORF file.
    if gene_info is None:
        logging.error(
            f"Gene info for '{gene_id}' not found in ORF file for MAG {g_mag_id} "
            f"at site {contig}:{position}. Please check file consistency."
        )
        return []

    # A single site can be significant in multiple subjectID. Loop through each one.
    for _, freq_row in site_freq_data.iterrows():
        # Determine the major allele at both timepoints for this subjectID.
        major_allele_before = get_major_allele(freq_row, g_before_timepoint_suffix)
        major_allele_after = get_major_allele(freq_row, g_args.focus_timepoint)

        # We only care about sites where the major allele has clearly switched.
        if (
            major_allele_before is None
            or major_allele_after is None
            or major_allele_before == major_allele_after
        ):
            continue

        # Perform the detailed S/NS analysis.
        mutation_info = analyze_mutation_effect(
            gene_info, position, major_allele_before, major_allele_after
        )

        # If the analysis was successful, aggregate the results.
        if mutation_info:
            site_results.append(
                {
                    "mag_id": g_mag_id,
                    "subjectID": freq_row["subjectID"],
                    "replicate": freq_row["replicate"],
                    "group": freq_row["group"],
                    "contig": contig,
                    "position": position,
                    "gene_id": gene_id,
                    "major_allele_before": major_allele_before,
                    "major_allele_after": major_allele_after,
                    "potential_N_sites": gene_info.get("potential_N_sites", np.nan),
                    "potential_S_sites": gene_info.get("potential_S_sites", np.nan),
                    **mutation_info,  # Unpack all results from the analysis function.
                }
            )
    return site_results


def summarize_results(results_df):
    """
    Generates several summary tables from the detailed mutation event data.

    This function creates three levels of summary:
    1. Position-level: Aggregates data by site location (contig + position)
    2. Gene-level: Aggregates data by gene and calculates dN/dS ratios
    3. MAG-level: Aggregates data by MAG and calculates overall dN/dS ratios

    All summaries are stratified by experimental group to enable comparisons.

    Args:
        results_df: DataFrame containing the full, detailed output from all mutation events.

    Returns:
        A dictionary of DataFrames, where each key is a summary type
        (e.g., 'position', 'gene', 'mag') and the value is the summary DataFrame.
    """

    # --- Helper function to show all unique values in a group ---
    def summarize_variants(series):
        """Joins unique values with a ',' to show all variants at a site."""
        # Use a comma as a separator.
        return ",".join(series.astype(str).unique())

    # --- 1. Position-level Summary (Stratified by Group) ---
    # Group by the site of the mutation and the experimental group.
    # This tells us how many subjects had switches at each position in each group.
    grouped = results_df.groupby(["mag_id", "group", "contig", "position", "gene_id"])

    # Perform aggregations to get detailed stats for each site
    position_stats = grouped.agg(
        num_subjects_changed=("subjectID", "nunique"),
        # Count "S" and "NS" instead of full words.
        s_count=("mutation_type", lambda s: (s == "S").sum()),
        ns_count=("mutation_type", lambda s: (s == "NS").sum()),
        major_allele_before=("major_allele_before", summarize_variants),
        major_allele_after=("major_allele_after", summarize_variants),
        codon_before=("codon_before", summarize_variants),
        codon_after=("codon_after", summarize_variants),
        aa_before=("aa_before", summarize_variants),
        aa_after=("aa_after", summarize_variants),
    ).reset_index()

    # Calculate the percentage of subjects with synonymous vs. non-synonymous changes.
    position_stats["percent_synonymous"] = (
        (position_stats["s_count"] / position_stats["num_subjects_changed"]) * 100
    ).round(2)
    position_stats["percent_nonsynonymous"] = (
        (position_stats["ns_count"] / position_stats["num_subjects_changed"]) * 100
    ).round(2)

    # --- 2. Gene-level Statistics (Stratified by Group, with dN/dS) ---
    gene_mutation_counts = (
        results_df.groupby(["mag_id", "group", "gene_id"])
        .agg(
            total_mutations=("mutation_type", "size"),  # Total number of mutations
            s_count=("mutation_type", lambda s: (s == "S").sum()),
            ns_count=("mutation_type", lambda s: (s == "NS").sum()),
        )
        .reset_index()
    )

    # Get the potential synonymous and nonsynonymous site counts for each gene.
    # These counts are group-independent (they depend only on the gene sequence).
    gene_site_counts = results_df[
        ["mag_id", "gene_id", "potential_N_sites", "potential_S_sites"]
    ].drop_duplicates()
    # Merge the site counts, broadcasting them to each group-specific row for a gene.
    gene_stats = pd.merge(
        gene_mutation_counts, gene_site_counts, on=["mag_id", "gene_id"]
    )

    # Calculate dN and dS for each group within each gene.
    # dN = (# non-synonymous mutations) / (# potential non-synonymous sites)
    # dS = (# synonymous mutations) / (# potential synonymous sites)
    gene_stats["dN"] = gene_stats.apply(
        lambda row: (
            row["ns_count"] / row["potential_N_sites"]
            if row["potential_N_sites"] > 0
            else 0
        ),
        axis=1,
    )
    gene_stats["dS"] = gene_stats.apply(
        lambda row: (
            row["s_count"] / row["potential_S_sites"]
            if row["potential_S_sites"] > 0
            else 0
        ),
        axis=1,
    )
    # Safer dN/dS calculation: return NaN if dS is 0 to avoid division errors.
    gene_stats["dN_dS_ratio"] = gene_stats.apply(
        lambda row: row["dN"] / row["dS"] if row["dS"] > 0 else np.nan, axis=1
    )

    # --- 3. MAG-level Statistics (Stratified by Group, with dN/dS) ---
    mag_stats = (
        gene_stats.groupby(["mag_id", "group"])
        .agg(
            total_ns_count=("ns_count", "sum"),
            total_s_count=("s_count", "sum"),
            potential_N_sites=("potential_N_sites", "sum"),
            potential_S_sites=("potential_S_sites", "sum"),
        )
        .reset_index()
    )

    # Calculate overall dN, dS, and dN/dS for each group within each MAG.
    mag_stats["dN"] = mag_stats.apply(
        lambda row: (
            row["total_ns_count"] / row["potential_N_sites"]
            if row["potential_N_sites"] > 0
            else 0
        ),
        axis=1,
    )
    mag_stats["dS"] = mag_stats.apply(
        lambda row: (
            row["total_s_count"] / row["potential_S_sites"]
            if row["potential_S_sites"] > 0
            else 0
        ),
        axis=1,
    )
    mag_stats["dN_dS_ratio"] = mag_stats.apply(
        lambda row: row["dN"] / row["dS"] if row["dS"] > 0 else np.nan, axis=1
    )

    return {"position": position_stats, "gene": gene_stats, "mag": mag_stats}


def write_output_files(summaries, results_df, args):
    """
    Writes the main events table and all summary tables to disk.

    This function handles the final output step, writing:
    1. All individual mutation events to a detailed TSV file
    2. Position-level summary statistics
    3. Gene-level summary statistics (including dN/dS ratios)
    4. MAG-level summary statistics (including dN/dS ratios)

    Args:
        summaries (dict): A dictionary of summary DataFrames from summarize_results().
        results_df (pd.DataFrame): The DataFrame of all mutation events.
        args (argparse.Namespace): The parsed command-line arguments.
    """
    # Write the detailed events file if we have results
    if not results_df.empty:
        path_all_events = args.outdir / f"{args.prefix}_all_events.tsv"
        # Save the final DataFrame to the specified output TSV file.
        results_df.to_csv(path_all_events, sep="\t", index=False)
        logging.info(
            f"Successfully wrote {len(results_df)} mutation events to {path_all_events}"
        )

    # Check if we have summaries to write
    if not summaries:
        logging.warning("Could not generate summaries.")
        return

    # Write each summary type to its own file
    logging.info("Writing summary files...")
    for summary_type, df in summaries.items():
        if df is not None and not df.empty:
            # Replace infinite values with NaN for cleaner output
            # This can happen if dS = 0 leading to infinite dN/dS ratios
            # df.replace([np.inf, -np.inf], np.nan, inplace=True)
            path = args.outdir / f"{args.prefix}_{summary_type}_summary.tsv"
            df.to_csv(path, sep="\t", index=False)
            logging.info(f"Wrote {summary_type} summary to {path}")


def main():
    """
    Main execution function to set up, run, and save the analysis.

    OVERALL WORKFLOW:
    1. Parse command-line arguments and set up logging
    2. Load and filter significant sites data
    3. For each MAG with significant sites:
       - Load frequency and ORF data
       - Calculate potential S/N sites for all genes
       - Run parallel analysis across sites
    4. Aggregate results and generate summaries
    5. Write output files with detailed events and summaries
    """
    # --- Setup logging and command-line arguments ---
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="Analyze functional impact of significant SNVs using Biopython.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--significant_sites",
        required=True,
        type=Path,
        help="Path to the *_summary_significant_rows.tsv file containing significant SNV sites.",
    )
    parser.add_argument(
        "--frequency_dir",
        required=True,
        type=Path,
        help="Directory with *_allele_frequency_changes.tsv.gz files containing allele frequencies.",
    )
    parser.add_argument(
        "--orf_dir",
        required=True,
        type=Path,
        help="Directory with Prodigal ORF files (.fna or .fna.gz) containing gene sequences.",
    )
    parser.add_argument(
        "--focus_timepoint",
        required=True,
        help="The 'after' timepoint for comparison (e.g., 'end', 'post'). Must match suffix in frequency files.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory to save the output files. Will be created if it doesn't exist.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="allele_switch_summary",
        help="Prefix for the output file names.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=cpu_count(),
        help="Number of threads to use for parallel processing.",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # --- 1. Setup and Initial Data Loading ---
    # Load the significant sites file and apply initial filtering
    sites_to_process = setup_and_load_data(args)
    if sites_to_process is None or sites_to_process.empty:
        logging.warning("No significant sites to process after filtering. Exiting.")
        return

    logging.info(
        f"Processing {len(sites_to_process):,} significant sites across {sites_to_process['mag_id'].nunique()} MAGs"
    )

    # --- 2. Main Processing Loop ---
    # Process each MAG independently and collect results
    all_results = []
    mags_to_process = sorted(sites_to_process["mag_id"].unique())

    for mag_id in tqdm(
        mags_to_process, desc="Processing MAGs", unit="MAG", total=len(mags_to_process)
    ):
        # Filter the dataframe to get only the sites for the current MAG.
        sites_for_mag = sites_to_process[sites_to_process["mag_id"] == mag_id]
        if sites_for_mag.empty:
            logging.warning(f"No sites found for MAG {mag_id}. Skipping.")
            continue

        # Process this MAG and collect mutation events
        mag_results = process_single_mag(mag_id, sites_for_mag, args)
        if mag_results:
            all_results.extend(mag_results)

    # --- 3. Final Aggregation and Output ---
    if not all_results:
        logging.warning(
            "No major allele switches resulting in mutations were found across all MAGs."
        )
        return

    # Convert the list of result dictionaries into a pandas DataFrame.
    results_df = pd.DataFrame(all_results)

    # Generate summary statistics at different levels
    summaries = summarize_results(results_df)

    # Write all output files
    write_output_files(summaries, results_df, args)

    logging.info("Analysis completed successfully!")


if __name__ == "__main__":
    main()
