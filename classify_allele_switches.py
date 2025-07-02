#!/usr/bin/env python3
"""
Analyzes the functional impact (synonymous vs. non-synonymous) of significant
single-nucleotide variants (SNVs) where a major allele switch has occurred.

This script takes a list of significant SNV sites and, for each site, checks
for changes in the major allele between two timepoints across different replicates.
For each detected switch, it identifies the affected codon, translates it to an
amino acid, and determines if the change was synonymous or non-synonymous.

This version uses the Biopython library for robust sequence manipulation and
is parallelized to process multiple MAGs simultaneously.
"""

import argparse
import gzip
import logging
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from tqdm import tqdm


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


def parse_orf_file(orf_path: Path) -> dict:
    """
    Parses a Prodigal ORF file (.fna or .fna.gz) into a dictionary using Biopython.
    This function pre-loads all gene data into memory for a given MAG, allowing for
    very fast lookups later during the site-by-site analysis.

    Args:
        orf_path: The file path to the Prodigal-generated ORF file.

    Returns:
        A dictionary mapping each gene_id to its sequence record and metadata.
        Example:
        {
            'k141_1_24': {
                'record': SeqRecord(...), # The full Biopython object
                'contig': 'k141_1',
                'start': 24124,         # 1-based start on contig
                'end': 24882,           # 1-based end on contig
                'strand': 1             # +1 for forward, -1 for reverse
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
            contig = gene_id.rsplit("_", 1)[0]

            # Store all relevant information in the dictionary.
            gene_data[gene_id] = {
                "record": record,
                "contig": contig,
                "start": int(header_parts[1]),
                "end": int(header_parts[2]),
                "strand": int(header_parts[3]),
            }

    return gene_data


def get_major_allele(row: pd.Series, timepoint_suffix: str) -> str:
    """
    Determines the major allele for a given timepoint from a row of frequency data.
    The major allele is the one with the highest frequency.

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
    valid_freqs = {a: f for a, f in freqs.items() if pd.notna(f)}
    if not valid_freqs:
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
    them, and classifies the mutation.

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
    gene_record = gene_info["record"]
    gene_seq = gene_info["record"].seq  # The 5'-3' coding sequence from Prodigal.
    gene_start = gene_info["start"]  # 1-based start on forward contig.
    gene_end = gene_info["end"]  # 1-based end on forward contig.
    strand = gene_info["strand"]

    # --- CRITICAL LOGIC: Calculate SNV index within the gene sequence ---
    # This section handles the different coordinate systems for forward vs. reverse strands.
    # We substract 1 from gene's position to make it 0-based index. To match with the
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
        complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
        effective_allele_before = complement_map.get(allele_before)
        effective_allele_after = complement_map.get(allele_after)

        # If the allele isn't one of the four bases (e.g., 'N'), we can't analyze it.
        if effective_allele_before is None or effective_allele_after is None:
            return {}  # Return empty dict to skip this site.
    else:
        logging.warning(
            f"Unknown strand '{strand}' for gene {gene_record.id}. Skipping."
        )
        return {}

    # --- Codon Identification ---
    # Find the 0-based start of the 3-base codon that contains our mutation.
    # '//' returns the quotient (integer division) effectively rounding down.
    # 15//3 gives, 5, so does 16//3 and 17//3
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
    # '%' returns the remainder.
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
    # `cds=True` tells Biopython this is a coding sequence, which helps with error checking.
    aa_before = codon_before.translate(table=11, cds=False)
    aa_after = codon_after.translate(table=11, cds=False)

    # Classify based on whether the amino acid changed.
    mutation_type = "S" if aa_before == aa_after else "NS"

    # Return a dictionary containing all analysis results.
    return {
        "codon_before": str(codon_before),
        "codon_after": str(codon_after),
        "aa_before": str(aa_before),
        "aa_after": str(aa_after),
        "mutation_type": mutation_type,
        "strand": strand,
    }


def analyze_mag(args_tuple):
    """
    Wrapper function to analyze all significant sites for a single MAG.
    This function is designed to be called by a multiprocessing pool, where it
    receives all its necessary arguments in a single tuple.

    Args:
        args_tuple: A tuple containing (mag_id, sites_dataframe, command_line_args).
    """
    mag_id, sites, args = args_tuple
    mag_results = []

    # --- File Path and Data Loading ---
    freq_path = (
        args.frequency_dir / f"{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
    )
    orf_path_fna = args.orf_dir / f"{mag_id}.fna"
    orf_path_fnagz = args.orf_dir / f"{mag_id}.fna.gz"
    orf_path = orf_path_fna if orf_path_fna.exists() else orf_path_fnagz

    if not freq_path.exists():
        logging.warning(f"Frequency file not found for {mag_id}. Skipping.")
        return mag_results

    # OPTIMIZATION: Read only necessary columns to save memory.
    header = pd.read_csv(freq_path, sep="\t", compression="gzip", nrows=0).columns
    cols_to_use = [c for c in header if not c.endswith("_frequency_diff")]
    freq_df = pd.read_csv(freq_path, sep="\t", compression="gzip", usecols=cols_to_use)
    freq_df.set_index(["contig", "position"], inplace=True)

    # Pre-load all gene data for this MAG into a dictionary for fast access.
    gene_data = parse_orf_file(orf_path)
    if not gene_data:
        logging.warning(f"No ORF data loaded for {mag_id}. Skipping.")
        return mag_results

    # --- Timepoint Determination ---
    # Dynamically find the timepoint suffixes (e.g., 'pre', 'post') from the column names.
    timepoint_cols = [c for c in freq_df.columns if c.startswith("A_frequency_")]
    suffixes = {c.split("_")[-1] for c in timepoint_cols}
    if args.focus_timepoint not in suffixes:
        raise ValueError(
            f"Focus timepoint '{args.focus_timepoint}' not found in columns of {freq_path}. Skipping MAG."
        )

    # Find the other timepoint to use as the reference ("before").
    other_suffixes = list(suffixes - {args.focus_timepoint})
    if len(other_suffixes) != 1:
        raise ValueError(
            f"Expected one other timepoint besides '{args.focus_timepoint}', found {len(other_suffixes)}. Skipping MAG."
        )
    before_timepoint_suffix = other_suffixes[0]
    logging.info(
        f"Starting analysis for {mag_id} with before timepoint: {before_timepoint_suffix} and after timepoint: {args.focus_timepoint}"
    )
    # --- Site-by-Site Analysis ---
    # Iterate through each significant site passed to this function for this MAG.
    for _, site_row in sites.iterrows():
        contig, position, gene_id = (
            site_row["contig"],
            site_row["position"],
            site_row["gene_id"],
        )

        # Skip rows where gene_id might be missing or is not a valid string.
        if pd.isna(gene_id):
            # logging.warning(
            #     f"Skipping site with missing gene_id ({gene_id}) at contig {contig}, position {position} for mag {mag_id}."
            # )
            continue
        gene_id = str(gene_id).strip()

        # Skip rows with multiple gene_ids, which are comma-separated.
        if "," in gene_id:
            logging.warning(
                f"Skipping site with multiple gene_ids ({gene_id}) at contig {contig}, position {position} for mag {mag_id}."
            )
            continue

        # Filter the large frequency dataframe to just the rows for this specific site.
        # site_freq_data = freq_df[
        #     (freq_df["contig"] == contig) & (freq_df["position"] == position)
        # ]
        # Using [[(contig, position)]] ensures the result is always a DataFrame.
        site_freq_data = freq_df.loc[[(contig, position)]]

        # Get the pre-parsed gene information using the gene_id.
        gene_info = gene_data.get(gene_id)

        # Skip if the gene_id from the significant sites file is not found in the ORF file.
        if gene_info is None:
            raise ValueError(
                f"Gene info not found for significant site in {gene_id} at position {position} in {mag_id}.\n"
                "Make sure that the correct ORF file and significant sites file are used."
            )

        # A single site can be significant in multiple replicates and subjectIDs. Loop through each.
        for _, freq_row in site_freq_data.iterrows():
            # Determine the major allele at both timepoints.
            major_allele_before = get_major_allele(freq_row, before_timepoint_suffix)
            major_allele_after = get_major_allele(freq_row, args.focus_timepoint)

            # We only care about sites where the major allele has switched.
            if (
                major_allele_before is None
                or major_allele_after is None
                or major_allele_before == major_allele_after
            ):
                continue

            # Perform the detailed S/NS analysis using the core logic function.
            mutation_info = analyze_mutation_effect(
                gene_info, position, major_allele_before, major_allele_after
            )

            # If the analysis was successful (returned a non-empty dict), aggregate the results.
            if mutation_info:
                mag_results.append(
                    {
                        "mag_id": mag_id,
                        "subjectID": freq_row["subjectID"],
                        "replicate": freq_row["replicate"],
                        "contig": contig,
                        "position": position,
                        "gene_id": gene_id,
                        f"major_allele_before ({before_timepoint_suffix})": major_allele_before,
                        f"major_allele_after ({args.focus_timepoint})": major_allele_after,
                        **mutation_info,  # Unpack all results from the analysis function.
                    }
                )
    return mag_results


def analyze_site_for_mag(site_row: pd.Series) -> list:
    """
    Analyzes a single significant site for the currently loaded MAG.

    This function is called by a multiprocessing pool and relies on the global
    variables set by the `init_worker` function.

    Args:
        site_row: A pandas Series representing one row from the significant sites file.

    Returns:
        A list of result dictionaries, one for each replicate with an allele switch.
    """
    # These global variables are guaranteed to exist in the worker process
    # because they were set by the `init_worker` function.
    global g_freq_df, g_gene_data, g_args, g_before_timepoint_suffix, g_mag_id

    site_results = []
    contig, position, gene_id = (
        site_row["contig"],
        site_row["position"],
        site_row["gene_id"],
    )

    try:
        site_freq_data = g_freq_df.loc[[(contig, position)]]
        gene_id = str(gene_id).strip()
        gene_info = g_gene_data.get(gene_id)
    except KeyError:
        return []

    if gene_info is None:
        logging.error(
            f"Gene info for '{gene_id}' not found in ORF file for MAG {g_mag_id} "
            f"at site {contig}:{position}. Please check file consistency."
        )
        return []

    for _, freq_row in site_freq_data.iterrows():
        major_allele_before = get_major_allele(freq_row, g_before_timepoint_suffix)
        major_allele_after = get_major_allele(freq_row, g_args.focus_timepoint)

        if (
            major_allele_before is None
            or major_allele_after is None
            or major_allele_before == major_allele_after
        ):
            continue

        mutation_info = analyze_mutation_effect(
            gene_info, position, major_allele_before, major_allele_after
        )

        if mutation_info:
            site_results.append(
                {
                    "mag_id": g_mag_id,
                    "subjectID": freq_row["subjectID"],
                    "replicate": freq_row["replicate"],
                    "contig": contig,
                    "position": position,
                    "gene_id": gene_id,
                    f"major_allele_before ({g_before_timepoint_suffix})": major_allele_before,
                    f"major_allele_after ({g_args.focus_timepoint})": major_allele_after,
                    **mutation_info,
                }
            )
    return site_results


def main():
    """Main execution function to set up, run, and save the analysis."""
    # --- Setup logging and command-line arguments ---
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] - %(message)s"
    )

    parser = argparse.ArgumentParser(
        description="Analyze functional impact of significant SNVs using Biopython.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--significant_sites",
        required=True,
        type=Path,
        help="Path to the *_summary_significant_rows.tsv file.",
    )
    parser.add_argument(
        "--frequency_dir",
        required=True,
        type=Path,
        help="Directory with *_allele_frequency_changes.tsv.gz files.",
    )
    parser.add_argument(
        "--orf_dir",
        required=True,
        type=Path,
        help="Directory with Prodigal ORF files (.fna or .fna.gz).",
    )
    parser.add_argument(
        "--focus_timepoint",
        required=True,
        help="The 'after' timepoint for comparison (e.g., 'end', 'post').",
    )
    parser.add_argument(
        "--out_file", required=True, type=Path, help="Path for the output TSV file."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=cpu_count(),
        help="Number of threads to use for parallel processing.",
    )
    args = parser.parse_args()

    # --- Load initial data ---
    logging.info(f"Loading significant sites from {args.significant_sites}")
    try:
        sig_sites_df = pd.read_csv(args.significant_sites, sep="\t")
    except FileNotFoundError:
        logging.error(f"Input file not found: {args.significant_sites}. Exiting.")
        return

    # Filter the significant sites dataframe to keep only rows with a valid, single gene_id.
    logging.info(
        "Pre-filtering sites to include only those with a valid, single gene_id."
    )
    sites_to_process = sig_sites_df.dropna(subset=["gene_id"]).copy()
    sites_to_process["gene_id"] = sites_to_process["gene_id"].astype(str)
    # Exclude sites that map to multiple genes (comma-separated).
    sites_to_process = sites_to_process[~sites_to_process["gene_id"].str.contains(",")]

    # Get a sorted, unique list of MAGs to process.
    mags_to_process = sorted(sites_to_process["mag_id"].unique())
    final_results = []

    if not mags_to_process:
        logging.warning(
            "No significant sites found in any MAGs after filtering. Exiting."
        )
        return

    # --- Sequential Loop Over MAGs ---
    # The outer loop processes one MAG at a time. Parallelization happens on the sites *within* this loop.
    for mag_id in tqdm(mags_to_process, desc="Processing MAGs", unit="MAG"):
        freq_path = (
            args.frequency_dir
            / f"{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
        )
        # Use glob to find the ORF file, whether it's gzipped or not.
        orf_path = next(args.orf_dir.glob(f"{mag_id}.fna*"), None)
        # Skip this MAG if essential files are missing.
        if not freq_path.exists() or not orf_path:
            logging.warning(f"Missing frequency or ORF file for {mag_id}. Skipping.")
            continue

        # --- 1. Load data for the current MAG ---
        # Read only necessary columns to save memory.
        header = pd.read_csv(freq_path, sep="\t", compression="gzip", nrows=0).columns
        cols_to_use = [c for c in header if not c.endswith("_frequency_diff")]
        freq_df = pd.read_csv(
            freq_path, sep="\t", compression="gzip", usecols=cols_to_use
        )
        freq_df.set_index(
            ["contig", "position"], inplace=True
        )  # Set index for fast lookups.

        # Parse the entire ORF file into a dictionary for fast access.
        gene_data = parse_orf_file(orf_path)
        if not gene_data:
            logging.warning(f"No ORF data loaded for {mag_id}. Skipping.")
            continue
        # --- Timepoint Determination ---
        # Dynamically find the timepoint suffixes (e.g., 'pre', 'post') from the column names.
        timepoint_cols = [c for c in freq_df.columns if c.startswith("A_frequency_")]
        suffixes = {c.split("_")[-1] for c in timepoint_cols}
        if args.focus_timepoint not in suffixes:
            raise ValueError(
                f"Focus timepoint '{args.focus_timepoint}' not found in columns of {freq_path}. Skipping MAG."
            )

        # Find the other timepoint to use as the reference ("before").
        other_suffixes = list(suffixes - {args.focus_timepoint})
        if len(other_suffixes) != 1:
            raise ValueError(
                f"Expected one other timepoint besides '{args.focus_timepoint}', found {len(other_suffixes)}. Skipping MAG."
            )
        before_timepoint_suffix = other_suffixes[0]
        # --- 2. Prepare site-level tasks and initializer arguments ---
        # Filter the dataframe to get only the sites for the current MAG.
        mag_sites_df = sites_to_process[sites_to_process["mag_id"] == mag_id]
        # Convert the dataframe rows into a list of tasks for the pool.
        tasks = [row for _, row in mag_sites_df.iterrows()]
        if not tasks:
            logging.warning(f"No sites found for MAG {mag_id}. Skipping.")
            continue
        # Package all the data needed by the workers into a tuple.
        init_args = (freq_df, gene_data, mag_id, before_timepoint_suffix, args)
        logging.info(
            f"Starting analysis for {mag_id} with {len(mag_sites_df):,} sites across using 'before' timepoint: {before_timepoint_suffix} and 'after' timepoint: {args.focus_timepoint}"
        )
        # --- 3. Analyze all sites for this MAG in parallel ---
        # Don't use more threads than there are tasks.
        threads_to_use = min(args.threads, len(tasks))
        # Create a pool of worker processes.
        with Pool(
            processes=threads_to_use, initializer=init_worker, initargs=init_args
        ) as pool:
            # `imap_unordered` is memory-efficient. It applies the function to each task
            # and yields results as they complete, which is great for progress bars.
            for site_results_list in pool.imap_unordered(analyze_site_for_mag, tasks):
                if site_results_list:
                    final_results.extend(site_results_list)
    # --- 4. Write final output ---
    if not final_results:
        logging.warning("No major allele switches resulting in mutations were found.")
    else:
        # Convert the list of result dictionaries into a pandas DataFrame.
        results_df = pd.DataFrame(final_results)
        # Save the final DataFrame to the specified output TSV file.
        results_df.to_csv(args.out_file, sep="\t", index=False)
        logging.info(
            f"Successfully wrote {len(results_df)} mutation events to {args.out_file}"
        )


if __name__ == "__main__":
    main()
