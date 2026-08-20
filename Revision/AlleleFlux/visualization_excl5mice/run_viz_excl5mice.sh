#!/bin/bash
#SBATCH --job-name=viz_excl5mice
#SBATCH --output=logs/viz_excl5mice_%j.out
#SBATCH --error=logs/viz_excl5mice_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1

# Re-run the full AlleleFlux visualization workflow (metadata prep -> terminal
# nucleotides -> frequency tracking -> line plots at bin widths 10/20/30/40/50,
# each output tagged _bin<W>d) with mice 534/538/539/540/541 excluded at the
# metadata level. Multi-width plotting needs the AlleleFlux-dev checkout with
# the 2026-08-19 multi-bin support (nargs list on --bin_width_days).
#
# PREREQUISITE (run once, interactively, before submitting):
#   python make_filtered_hackflex_metadata.py
# The config points at the filtered sheet this writes; snakemake will fail
# fast on a missing input if the prerequisite was skipped.
#
# Pattern mirrors the July 2026 per-replicate run: raw snakemake against the
# viz Snakefile in AlleleFlux-dev, executed locally inside one SLURM job
# (16 cpus / 150G per the config's resources block), from a dedicated workdir
# on scratch so .snakemake state never collides with other runs.

SNAKEFILE=/home/su2806/AlleleFlux-dev/alleleflux/smk_workflow/visualization/visualization.smk
CONFIG=/home/su2806/DietDrivenMicrobiota/Revision/AlleleFlux/visualization_excl5mice/alleleflux_visualization_config_excl5mice.yaml

args=(
    -s "$SNAKEFILE"
    --configfile "$CONFIG"     # overrides the template configfile in the Snakefile
    --profile /home/su2806/AlleleFlux-dev/alleleflux/smk_workflow/slurm_profile_native
)

# Clear any stale lock from a previously killed job, then run for real.
snakemake "${args[@]}" --unlock || true
snakemake "${args[@]}"
