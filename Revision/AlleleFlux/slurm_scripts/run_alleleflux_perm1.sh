#!/bin/bash
#SBATCH --job-name=run_alleleflux_perm1
#SBATCH --output=logs/run_alleleflux_perm1_%j.out
#SBATCH --error=logs/run_alleleflux_perm1_%j.err
#SBATCH --time=6-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

# Permuted (null) run 1 — BYO mode.
# Uses the pre-made Fig1 swap sheet (perm_group_swap_set1.tsv, staged as
# <root_dir>/permuted_metadata_fat_control.tsv) and reuses the completed real
# run's profiles/QC/allele-freq cache via input.reuse_from, so only the
# group-dependent tail is recomputed on the permuted labels.

CONFIG=/home/su2806/DietDrivenMicrobiota/Revision/AlleleFlux/alleleflux_config_perm1.yaml

# Snakemake working dir (.snakemake lock + metadata). Deliberately the perm run's
# OWN output root — NOT the real run's workdir — so this run never shares state
# with, or re-touches the reuse sentinels of, the completed real run.
WORKDIR=/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/AlleleFlux_revision/AlleleFlux/permuted/perm1

# Tag every child SLURM job with a WCKey so sacct_stats.py can isolate this
# run via --wckey filter. sbatch reads SBATCH_WCKEY as the default for --wckey
# and the env is inherited by all child sbatch calls plugin-slurm makes.
RUN_NAME=$(grep -E '^run_name:' "$CONFIG" | sed -E 's/.*"([^"]+)".*/\1/')
export SBATCH_WCKEY="alleleflux_${RUN_NAME}"

alleleflux run \
    -w "$WORKDIR" \
    -c "$CONFIG" \
    -p /home/su2806/AlleleFlux-dev/alleleflux/smk_workflow/slurm_profile_native/ --unlock

alleleflux run \
    -w "$WORKDIR" \
    -c "$CONFIG" \
    -p /home/su2806/AlleleFlux-dev/alleleflux/smk_workflow/slurm_profile_native/
