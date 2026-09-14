#!/bin/bash
#SBATCH --job-name=run_alleleflux
#SBATCH --output=logs/run_alleleflux_%j.out
#SBATCH --error=logs/run_alleleflux_%j.err
#SBATCH --time=6-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

CONFIG=/home/su2806/DietDrivenMicrobiota/Revision/AlleleFlux/alleleflux_config.yaml

# Tag every child SLURM job with a WCKey so sacct_stats.py can isolate this
# run via --wckey filter. sbatch reads SBATCH_WCKEY as the default for --wckey
# and the env is inherited by all child sbatch calls plugin-slurm makes.
RUN_NAME=$(grep -E '^run_name:' "$CONFIG" | sed -E 's/.*"([^"]+)".*/\1/')
export SBATCH_WCKEY="alleleflux_${RUN_NAME}"

alleleflux run \
    -w /scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/AlleleFlux_revision \
    -c "$CONFIG" \
    -p /home/su2806/AlleleFlux-dev/alleleflux/smk_workflow/slurm_profile_native/ --unlock

alleleflux run \
    -w /scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/AlleleFlux_revision \
    -c "$CONFIG" \
    -p /home/su2806/AlleleFlux-dev/alleleflux/smk_workflow/slurm_profile_native/
