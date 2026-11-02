#!/bin/bash
#SBATCH --job-name=run_alleleflux_midas_ratio10
#SBATCH --output=logs/run_alleleflux_midas_ratio10_%j.out
#SBATCH --error=logs/run_alleleflux_midas_ratio10_%j.err
#SBATCH --time=6-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

# Driver for the MIDAS-filter site_ratio-10 run (profiles reused from the
# ratio-3 run via input.profiles_path; see the config header): Snakemake runs here and submits one SLURM job
# per rule.  Code = frozen clone of AlleleFlux at commit 48e3db6, env =
# alleleflux-run-dm (cloned from alleleflux-dev, pinned to the clone), so Phase 2
# edits in the development tree never reach a running child job.
# -u is switched on only AFTER activation: conda's activate hook reads $PS1,
# which is unset when the job is submitted from a non-interactive shell
# (sbatch forwards the submitter's environment), and under -u that kills the
# driver before it does anything (job 16323991, 2026-10-27).
set -eo pipefail
module load anaconda3/2025.12
conda activate alleleflux-run-dm
set -u

CODE=/home/su2806/AlleleFlux-run-2026-10
CONFIG=/home/su2806/DietDrivenMicrobiota/additional_analyses/AlleleFlux/alleleflux_config_midas_filter_ratio10.yaml
WORKDIR=/scratch/gpfs/AMOELLER/sidd/diet_manip/midas_filter_Oct_2026/AlleleFlux_midas_filter_ratio10

# Refuse to start unless the env is still pinned to the frozen clone.
MAPPED=$(grep -ohP "MAPPING.*?\{'alleleflux': '\K[^']+" "$CONDA_PREFIX"/lib/python*/site-packages/__editable___alleleflux_*_finder.py | head -1)
[[ "$MAPPED" == "$CODE/alleleflux" ]] || { echo "env points at '$MAPPED', not the frozen clone; refusing" >&2; exit 1; }

export SBATCH_WCKEY="alleleflux_AlleleFlux_midas_filter_ratio10"   # tags every child job for sacct
mkdir -p "$WORKDIR"

alleleflux run -w "$WORKDIR" \
    -c "$CONFIG" \
    -p "$CODE/alleleflux/smk_workflow/slurm_profile_native/" --unlock

alleleflux run -w "$WORKDIR" \
    -c "$CONFIG" \
    -p "$CODE/alleleflux/smk_workflow/slurm_profile_native/"
