#!/bin/bash
#SBATCH --job-name=baseline_presence
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/home/su2806/DietDrivenMicrobiota/Revision/baseline_presence/logs/%x_%j.out
#SBATCH --error=/home/su2806/DietDrivenMicrobiota/Revision/baseline_presence/logs/%x_%j.err

COMPARISON=pre_end-fat_control

# ---- inputs (the mapq20 run)
RUN_DIR=/scratch/gpfs/AMOELLER/diet_manip/AlleleFlux_mapq20/longitudinal
PROFILES_DIR=/scratch/gpfs/AMOELLER/diet_manip/AlleleFlux_mapq20/longitudinal/profiles
FASTA=/scratch/gpfs/AMOELLER/diet_manip/redo_representative_megamag.fa       # contig reference (NOT the prodigal .fna)
MAG_MAPPING=/scratch/gpfs/AMOELLER/diet_manip/megamag_mapping.tsv
METADATA=/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/metadata_md_bam_fixed.tsv

# ---- outputs
OUTPUT_DIR=/scratch/gpfs/AMOELLER/sidd/diet_manip/revision_July_2026/baseline_presence
mkdir -p "$OUTPUT_DIR"

# ---- arguments (comments live here; a comment after a trailing backslash breaks the command)
ARGS=(
    --run_dir      "$RUN_DIR"
    --comparison   "$COMPARISON"
    --summary      two_sample_paired          # divergence: the paper's canonical test family
    --test_type    two_sample_paired_tTest    # exactly as the p_value_summary spells it
    --threshold_column q_value                # BH-corrected
    --threshold    0.05
    --profiles_dir "$PROFILES_DIR"
    --metadata     "$METADATA"
    --fasta        "$FASTA"
    --mag_mapping  "$MAG_MAPPING"
    --output_dir   "$OUTPUT_DIR"
    --min_cov      5                          # reads at the position before any verdict (1 = off)
    --min_freq     0.05                       # presence floor on top of the null-model bar
    --cpus         "$SLURM_CPUS_PER_TASK"
    # --turnover_dir <strain_turnover outputs>   # optional: adds the strain_background column
)

echo "[$(date)] $COMPARISON  metadata=$METADATA  cpus=$SLURM_CPUS_PER_TASK"
alleleflux-baseline-presence "${ARGS[@]}"
echo "[$(date)] done: $(ls "$OUTPUT_DIR"/${COMPARISON}_two_sample_paired_tTest_baseline_presence*)"
