#!/bin/bash
#SBATCH --job-name=af2_${JOB_NAME}
#SBATCH --output=/rds/user/jm2564/hpc-work/chimeric_biosecurity_benchmark/logs/af2_${JOB_NAME}_%j.out
#SBATCH --error=/rds/user/jm2564/hpc-work/chimeric_biosecurity_benchmark/logs/af2_${JOB_NAME}_%j.err
#SBATCH --partition=icelake
#SBATCH --account=BRYANT-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00
# No GPU directives — icelake is CPU-only

# ── sanity checks ────────────────────────────────────────────────────────────
if [[ -z "${FASTA_PATH}" || -z "${OUTPUT_DIR}" || -z "${JOB_NAME}" ]]; then
    echo "ERROR: FASTA_PATH, OUTPUT_DIR, and JOB_NAME must be set via --export"
    exit 1
fi

# ── environment ──────────────────────────────────────────────────────────────
source /etc/profile.d/modules.sh
module purge
module load alphafold/2.3.2-singularity

WORKDIR=/rds/user/jm2564/hpc-work/chimeric_biosecurity_benchmark
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${WORKDIR}/logs"

echo "========================================"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Job name:   ${JOB_NAME}"
echo "FASTA:      ${FASTA_PATH}"
echo "Output:     ${OUTPUT_DIR}"
echo "Node:       $(hostname)"
echo "CPUs:       ${SLURM_CPUS_PER_TASK}"
echo "Started:    $(date)"
echo "========================================"

# ── AlphaFold2 ───────────────────────────────────────────────────────────────
# --max_template_date=1900-01-01 is essential: skips missing CIF files
# --use_gpu_relax=false + --models_to_relax=none: safe for CPU-only icelake
run_alphafold \
    --fasta_paths="${FASTA_PATH}" \
    --output_dir="${OUTPUT_DIR}" \
    --bfd_database_path=/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb70_database_path=/data/pdb70/pdb70 \
    --template_mmcif_dir=/rds/user/jm2564/hpc-work/chimeric_biosecurity_benchmark/empty_mmcif \
    --model_preset=monomer \
    --db_preset=full_dbs \
    --models_to_relax=none \
    --use_gpu_relax=false \
    --max_template_date=1900-01-01

EXIT_CODE=$?
echo "========================================"
echo "AF2 finished: exit code ${EXIT_CODE}"
echo "Ended:        $(date)"
echo "========================================"
exit ${EXIT_CODE}
