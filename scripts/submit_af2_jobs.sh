#!/bin/bash
# submit_af2_jobs.sh
# Submits 3 AlphaFold2 jobs in parallel on CSD3.
# Run from: /rds/user/jm2564/hpc-work/chimeric_biosecurity_benchmark
#
# Usage:  bash submit_af2_jobs.sh
#         bash submit_af2_jobs.sh --dry-run   # print sbatch commands without submitting

set -euo pipefail

WORKDIR=/rds/user/jm2564/hpc-work/chimeric_biosecurity_benchmark
FASTA_DIR="${WORKDIR}/chimeras"
OUTPUT_BASE="${WORKDIR}/af2_results"
SCRIPT="${WORKDIR}/scripts/af2_single.sh"
DRY_RUN=false

[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

mkdir -p "${OUTPUT_BASE}"
mkdir -p "${WORKDIR}/logs"

# ── job definitions ──────────────────────────────────────────────────────────
# Format: "JOB_NAME:FASTA_FILE"
declare -a JOBS=(
    "GFP_alone:GFP_alone.fasta"
    "IHFbeta_alone:IHFbeta_alone.fasta"
    "F002_chimera:F002.fasta"
)

# ── submit ───────────────────────────────────────────────────────────────────
declare -a JOB_IDS=()

for entry in "${JOBS[@]}"; do
    JOB_NAME="${entry%%:*}"
    FASTA_FILE="${entry##*:}"
    FASTA_PATH="${FASTA_DIR}/${FASTA_FILE}"
    OUTPUT_DIR="${OUTPUT_BASE}/${JOB_NAME}"

    if [[ ! -f "${FASTA_PATH}" ]]; then
        echo "WARNING: FASTA not found — skipping ${JOB_NAME}: ${FASTA_PATH}"
        continue
    fi

    CMD=(sbatch
        --job-name="af2_${JOB_NAME}"
        --export="ALL,JOB_NAME=${JOB_NAME},FASTA_PATH=${FASTA_PATH},OUTPUT_DIR=${OUTPUT_DIR}"
        "${SCRIPT}"
    )

    if $DRY_RUN; then
        echo "[DRY-RUN] ${CMD[*]}"
    else
        JOB_ID=$(${CMD[@]})
        JOB_ID="${JOB_ID##* }"   # extract numeric ID from "Submitted batch job XXXXXX"
        JOB_IDS+=("${JOB_ID}")
        echo "Submitted ${JOB_NAME} → job ${JOB_ID}  (output: ${OUTPUT_DIR})"
    fi
done

# ── monitoring hint ──────────────────────────────────────────────────────────
if ! $DRY_RUN && [[ ${#JOB_IDS[@]} -gt 0 ]]; then
    IDS_CSV=$(IFS=,; echo "${JOB_IDS[*]}")
    echo ""
    echo "All jobs submitted. Monitor with:"
    echo "  squeue -u jm2564"
    echo "  squeue -j ${IDS_CSV}"
    echo ""
    echo "Watch live:"
    echo "  watch -n 30 'squeue -u jm2564'"
    echo ""
    echo "Tail logs (example for GFP_alone):"
    echo "  tail -f ${WORKDIR}/logs/af2_GFP_alone_*.out"
fi
