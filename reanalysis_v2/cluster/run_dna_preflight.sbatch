#!/bin/bash
#SBATCH --job-name=bk_dna_preflight
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=dna_preflight_%j.out
#SBATCH --error=dna_preflight_%j.err
set -euo pipefail
: "${PROJECT:?Export PROJECT before submitting this job}"
: "${BIOINFO:?Export BIOINFO before submitting this job}"
module load singularityce || true
export SINGULARITY_BINDPATH=/shares
SCRIPT_DIR="${SCRIPT_DIR:-$PROJECT/02_scripts}"
OUT="$PROJECT/results/revision/DNA_preflight_${SLURM_JOB_ID}"
python3 "$SCRIPT_DIR/dna_preflight.py" --project "$PROJECT" --bioinfo "$BIOINFO" --out "$OUT"
