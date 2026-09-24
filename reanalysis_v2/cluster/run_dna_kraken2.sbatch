#!/bin/bash
# Run one task as a pilot with: sbatch --array=3 run_dna_kraken2.sbatch
# After reviewing the pilot, submit all 22 with: sbatch run_dna_kraken2.sbatch
# Array index 3 is SRR33231997 (moderate-sized input).
#SBATCH --job-name=bk_dna_kraken2
#SBATCH --array=1-22%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=dna_kraken2_%A_%a.out
#SBATCH --error=dna_kraken2_%A_%a.err
set -euo pipefail

PROJECT=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk
DB=${KRAKEN_DB:-/shares/amr.imm.uzh/bioinfo/databases/kraken2_PlusPF}
KRAKEN_MEMORY_MAPPING=${KRAKEN_MEMORY_MAPPING:-0}
SIF=/shares/amr.imm.uzh/bioinfo/singularity/kraken2_2.1.2--pl5321h9f5acd7_3.sif
ROOT="$PROJECT/results/revision/DNA_kraken2_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
TASK_ID=${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array job}

# These are the 22 accession directories in the original shotgun workflow.
mapfile -t SAMPLES < <(find "$PROJECT/results/metaphlan" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' -printf '%f\n' | sort)
if [[ ${#SAMPLES[@]} -ne 22 || $TASK_ID -lt 1 || $TASK_ID -gt ${#SAMPLES[@]} ]]; then
  echo "Expected exactly 22 accession directories and a valid array index; found ${#SAMPLES[@]} directories, index $TASK_ID" >&2
  exit 1
fi
SAMPLE=${SAMPLES[$((TASK_ID - 1))]}
R1="$PROJECT/results/metaphlan/$SAMPLE/${SAMPLE}_UC1.fastq"
R2="$PROJECT/results/metaphlan/$SAMPLE/${SAMPLE}_UC2.fastq"
OUT="$ROOT/$SAMPLE"
mkdir -p "$OUT" "$ROOT/small_reports"

for f in "$R1" "$R2"; do
  if [[ ! -s $f ]]; then echo "Missing or empty input: $f" >&2; exit 1; fi
done
for f in hash.k2d opts.k2d taxo.k2d; do
  if [[ ! -s "$DB/$f" ]]; then echo "Incomplete Kraken2 database: missing $DB/$f" >&2; exit 1; fi
done

if command -v module >/dev/null 2>&1; then
  module load singularityce 2>/dev/null || module load apptainer/1.4.1 2>/dev/null || true
fi
if command -v singularity >/dev/null 2>&1 && [[ -s $SIF ]]; then
  KRAKEN=(singularity exec --bind /shares "$SIF" kraken2)
elif command -v apptainer >/dev/null 2>&1 && [[ -s $SIF ]]; then
  KRAKEN=(apptainer exec --bind /shares "$SIF" kraken2)
elif command -v kraken2 >/dev/null 2>&1; then
  KRAKEN=(kraken2)
else
  echo "No usable Kraken2 executable or Singularity/Apptainer runtime. Ask cluster support about module load singularityce." >&2
  exit 1
fi

case "$KRAKEN_MEMORY_MAPPING" in
  0) MEMORY_OPTION=() ;;
  1) MEMORY_OPTION=(--memory-mapping) ;;
  *) echo "KRAKEN_MEMORY_MAPPING must be 0 or 1" >&2; exit 1 ;;
esac

{
  printf 'sample\t%s\n' "$SAMPLE"
  printf 'array_job_id\t%s\n' "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
  printf 'array_task_id\t%s\n' "$TASK_ID"
  printf 'host\t%s\n' "$(hostname)"
  printf 'start_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'database\t%s\n' "$DB"
  printf 'container\t%s\n' "$SIF"
  printf 'method\tpaired Kraken2; default confidence; minimizer-data report\n'
  printf 'memory_mapping\t%s\n' "$KRAKEN_MEMORY_MAPPING"
  printf 'slurm_memory\t%s\n' "${SLURM_MEM_PER_NODE:-unavailable}"
  printf 'database_files\n'
  stat -c '%n\t%s bytes\t%y' "$DB"/{hash,opts,taxo}.k2d
  printf 'input_files\n'
  stat -c '%n\t%s bytes\t%y' "$R1" "$R2"
  printf 'kraken2_version\n'
  "${KRAKEN[@]}" --version
} > "$OUT/run_metadata.txt"

echo "Classifying $SAMPLE against $DB; results: $OUT"
"${KRAKEN[@]}" \
  --db "$DB" \
  --paired \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  "${MEMORY_OPTION[@]}" \
  --report-minimizer-data \
  --report "$OUT/report.tsv" \
  --output "$OUT/assignments.kraken" \
  "$R1" "$R2" \
  2> "$OUT/kraken_stderr.txt"

gzip -1 "$OUT/assignments.kraken"
printf 'finish_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$OUT/run_metadata.txt"
printf 'Complete: %s\n' "$SAMPLE" > "$OUT/SUCCESS"
tar -czf "$ROOT/small_reports/${SAMPLE}_small_outputs.tar.gz" \
  -C "$OUT" report.tsv run_metadata.txt kraken_stderr.txt SUCCESS
echo "Done: $ROOT/small_reports/${SAMPLE}_small_outputs.tar.gz"
