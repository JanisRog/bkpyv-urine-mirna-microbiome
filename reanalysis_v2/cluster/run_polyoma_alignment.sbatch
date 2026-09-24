#!/bin/bash
# Pilot the two key samples: sbatch --array=15,20 run_polyoma_alignment.sbatch
# After reviewing that pilot: sbatch run_polyoma_alignment.sbatch
#SBATCH --job-name=bk_jc_align
#SBATCH --array=1-22%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=polyoma_align_%A_%a.out
#SBATCH --error=polyoma_align_%A_%a.err
set -euo pipefail

PROJECT=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk
SCRIPT_DIR="$PROJECT/02_scripts"
REF_SOURCE="$SCRIPT_DIR/BK_JC_RefSeq.fasta"
BWA_SIF=/shares/amr.imm.uzh/bioinfo/singularity/quay.io-biocontainers-bwa-0.7.17--h5bf99c6_8.img
SAM_SIF=/shares/amr.imm.uzh/bioinfo/singularity/samtools_1.15.1--h1170115_0.sif
TASK_ID=${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array job}
ROOT="$PROJECT/results/revision/DNA_polyoma_alignment_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"

mapfile -t SAMPLES < <(find "$PROJECT/results/metaphlan" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' -printf '%f\n' | sort)
if [[ ${#SAMPLES[@]} -ne 22 || $TASK_ID -lt 1 || $TASK_ID -gt ${#SAMPLES[@]} ]]; then
  echo "Expected exactly 22 accession directories and valid array index; found ${#SAMPLES[@]}, index $TASK_ID" >&2
  exit 1
fi
SAMPLE=${SAMPLES[$((TASK_ID - 1))]}
R1="$PROJECT/results/metaphlan/$SAMPLE/${SAMPLE}_UC1.fastq"
R2="$PROJECT/results/metaphlan/$SAMPLE/${SAMPLE}_UC2.fastq"
OUT="$ROOT/$SAMPLE"
mkdir -p "$OUT" "$ROOT/small_reports"
for f in "$REF_SOURCE" "$R1" "$R2" "$BWA_SIF" "$SAM_SIF" "$SCRIPT_DIR/summarize_polyoma_alignment.py"; do
  if [[ ! -s $f ]]; then echo "Missing or empty input: $f" >&2; exit 1; fi
done
if command -v module >/dev/null 2>&1; then
  module load singularityce 2>/dev/null || module load apptainer/1.4.1 2>/dev/null || true
fi
if command -v singularity >/dev/null 2>&1; then
  ENGINE=(singularity exec --bind /shares)
elif command -v apptainer >/dev/null 2>&1; then
  ENGINE=(apptainer exec --bind /shares)
else
  echo "No Singularity/Apptainer runtime available" >&2
  exit 1
fi
BWA=("${ENGINE[@]}" "$BWA_SIF" bwa)
SAM=("${ENGINE[@]}" "$SAM_SIF" samtools)

REF="$OUT/BK_JC_RefSeq.fasta"
cp "$REF_SOURCE" "$REF"
if [[ $(sha256sum "$REF" | cut -d' ' -f1) != 26dc1b8c7117a6e52b72cc17dba547887f7e660b4fb2aeae7e38076fdfafdb23 ]]; then
  echo "Reference FASTA checksum mismatch" >&2
  exit 1
fi
{
  printf 'sample\t%s\n' "$SAMPLE"
  printf 'array_job_id\t%s\n' "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
  printf 'array_task_id\t%s\n' "$TASK_ID"
  printf 'start_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'reference_sha256\t%s\n' "$(sha256sum "$REF" | cut -d' ' -f1)"
  printf 'bwa_version\n'
  "${BWA[@]}" 2>&1 || true
  printf 'samtools_version\n'
  "${SAM[@]}" --version
  printf 'input_files\n'
  stat -c '%n %s bytes %y' "$R1" "$R2"
} > "$OUT/run_metadata.txt"

"${BWA[@]}" index "$REF" > "$OUT/bwa_index_stdout.txt" 2> "$OUT/bwa_index_stderr.txt"
"${BWA[@]}" mem -t "${SLURM_CPUS_PER_TASK:-8}" "$REF" "$R1" "$R2" 2> "$OUT/bwa_mem_stderr.txt" |
  "${SAM[@]}" view -@ 2 -u -F 0x904 - |
  "${SAM[@]}" sort -@ 2 -m 1G -o "$OUT/primary_mapped.bam" -
"${SAM[@]}" index -@ 2 "$OUT/primary_mapped.bam"
"${SAM[@]}" idxstats "$OUT/primary_mapped.bam" > "$OUT/idxstats.tsv"

{
  printf 'reference\tmapq30_first_mates\n'
  for ref in NC_001538.1 NC_001699.1; do
    count=$("${SAM[@]}" view -c -f 0x40 -q 30 "$OUT/primary_mapped.bam" "$ref")
    printf '%s\t%s\n' "$ref" "$count"
  done
} > "$OUT/mapq30_first_mates.tsv"

# In samtools depth, -q is base quality and -Q is mapping quality.
"${SAM[@]}" depth -aa -q 20 -Q 30 "$OUT/primary_mapped.bam" > "$OUT/depth_mapq30_baseq20.tsv"
python3 "$SCRIPT_DIR/summarize_polyoma_alignment.py" \
  --sample "$SAMPLE" --reference "$REF" \
  --depth "$OUT/depth_mapq30_baseq20.tsv" \
  --idxstats "$OUT/idxstats.tsv" \
  --counts "$OUT/mapq30_first_mates.tsv" \
  --out "$OUT/summary.tsv"
printf 'finish_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$OUT/run_metadata.txt"
printf 'Complete: %s\n' "$SAMPLE" > "$OUT/SUCCESS"
tar -czf "$ROOT/small_reports/${SAMPLE}_polyoma_summary.tar.gz" \
  -C "$OUT" summary.tsv run_metadata.txt idxstats.tsv mapq30_first_mates.tsv \
  depth_mapq30_baseq20.tsv bwa_index_stderr.txt bwa_mem_stderr.txt SUCCESS
echo "Done: $ROOT/small_reports/${SAMPLE}_polyoma_summary.tar.gz"
