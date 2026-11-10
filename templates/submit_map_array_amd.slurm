#!/bin/bash
#SBATCH -J map_array
#SBATCH -p AMD_9A14,AMD_9654
#SBATCH -c 12
#SBATCH -t 2-00:00:00

# Parallel mapping: one array task = one sample.
# Required env:
#   OUTDIR, FASTQ_DIR, GENOME, INDEX, SAMPLE_LIST
# Optional env:
#   THREADS (default 12), SAMTOOLS_THREADS (default 4), SAMTOOLS_SORT_MEM (default 1G)
set -euo pipefail
OUTDIR="${OUTDIR:?set OUTDIR}"
FASTQ_DIR="${FASTQ_DIR:?set FASTQ_DIR}"
GENOME="${GENOME:?set GENOME}"
INDEX="${INDEX:?set INDEX}"
SAMPLE_LIST="${SAMPLE_LIST:?set SAMPLE_LIST}"
THREADS="${THREADS:-12}"
export SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-4}"
export SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-1G}"

source /public/home/xuyuchen/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

PIPE=/public/home/xuyuchen/local/RNA_seq_pipeline
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")
[[ -n "$SAMPLE" ]] || { echo "ERROR: empty sample at index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

mkdir -p "$OUTDIR/logs"
cd "$OUTDIR"
echo "task=${SLURM_ARRAY_TASK_ID} sample=${SAMPLE} node=${SLURMD_NODENAME:-?} threads=${THREADS}"
bash "$PIPE/map_only.sh" -G "$GENOME" -T "$FASTQ_DIR" -p "$THREADS" -x "$INDEX" --only "$SAMPLE"
