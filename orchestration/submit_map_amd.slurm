#!/bin/bash
# Thin Slurm wrapper: HISAT2 map_only (one-write). No #SBATCH --mem.
# AMD_9A14 DefMemPerCPU≈3600 MB → -c 48 ≈ 168G（按样本线程另限 PIPELINE_THREADS）。
#SBATCH -J rnaseq_map
#SBATCH -p AMD_9A14
#SBATCH -c 48
#SBATCH -t 3-00:00:00
#SBATCH -o logs/rnaseq_map_%j.out
#SBATCH -e logs/rnaseq_map_%j.err
set -euo pipefail

# --- edit these ---
OUTDIR="${OUTDIR:?set OUTDIR to rnaseq_quant output root}"
GENOME="${GENOME:?set GENOME fasta}"
FASTQ_DIR="${FASTQ_DIR:?set FASTQ_DIR}"
PIPE="${PIPE:-/public/home/xuyuchen/local/RNA_seq_pipeline}"
CONDA_ENV="${CONDA_ENV:-/public/home/xuyuchen/miniconda3/envs/RNA-seq}"
# Optional: existing HISAT2 index prefix (reuse; do not rebuild/delete)
HISAT2_INDEX="${HISAT2_INDEX:-}"
PIPELINE_THREADS="${PIPELINE_THREADS:-12}"
export SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-1G}"
# ------------------

mkdir -p logs
cd "$OUTDIR"
set +u
source /public/home/xuyuchen/miniconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"
set -u

map_args=(-G "$GENOME" -T "$FASTQ_DIR" -p "$PIPELINE_THREADS")
[[ -n "$HISAT2_INDEX" ]] && map_args+=(-x "$HISAT2_INDEX")
# Wheat / IGV: enable CSI
map_args+=(--index-bam)

echo "=== map start $(date) job=$SLURM_JOB_ID host=$(hostname) ==="
bash "$PIPE/map_only.sh" "${map_args[@]}"
echo "=== map done $(date) ==="
