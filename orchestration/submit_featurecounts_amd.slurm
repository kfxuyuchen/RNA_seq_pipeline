#!/bin/bash
# Thin Slurm wrapper: featureCounts post on existing 03.stringtie_sort BAMs.
# No #SBATCH --mem. AMD_9A14 DefMemPerCPU≈3600 MB → -c 16 ≈ 57G.
#SBATCH -J rnaseq_fc
#SBATCH -p AMD_9A14
#SBATCH -c 16
#SBATCH -t 1-00:00:00
#SBATCH -o logs/rnaseq_fc_%j.out
#SBATCH -e logs/rnaseq_fc_%j.err
set -euo pipefail

# --- edit these ---
OUTDIR="${OUTDIR:?set OUTDIR to rnaseq_quant output root}"
GTF="${GTF:?set GTF annotation}"
PIPE="${PIPE:-/public/home/xuyuchen/local/RNA_seq_pipeline}"
CONDA_ENV="${CONDA_ENV:-/public/home/xuyuchen/miniconda3/envs/RNA-seq}"
THREADS="${THREADS:-16}"
STRAND="${STRAND:-0}"
MAPQ="${MAPQ:-10}"
# ------------------

mkdir -p logs
cd "$OUTDIR"
set +u
source /public/home/xuyuchen/miniconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"
set -u

echo "=== featureCounts start $(date) job=$SLURM_JOB_ID ==="
bash "$PIPE/post_featurecounts.sh" -a "$GTF" -p "$THREADS" -s "$STRAND" -Q "$MAPQ"
bash "$PIPE/run_qc.sh" --title "featureCounts QC" --min-assigned 50 || true
echo "=== featureCounts done $(date) ==="
