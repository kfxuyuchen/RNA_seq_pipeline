#!/bin/bash
# Submit one SLURM array task per sample (parallel mapping).
# Usage: submit_map_array.sh <sample_list> [max_parallel]
set -euo pipefail
LIST=${1:?sample list}
MAXJ=${2:-8}
N=$(sed '/^[[:space:]]*$/d' "$LIST" | wc -l)
[[ "$N" -gt 0 ]] || { echo "empty list" >&2; exit 1; }
sbatch --array="0-$((N-1))%${MAXJ}" --export=ALL orchestration/submit_map_array_amd.slurm
echo "submitted N=${N} max_parallel=${MAXJ}"
