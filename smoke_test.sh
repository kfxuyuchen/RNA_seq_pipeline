#!/bin/bash
# Official smoke (≤5 min): toy genome + PE reads → map_only onewrite + rate table + QC.
# Exit 0 + SMOKE_PASS means the core map/QC contract is green.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
E="${RNASEQ_BIN:-/public/home/xuyuchen/miniconda3/envs/RNA-seq/bin}"
export PATH="$E:$PATH"

ROOT="${TMPDIR:-/tmp}/rnaseq_pipeline_smoke_$$"
rm -rf "$ROOT"
mkdir -p "$ROOT/fq" "$ROOT/out"

python3 - <<PY
from pathlib import Path
import gzip
root = Path("$ROOT")
seq = ("ACGT" * 250)
(root / "toy.fa").write_text(">chr1\n" + seq + "\n")
r1, r2 = [], []
for i in range(40):
    s = i * 8
    read = seq[s : s + 40]
    if len(read) < 40:
        break
    r1 += [f"@r{i}/1", read, "+", "I" * 40]
    r2 += [f"@r{i}/2", read, "+", "I" * 40]
with gzip.open(root / "fq" / "toy_clean_1.fastq.gz", "wt") as f:
    f.write("\n".join(r1) + "\n")
with gzip.open(root / "fq" / "toy_clean_2.fastq.gz", "wt") as f:
    f.write("\n".join(r2) + "\n")
print(f"smoke data: {root}")
PY

cd "$ROOT/out"
bash "$PIPELINE_DIR/map_only.sh" -G "$ROOT/toy.fa" -T "$ROOT/fq" -p 2 --index-bam

# Continuous dirs present
for d in 01.ref 02.hisat2-mapping 03.stringtie_sort 08.align_stats; do
  [[ -d "$d" ]] || { echo "FAIL missing $d"; exit 1; }
done
# One-write: no intermediate BAMs
n_inter=$(find 02.hisat2-mapping -name '*_raw.bam' -o -name '*_unique.bam' | wc -l)
[[ "$n_inter" -eq 0 ]] || { echo "FAIL unexpected intermediate BAM"; exit 1; }
[[ -s 02.hisat2-mapping/toy_clean.hisat2.log ]]
bam=03.stringtie_sort/toy_clean_sorted.bam
samtools quickcheck "$bam"
samtools view -H "$bam" | grep -q 'SO:coordinate'
[[ -s 08.align_stats/hisat2_alignment_rates.tsv ]]
grep -qE $'\t[0-9.]+\t' 08.align_stats/hisat2_alignment_rates.tsv
[[ -s "${bam}.csi" ]] || { echo "FAIL missing CSI"; exit 1; }

bash "$PIPELINE_DIR/run_qc.sh" --no-multiqc --min-hisat2 50
[[ -s 10.qc/hisat2_alignment_rates.tsv ]]
[[ -s 10.qc/qc_overview.tsv ]]
grep -q hisat2_overall_pct 10.qc/qc_overview.tsv

echo "SMOKE_PASS root=$ROOT"
rm -rf "$ROOT"
echo "SMOKE_CLEANED"
