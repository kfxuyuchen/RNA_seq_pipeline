#!/bin/bash
# featureCounts post-quant with resume.
# Skips when final matrices exist and BAM set fingerprint unchanged (unless --force).
set -euo pipefail

VERSION="2.1.0-fc-resume"
threads=16
annot=""
mapq=10
strand=0
feature_type="exon"
attr="gene_id"
bam_dir="03.stringtie_sort"
outdir="09.featureCounts"
force=0
counts_to_fpkm_tpm="$(cd "$(dirname "$0")" && pwd)/counts_to_fpkm_tpm.py"

usage() {
  cat <<EOF
POST-featureCounts v${VERSION} (断点续传)
Usage: $(basename "$0") -a <genes.gtf> [options]
  -a|--annot     GTF annotation (required)
  -p|--threads   threads (default: ${threads})
  -Q|--mapq      min MAPQ (default: ${mapq})
  -s|--strand    0/1/2 (default: ${strand})
  -t|--feature   feature type (default: ${feature_type})
  -g|--attr      gene id attr (default: ${attr})
  -b|--bam-dir   BAM dir (default: ${bam_dir})
  -o|--outdir    output dir (default: ${outdir})
  --force        强制重跑（忽略已有结果）
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--annot) annot="$2"; shift 2 ;;
    -p|--threads) threads="$2"; shift 2 ;;
    -Q|--mapq) mapq="$2"; shift 2 ;;
    -s|--strand) strand="$2"; shift 2 ;;
    -t|--feature) feature_type="$2"; shift 2 ;;
    -g|--attr) attr="$2"; shift 2 ;;
    -b|--bam-dir) bam_dir="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage ;;
    *) echo "未知选项: $1" >&2; usage ;;
  esac
done

[[ -n "$annot" && -f "$annot" ]] || usage
command -v featureCounts >/dev/null || {
  echo "ERROR: featureCounts not found. Install subread." >&2
  exit 1
}
command -v python3 >/dev/null
[[ -f "$counts_to_fpkm_tpm" ]] || { echo "ERROR: missing $counts_to_fpkm_tpm" >&2; exit 1; }

annot=$(realpath "$annot")
mkdir -p "$outdir"
rm -f "${outdir}"/*.tmp 2>/dev/null || true

mapfile -t bams < <(find "$bam_dir" -maxdepth 1 -name '*_sorted.bam' | sort)
n_bam=${#bams[@]}
[[ "$n_bam" -gt 0 ]] || { echo "ERROR: no *_sorted.bam in $bam_dir" >&2; exit 1; }

# fingerprint: bam paths + sizes + annot path/mtime
fp_file="${outdir}/.bam_fingerprint"
{
  echo "annot=$(realpath "$annot")"
  stat -c '%n %s %Y' "$annot" 2>/dev/null || stat -f '%N %z %m' "$annot"
  printf '%s\n' "${bams[@]}" | while read -r b; do
    stat -c '%n %s %Y' "$b" 2>/dev/null || stat -f '%N %z %m' "$b"
  done
} > "${fp_file}.new"
new_fp=$(cksum "${fp_file}.new" | awk '{print $1"_"$2}')

raw="${outdir}/gene_counts.txt"
count_tsv="${outdir}/gene_count.tsv"
tpm_tsv="${outdir}/gene_TPM.tsv"
fpkm_tsv="${outdir}/gene_FPKM.tsv"
stamp="${outdir}/.fc_stamp"

echo "=== featureCounts v${VERSION} ==="
echo "cwd=$(pwd)  bams=$n_bam  force=$force"
echo "annot=$annot  out=$outdir"

if [[ "$force" -eq 0 \
   && -f "$raw" && -s "$raw" \
   && -f "$count_tsv" && -s "$count_tsv" \
   && -f "$tpm_tsv" && -s "$tpm_tsv" \
   && -f "$fpkm_tsv" && -s "$fpkm_tsv" \
   && -f "$stamp" && "$(cat "$stamp")" == "$new_fp" ]]; then
  echo "结果已是最新（fingerprint=$new_fp），跳过 featureCounts / FPKM-TPM"
  rm -f "${fp_file}.new"
  echo "=== featureCounts POST 完成（续传跳过）==="
  exit 0
fi

if featureCounts -h 2>&1 | grep -q -- '--countReadPairs'; then
  fc_extra=(-p --countReadPairs)
else
  fc_extra=(-p)
  echo "NOTE: no --countReadPairs; using -p only"
fi

echo -e "\n=== 运行 featureCounts ==="
featureCounts \
  -T "$threads" \
  -a "$annot" \
  -o "${raw}.tmp" \
  -t "$feature_type" \
  -g "$attr" \
  -Q "$mapq" \
  -s "$strand" \
  "${fc_extra[@]}" \
  --primary \
  -C \
  "${bams[@]}"

# featureCounts also writes ${raw}.tmp.summary
mv -f "${raw}.tmp" "$raw"
[[ -f "${raw}.tmp.summary" ]] && mv -f "${raw}.tmp.summary" "${raw}.summary"

echo -e "\n=== FPKM / TPM ==="
python3 "$counts_to_fpkm_tpm" -i "$raw" -o "${outdir}.work" --prefix gene
mv -f "${outdir}.work/gene_count.tsv" "$count_tsv"
mv -f "${outdir}.work/gene_FPKM.tsv" "$fpkm_tsv"
mv -f "${outdir}.work/gene_TPM.tsv" "$tpm_tsv"
mv -f "${outdir}.work/gene_length.tsv" "${outdir}/gene_length.tsv"
rmdir "${outdir}.work" 2>/dev/null || rm -rf "${outdir}.work"

mv -f "${fp_file}.new" "$fp_file"
echo "$new_fp" > "$stamp"

echo -e "\n=== featureCounts POST 完成 ==="
echo "count: $count_tsv"
echo "FPKM:  $fpkm_tsv"
echo "TPM:   $tpm_tsv"
echo "raw:   $raw"
