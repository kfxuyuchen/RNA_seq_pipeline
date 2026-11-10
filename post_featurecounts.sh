#!/bin/bash
# featureCounts post-quant with resume + preflight checks.
# Skips when final matrices exist and BAM set fingerprint unchanged (unless --force).
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERSION="2.2.0-fc-resume"
threads=16
annot=""
mapq=10
strand=0
feature_type="exon"
attr="gene_id"
bam_dir="03.stringtie_sort"
outdir="09.featureCounts"
force=0
skip_preflight=0
min_assigned_pct=50
counts_to_fpkm_tpm="${PIPELINE_DIR}/counts_to_fpkm_tpm.py"
preflight_script="${PIPELINE_DIR}/preflight_check.py"

usage() {
  cat <<EOF
POST-featureCounts v${VERSION}（断点续传 + 开跑前自检）

从 03.stringtie_sort/*_sorted.bam 生成基因 count / FPKM / TPM。
PE 默认 -p --countReadPairs（Subread≥2.0.2）；结果行按 gene_id 排序。

Usage:
  $(basename "$0") -a <genes.gtf> [options]

必需:
  -a, --annot PATH       GTF 注释（需含 -g 属性与 -t 特征类型）

常用:
  -p, --threads N        线程 (默认: ${threads})
  -Q, --mapq N           最低 MAPQ (默认: ${mapq})
  -s, --strand 0|1|2     链特异: 0=非链 1=正向 2=反向 (默认: ${strand})
  -t, --feature TYPE     计数特征 (默认: ${feature_type})
  -g, --attr ATTR        基因 ID 属性 (默认: ${attr})
  -b, --bam-dir DIR      BAM 目录 (默认: ${bam_dir})
  -o, --outdir DIR       输出目录 (默认: ${outdir})
  --force                强制重跑
  --skip-preflight       跳过 BAM/注释染色体与属性检查
  --min-assigned PCT     summary Assigned% 低于此值告警 (默认: ${min_assigned_pct})

开跑前自检 (preflight_check.py):
  - 注释中是否存在 -t 特征与 -g 属性
  - BAM 与注释染色体名是否有交集（避免 chr1 vs 1 全零矩阵）

输出:
  ${outdir}/gene_counts.txt[.summary]
  ${outdir}/gene_count.tsv
  ${outdir}/gene_FPKM.tsv
  ${outdir}/gene_TPM.tsv
  ${outdir}/gene_length.tsv

示例:
  bash $(basename "$0") -a genes.gtf -p 32
  bash $(basename "$0") -a genes.gtf -s 2 -Q 10 -g gene_id -t exon -p 16
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
    --skip-preflight) skip_preflight=1; shift ;;
    --min-assigned) min_assigned_pct="$2"; shift 2 ;;
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
command -v samtools >/dev/null
[[ -f "$counts_to_fpkm_tpm" ]] || { echo "ERROR: missing $counts_to_fpkm_tpm" >&2; exit 1; }
[[ -f "$preflight_script" ]] || { echo "ERROR: missing $preflight_script" >&2; exit 1; }

annot=$(realpath "$annot")
mkdir -p "$outdir"
# only clear our own tmp markers; do not wipe featureCounts temp-core files mid-run
rm -f "${outdir}"/gene_counts.txt.tmp "${outdir}"/gene_counts.txt.tmp.summary 2>/dev/null || true

mapfile -t bams < <(find "$bam_dir" -maxdepth 1 -name '*_sorted.bam' | sort)
n_bam=${#bams[@]}
[[ "$n_bam" -gt 0 ]] || { echo "ERROR: no *_sorted.bam in $bam_dir" >&2; exit 1; }

fp_file="${outdir}/.bam_fingerprint"
{
  echo "annot=$(realpath "$annot")"
  echo "strand=$strand mapq=$mapq feature=$feature_type attr=$attr"
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
echo "strand=$strand mapq=$mapq feature=$feature_type attr=$attr"

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

if [[ "$skip_preflight" -eq 0 ]]; then
  echo -e "\n=== preflight ==="
  python3 "$preflight_script" \
    --annot "$annot" \
    --bam "${bams[0]}" \
    -g "$attr" \
    -t "$feature_type"
fi

# Detect library layout: SE BAMs must not use -p
has_pe=0; has_se=0
for b in "${bams[@]}"; do
  paired=$(samtools flagstat "$b" 2>/dev/null | awk '/paired in sequencing/{print $1; exit}')
  if [[ "${paired:-0}" -gt 0 ]]; then has_pe=1; else has_se=1; fi
done
if [[ "$has_pe" -eq 1 && "$has_se" -eq 1 ]]; then
  echo "ERROR: 混合了 PE 和 SE BAM；请分别建目录定量（per-sample 目录应为单一 layout）" >&2
  exit 1
fi

fc_help=$(featureCounts -h 2>&1 || true)
if [[ "$has_se" -eq 1 ]]; then
  fc_extra=()
  echo "NOTE: single-end BAM detected; 不使用 -p/--countReadPairs"
elif grep -q -- '--countReadPairs' <<<"$fc_help"; then
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

mv -f "${raw}.tmp" "$raw"
[[ -f "${raw}.tmp.summary" ]] && mv -f "${raw}.tmp.summary" "${raw}.summary"

# Assigned% gate (warn only)
if [[ -f "${raw}.summary" ]]; then
  python3 - "$raw.summary" "$min_assigned_pct" <<'PY'
import sys
from pathlib import Path
path, thr = Path(sys.argv[1]), float(sys.argv[2])
lines = path.read_text().splitlines()
hdr = lines[0].split("\t")
rows = {l.split("\t")[0]: l.split("\t")[1:] for l in lines[1:]}
assigned = rows.get("Assigned")
if not assigned:
    sys.exit(0)
low = []
for i, a in enumerate(assigned):
    vals = []
    for k, v in rows.items():
        try:
            vals.append(float(v[i]))
        except Exception:
            pass
    den = sum(vals) or 1.0
    pct = 100.0 * float(a) / den
    if pct < thr:
        name = Path(hdr[i+1]).name.replace("_sorted.bam","").replace(".bam","")
        low.append((name, pct))
if low:
    print(f"WARN: {len(low)} sample(s) Assigned% < {thr}")
    for n, p in low:
        print(f"  {n}: {p:.1f}%")
PY
fi

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
