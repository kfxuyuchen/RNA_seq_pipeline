#!/bin/bash
# LEGACY: StringTie 路线（可选）。主矩阵请用 09.featureCounts。
# 依赖 python2 的 prepDE/getTPM/getFPKM；详见 docs/LEGACY.md。
# StringTie post-quant with resume: skip samples that already have assemble.gtf;
# matrices rebuilt only when needed (or --force).
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERSION="2.2.0-st-resume"
threads=8
tiegff=""
force=0
prepDE_script="${PIPELINE_DIR}/prepDE.py"
getTPM_script="${PIPELINE_DIR}/getTPM.py"
getFPKM_script="${PIPELINE_DIR}/getFPKM.py"
sort_matrix_script="${PIPELINE_DIR}/sort_expr_matrix.py"

usage() {
  cat <<EOF
POST-StringTie v${VERSION}（断点续传）

对 03.stringtie_sort/*_sorted.bam 跑 StringTie -e，再汇总 count/TPM/FPKM。
已有非空 *.assemble.gtf 的样本跳过；矩阵在样本 GTF 集合变化时重算。
矩阵写出后按 gene/transcript ID 字典序排序。

Usage:
  $(basename "$0") -t <gff|gff3> [options]

必需:
  -t, --tiegff PATH      参考注释 GFF/GFF3

可选:
  -p, --threads N        总线程预算；内部并行数≈N/4 (默认: ${threads})
  --force                强制重跑全部样本与矩阵
  --prepDE/--getTPM/--getFPKM PATH   覆盖汇总脚本路径

依赖:
  stringtie, python2 (prepDE/getTPM/getFPKM), python3 (sort_expr_matrix.py)

输出:
  04.stringtie_out/<sample>/*.assemble.gtf
  05.count_out/gene_count.csv  transcript_count.csv
  06.TPM_out/gene_TPM.csv      transcript_TPM.csv
  07.FPKM_out/gene_FPKM.csv    transcript_FPKM.csv

示例:
  bash $(basename "$0") -t genes.gff3 -p 32
  bash $(basename "$0") -t genes.gff3 -p 16 --force
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--tiegff) tiegff="$2"; shift 2 ;;
    -p|--threads) threads="$2"; shift 2 ;;
    --force) force=1; shift ;;
    --prepDE) prepDE_script="$2"; shift 2 ;;
    --getTPM) getTPM_script="$2"; shift 2 ;;
    --getFPKM) getFPKM_script="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "未知选项: $1" >&2; usage ;;
  esac
done

[[ -n "$tiegff" && -f "$tiegff" ]] || usage
command -v stringtie >/dev/null
command -v python2 >/dev/null

st_parallel=$(( threads / 4 ))
[[ "$st_parallel" -lt 1 ]] && st_parallel=1
tiegff=$(realpath "$tiegff")

for dir in 04.stringtie_out 05.count_out 06.TPM_out 07.FPKM_out; do
  mkdir -p "$dir"
done

# clean killed-job temps
find 04.stringtie_out -name '*.tmp' -type f -delete 2>/dev/null || true

n_bam=$(find 03.stringtie_sort -name '*_sorted.bam' | wc -l)
echo "=== POST-StringTie v${VERSION} ==="
echo "cwd=$(pwd)"
echo "sorted BAMs: $n_bam  force=$force"
echo "gff=$tiegff threads=$threads stringtie_parallel=$st_parallel"
[[ "$n_bam" -gt 0 ]] || { echo "ERROR: no sorted BAMs" >&2; exit 1; }

echo -e "\n=== StringTie（可续传）==="
export FORCE_ST="$force"
find 03.stringtie_sort -name '*_sorted.bam' -print0 | xargs -0 -I {} -P "${st_parallel}" sh -c '
  id="$1"
  gff="$2"
  force="$3"
  sample_name=$(basename "$id" "_sorted.bam")
  output_dir="04.stringtie_out/${sample_name}"
  mkdir -p "$output_dir"
  gtf="${output_dir}/${sample_name}.assemble.gtf"
  abund="${output_dir}/${sample_name}_gene_abund.tab"
  if [[ "$force" -eq 0 && -f "$gtf" && -s "$gtf" ]]; then
    echo "定量结果已存在，跳过: $sample_name"
    exit 0
  fi
  echo "定量样本: $sample_name"
  rm -f "${gtf}.tmp" "${abund}.tmp"
  stringtie -e -A "${abund}.tmp" \
    -B -p 8 -G "$gff" -o "${gtf}.tmp" "$id"
  mv -f "${gtf}.tmp" "$gtf"
  [[ -f "${abund}.tmp" ]] && mv -f "${abund}.tmp" "$abund"
' sh {} "$tiegff" "$force"

# fingerprint of current GTF set
sample_lst="05.count_out/.sample_lst"
stamp="05.count_out/.matrix_stamp"
find 04.stringtie_out -name '*.assemble.gtf' -print0 | while IFS= read -r -d '' gtf; do
  sample_name=$(basename "$(dirname "$gtf")")
  echo -e "${sample_name}\t${gtf}"
done | sort > "${sample_lst}.tmp"
mv -f "${sample_lst}.tmp" "$sample_lst"

n_gtf=$(wc -l < "$sample_lst" | tr -d ' ')
new_stamp="$(cksum "$sample_lst" | awk '{print $1"_"$2}')_n${n_gtf}"

need_matrix=1
if [[ "$force" -eq 0 \
   && -f 05.count_out/gene_count.csv && -s 05.count_out/gene_count.csv \
   && -f 06.TPM_out/gene_TPM.csv && -s 06.TPM_out/gene_TPM.csv \
   && -f 07.FPKM_out/gene_FPKM.csv && -s 07.FPKM_out/gene_FPKM.csv \
   && -f "$stamp" && "$(cat "$stamp")" == "$new_stamp" ]]; then
  echo -e "\n=== 矩阵已是最新（stamp=$new_stamp），跳过重算 ==="
  need_matrix=0
fi

if [[ "$need_matrix" -eq 1 ]]; then
  echo -e "\n=== 重建矩阵（n_gtf=$n_gtf）==="
  rm -f 05.count_out/gene_count.csv 05.count_out/transcript_count.csv
  rm -f 06.TPM_out/gene_TPM.csv 06.TPM_out/transcript_TPM.csv
  rm -f 07.FPKM_out/gene_FPKM.csv 07.FPKM_out/transcript_FPKM.csv

  python2 "$prepDE_script" -i "$sample_lst" \
    -g 05.count_out/gene_count.csv.tmp \
    -t 05.count_out/transcript_count.csv.tmp
  mv -f 05.count_out/gene_count.csv.tmp 05.count_out/gene_count.csv
  mv -f 05.count_out/transcript_count.csv.tmp 05.count_out/transcript_count.csv

  python2 "$getTPM_script" -i "$sample_lst" \
    -g 06.TPM_out/gene_TPM.csv.tmp \
    -t 06.TPM_out/transcript_TPM.csv.tmp
  mv -f 06.TPM_out/gene_TPM.csv.tmp 06.TPM_out/gene_TPM.csv
  mv -f 06.TPM_out/transcript_TPM.csv.tmp 06.TPM_out/transcript_TPM.csv

  python2 "$getFPKM_script" -i "$sample_lst" \
    -g 07.FPKM_out/gene_FPKM.csv.tmp \
    -t 07.FPKM_out/transcript_FPKM.csv.tmp
  mv -f 07.FPKM_out/gene_FPKM.csv.tmp 07.FPKM_out/gene_FPKM.csv
  mv -f 07.FPKM_out/transcript_FPKM.csv.tmp 07.FPKM_out/transcript_FPKM.csv

  echo "$new_stamp" > "$stamp"
fi

# 按 gene/transcript ID 字典序排行（样本列已在 sample_lst 中 sort）
echo -e "\n=== 按 ID 排序矩阵行 ==="
command -v python3 >/dev/null || { echo "ERROR: python3 required for sort_expr_matrix.py" >&2; exit 1; }
[[ -f "$sort_matrix_script" ]] || { echo "ERROR: missing $sort_matrix_script" >&2; exit 1; }
python3 "$sort_matrix_script" \
  05.count_out/gene_count.csv \
  05.count_out/transcript_count.csv \
  06.TPM_out/gene_TPM.csv \
  06.TPM_out/transcript_TPM.csv \
  07.FPKM_out/gene_FPKM.csv \
  07.FPKM_out/transcript_FPKM.csv

echo -e "\n=== POST-StringTie 完成 ==="
