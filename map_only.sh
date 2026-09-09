#!/bin/bash
# Map-only: HISAT2 index + align/filter/sort. Supports resume (skip finished samples).
# Incomplete intermediates are cleaned; finals written via tmp → mv.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=pipeline_common.sh
source "${PIPELINE_DIR}/pipeline_common.sh"

VERSION="2.2.0-map-resume"
threads=20
debug=0
genome=""
transdb=""
tiegff=""
force=0

usage() {
  cat <<EOF
RNA-seq MAP-ONLY v${VERSION}（断点续传）

HISAT2 (--dta --mm) → unique (NH:i:1 + MAPQ≥10) → coordinate-sorted BAM。
已有合法 *_sorted.bam 的样本自动跳过；中断产生的 *.tmp / 坏 BAM 会清理后重做。

Usage:
  $(basename "$0") -G <genome.fa> -T <fastq_dir> [options]

必需:
  -G, --genome PATH      参考基因组 FASTA
  -T, --transdb DIR      PE FASTQ 目录

可选:
  -t, --tiegff PATH      注释 GFF（仅记录/兼容旧接口；比对本身不强制）
  -p, --threads N        线程数 (默认: ${threads})
  -d, --debug            保留 raw/unique 中间 BAM
  --force                忽略已有 sorted BAM，全部重跑

FASTQ 命名（提高跨项目兼容性）:
  --r1-suffix STR        R1 后缀 (默认: ${R1_SUFFIX})
  --r2-suffix STR        R2 后缀 (默认: ${R2_SUFFIX})
  --r1-glob  GLOB        在 -T 目录下查找 R1 的 glob (默认: ${R1_GLOB})

  默认兼容本实验室 *clean_1.fastq.gz / *clean_2.fastq.gz。
  其它常见写法示例:
    --r1-suffix _1.fq.gz --r2-suffix _2.fq.gz --r1-glob '*_1.fq.gz'
    --r1-suffix _R1.fastq.gz --r2-suffix _R2.fastq.gz --r1-glob '*_R1.fastq.gz'

环境变量:
  SAMTOOLS_SORT_MEM      samtools sort -m (默认: 1G)

输出:
  01.ref/                  HISAT2 索引
  03.stringtie_sort/       *_sorted.bam（下游定量共用）
  08.align_stats/          raw/unique flagstat
  02.hisat2-mapping/       debug 时保留 raw/unique BAM

示例:
  bash $(basename "$0") -G genome.fa -T clean_fastq -p 32
  bash $(basename "$0") -G genome.fa -T fq \\
    --r1-suffix _1.fq.gz --r2-suffix _2.fq.gz --r1-glob '*_1.fq.gz' -p 16
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -G|--genome) genome="$2"; shift 2 ;;
    -T|--transdb) transdb="$2"; shift 2 ;;
    -t|--tiegff) tiegff="$2"; shift 2 ;;
    -p|--threads) threads="$2"; shift 2 ;;
    -d|--debug) debug=1; shift ;;
    --force) force=1; shift ;;
    --r1-suffix) R1_SUFFIX="$2"; shift 2 ;;
    --r2-suffix) R2_SUFFIX="$2"; shift 2 ;;
    --r1-glob) R1_GLOB="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "未知选项: $1" >&2; usage ;;
  esac
done

[[ -n "$genome" && -n "$transdb" ]] || usage
command -v hisat2 >/dev/null
command -v hisat2-build >/dev/null
command -v samtools >/dev/null
[[ -f "$genome" && -d "$transdb" ]]
[[ -z "$tiegff" || -f "$tiegff" ]] || { echo "ERROR: -t file not found: $tiegff" >&2; exit 1; }

sam_threads=$(( threads > 4 ? 4 : threads ))
sort_threads=$(( threads > 4 ? 4 : threads ))
sort_mem="${SAMTOOLS_SORT_MEM:-1G}"

bam_ok() {
  local f="$1"
  [[ -f "$f" && -s "$f" ]] || return 1
  if samtools quickcheck -q "$f" 2>/dev/null; then
    return 0
  fi
  samtools view -H "$f" >/dev/null 2>&1
}

atomic_mv() {
  local src="$1" dst="$2"
  mv -f "$src" "$dst"
}

trans_files=$(find_r1_files "$transdb")
[[ -n "$trans_files" ]] || {
  echo "ERROR: no files matching R1_GLOB='$R1_GLOB' in $transdb" >&2
  echo "  hint: set --r1-glob / --r1-suffix to match your FASTQ names" >&2
  exit 1
}

# PE mate preflight
while IFS= read -r id; do
  [[ -n "$id" ]] || continue
  id2=$(fq_mate_path "$id") || exit 1
  [[ -f "$id2" ]] || { echo "ERROR: missing mate for $(basename "$id"): $id2" >&2; exit 1; }
done <<< "$trans_files"

for dir in 01.ref 02.hisat2-mapping 03.stringtie_sort 08.align_stats; do
  mkdir -p "$dir"
done

find 02.hisat2-mapping 03.stringtie_sort -name '*.tmp' -type f -delete 2>/dev/null || true

genome=$(realpath "$genome")
transdb=$(realpath "$transdb")
ref_dir=$(realpath "01.ref")
genome_basename=$(basename "$genome")

echo "=== MAP-ONLY v${VERSION} ==="
echo "genome=$genome"
echo "fastq=$transdb"
echo "R1_GLOB=$R1_GLOB  R1_SUFFIX=$R1_SUFFIX  R2_SUFFIX=$R2_SUFFIX"
echo "threads=$threads sort_mem=$sort_mem force=$force"
echo "cwd=$(pwd)"
n_r1=$(grep -c . <<<"$trans_files" || true)
echo "R1 samples: $n_r1"

echo -e "\n=== 步骤1: 参考索引 ==="
if [[ ! -f "${ref_dir}/${genome_basename}.1.ht2l" && ! -f "${ref_dir}/${genome_basename}.1.ht2" ]]; then
  ln -sf "$genome" "${ref_dir}/"
  hisat2-build -p "$threads" "${ref_dir}/${genome_basename}" "${ref_dir}/${genome_basename}"
else
  echo "索引已存在，跳过"
fi

echo -e "\n=== 步骤2: 比对 / 筛选 / 排序（可续传）==="
n_skip=0
n_done=0
n_fail_clean=0
while IFS= read -r id; do
  [[ -n "$id" ]] || continue
  sample_name=$(fq_sample_name "$id") || exit 1
  legacy_sam="02.hisat2-mapping/${sample_name}.sam"
  raw_bam="02.hisat2-mapping/${sample_name}_raw.bam"
  unique_bam="02.hisat2-mapping/${sample_name}_unique.bam"
  sorted_bam="03.stringtie_sort/${sample_name}_sorted.bam"
  stats_raw="08.align_stats/${sample_name}_raw.flagstat.txt"
  stats_unique="08.align_stats/${sample_name}_unique.flagstat.txt"

  if [[ "$force" -eq 0 ]] && bam_ok "$sorted_bam"; then
    echo "已完成，跳过: $sample_name"
    n_skip=$((n_skip + 1))
    continue
  fi

  if [[ -f "$sorted_bam" ]] && ! bam_ok "$sorted_bam"; then
    echo "  损坏的 sorted BAM，删除重做: $sorted_bam"
    rm -f "$sorted_bam"
    n_fail_clean=$((n_fail_clean + 1))
  fi

  echo "处理样本: $sample_name"
  rm -f "${raw_bam}.tmp" "${unique_bam}.tmp" "${sorted_bam}.tmp"

  if bam_ok "$raw_bam"; then
    echo "  原始BAM已存在且有效，跳过比对"
  else
    rm -f "$raw_bam"
    if [[ -f "$legacy_sam" ]]; then
      samtools view -bS -@"${sam_threads}" "$legacy_sam" > "${raw_bam}.tmp"
      atomic_mv "${raw_bam}.tmp" "$raw_bam"
    else
      id2=$(fq_mate_path "$id") || exit 1
      [[ -f "$id2" ]] || { echo "ERROR: missing mate $id2" >&2; exit 1; }
      echo "  hisat2 → samtools view (pipe BAM) [-p ${threads}, --dta --mm]"
      hisat2 -x "${ref_dir}/${genome_basename}" -p "$threads" \
        --dta --mm \
        --fr -1 "$id" -2 "$id2" \
        | samtools view -@"${sam_threads}" -bS - > "${raw_bam}.tmp"
      atomic_mv "${raw_bam}.tmp" "$raw_bam"
    fi
  fi

  [[ -f "$stats_raw" && -s "$stats_raw" ]] || samtools flagstat -@"${sam_threads}" "$raw_bam" > "$stats_raw"

  if bam_ok "$unique_bam"; then
    echo "  unique BAM 已存在且有效，跳过筛选"
  else
    rm -f "$unique_bam"
    echo "  筛选unique mapping → BAM"
    samtools view -h -q 10 -@"${sam_threads}" "$raw_bam" | \
      grep -E "^@|NH:i:1$|NH:i:1[^0-9]" | \
      samtools view -bS -@"${sam_threads}" - > "${unique_bam}.tmp"
    atomic_mv "${unique_bam}.tmp" "$unique_bam"
  fi

  [[ -f "$stats_unique" && -s "$stats_unique" ]] || samtools flagstat -@"${sam_threads}" "$unique_bam" > "$stats_unique"

  if [[ "$debug" -eq 0 ]]; then
    rm -f "$raw_bam" "$legacy_sam"
  fi

  echo "  排序BAM: $sorted_bam (sort -m ${sort_mem} -@${sort_threads})"
  samtools sort -m "${sort_mem}" -@"${sort_threads}" -o "${sorted_bam}.tmp" "$unique_bam"
  atomic_mv "${sorted_bam}.tmp" "$sorted_bam"
  bam_ok "$sorted_bam" || { echo "ERROR: sorted BAM invalid after write: $sorted_bam" >&2; exit 1; }

  if [[ "$debug" -eq 0 ]]; then
    rm -f "$unique_bam"
  fi

  echo "  完成: $sample_name"
  n_done=$((n_done + 1))
done <<< "$trans_files"

echo -e "\n=== MAP-ONLY 完成 ==="
echo "本任务新完成: $n_done  已存在跳过: $n_skip  清理坏文件: $n_fail_clean"
echo "（未跑定量；续跑同一命令即可跳过已完成样本）"
