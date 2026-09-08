#!/bin/bash
# Map-only: HISAT2 index + align/filter/sort. Supports resume (skip finished samples).
# Incomplete intermediates are cleaned; finals written via tmp → mv.
set -euo pipefail

VERSION="2.1.0-map-resume"
threads=20
debug=0
genome=""
transdb=""
tiegff=""
force=0

usage() {
  echo "RNA-seq MAP-ONLY v${VERSION} (断点续传)"
  echo "Usage: $(basename "$0") -G <genome.fa> -T <fastq_dir> [-t <gff>] [-p threads] [-d] [--force]"
  echo "  已有合法 *_sorted.bam 的样本会跳过；中断产生的 *.tmp / 坏 BAM 会清理后重做"
  echo "  --force  忽略已有 sorted BAM，全部重跑"
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
  # fallback when quickcheck unavailable/fails: header readable
  samtools view -H "$f" >/dev/null 2>&1
}

atomic_mv() {
  local src="$1" dst="$2"
  mv -f "$src" "$dst"
}

trans_files=$(find "$transdb" -maxdepth 1 -name "*clean_1.fastq.gz" -print | sort)
[[ -n "$trans_files" ]] || { echo "ERROR: no *clean_1.fastq.gz in $transdb" >&2; exit 1; }

for dir in 01.ref 02.hisat2-mapping 03.stringtie_sort 08.align_stats; do
  mkdir -p "$dir"
done

# drop leftovers from killed jobs
find 02.hisat2-mapping 03.stringtie_sort -name '*.tmp' -type f -delete 2>/dev/null || true

genome=$(realpath "$genome")
transdb=$(realpath "$transdb")
ref_dir=$(realpath "01.ref")
genome_basename=$(basename "$genome")

echo "=== MAP-ONLY v${VERSION} ==="
echo "genome=$genome"
echo "fastq=$transdb"
echo "threads=$threads sort_mem=$sort_mem force=$force"
echo "cwd=$(pwd)"

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
  sample_name=$(basename "$id" "_1.fastq.gz")
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

  # stale/corrupt final → remove and redo
  if [[ -f "$sorted_bam" ]] && ! bam_ok "$sorted_bam"; then
    echo "  损坏的 sorted BAM，删除重做: $sorted_bam"
    rm -f "$sorted_bam"
    n_fail_clean=$((n_fail_clean + 1))
  fi

  echo "处理样本: $sample_name"
  rm -f "${raw_bam}.tmp" "${unique_bam}.tmp" "${sorted_bam}.tmp"

  # raw BAM
  if bam_ok "$raw_bam"; then
    echo "  原始BAM已存在且有效，跳过比对"
  else
    rm -f "$raw_bam"
    if [[ -f "$legacy_sam" ]]; then
      samtools view -bS -@"${sam_threads}" "$legacy_sam" > "${raw_bam}.tmp"
      atomic_mv "${raw_bam}.tmp" "$raw_bam"
    else
      id2="$(dirname "$id")/$(basename "$id" "_1.fastq.gz")_2.fastq.gz"
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

  # unique BAM
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
