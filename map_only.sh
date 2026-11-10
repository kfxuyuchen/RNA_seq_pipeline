#!/bin/bash
# Map-only: HISAT2 index + align/filter/sort. Supports resume (skip finished samples).
# Default: single disk write → 03.stringtie_sort/*_sorted.bam (hisat2 | unique | sort).
# --debug: also keep raw/unique under 02.hisat2-mapping/ (legacy multi-write path).
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=pipeline_common.sh
source "${PIPELINE_DIR}/pipeline_common.sh"

VERSION="2.4.2"
threads=20
debug=0
genome=""
transdb=""
tiegff=""
force=0
index_bam=0
# Optional external HISAT2 index prefix (no .1.ht2 / .1.ht2l suffix). Prefer reuse over rebuild.
index_prefix=""

usage() {
  cat <<EOF
RNA-seq MAP-ONLY v${VERSION}（断点续传 / 默认一次落盘）

默认: HISAT2 (--dta --mm) | unique(NH:i:1 + MAPQ≥10) | sort
  → 只写 03.stringtie_sort/*_sorted.bam
  → 02.hisat2-mapping/ 只写 hisat2 日志（目录编号连续保留）
  → 08.align_stats/*_unique.flagstat.txt 来自最终 sorted BAM

--debug: 恢复旧式 raw→unique→sort 分步落盘，并保留 02 下中间 BAM。

已有合法 *_sorted.bam 的样本自动跳过；*.tmp / 坏 BAM 会清理后重做。

Usage:
  $(basename "$0") -G <genome.fa> -T <fastq_dir> [options]

必需:
  -G, --genome PATH      参考基因组 FASTA
  -T, --transdb DIR      PE FASTQ 目录

可选:
  -t, --tiegff PATH      注释 GFF（仅记录/兼容旧接口；比对本身不强制）
  -x, --index PREFIX     已有 HISAT2 索引前缀（优先复用；勿擅自删除旧索引）
  -p, --threads N        线程数 (默认: ${threads})
  -d, --debug            多步落盘并保留 raw/unique 中间 BAM
  --index-bam            对 sorted BAM 建 CSI 索引（samtools index -c；小麦超长染色体必须用 CSI，勿用默认 BAI）
  --force                忽略已有 sorted BAM，全部重跑

FASTQ 命名（提高跨项目兼容性）:
  --r1-suffix STR        R1 后缀 (默认: ${R1_SUFFIX})
  --r2-suffix STR        R2 后缀 (默认: ${R2_SUFFIX})
  --r1-glob  GLOB        在 -T 目录下查找 R1 的 glob (默认: ${R1_GLOB})

  默认兼容本实验室 *clean_1.fastq.gz / *clean_2.fastq.gz。

环境变量:
  SAMTOOLS_SORT_MEM      samtools sort -m (默认: 1G)
  MAPQ_MIN              唯一比对最低 MAPQ（默认: 10）
  UNIQUE_ONLY           1=只保留 NH:i:1（默认），0=保留多重比对
  HISAT2_EXTRA         追加给 hisat2 的额外参数（如 "--max-intronlen 5000"）
  HISAT2_INDEX           同 -x（命令行优先）

输出（编号连续）:
  01.ref/                HISAT2 索引（或链到已有索引）
  02.hisat2-mapping/     *.hisat2.log；--debug 时另有 *_raw/_unique.bam
  03.stringtie_sort/     *_sorted.bam（下游定量共用，默认唯一 BAM）；--index-bam 时另有 *.csi
  08.align_stats/        *_unique.flagstat.txt；*_hisat2_rate.tsv；hisat2_alignment_rates.tsv
                         （debug 时另有 *_raw.flagstat.txt）

示例:
  bash $(basename "$0") -G genome.fa -T clean_fastq -p 32
  bash $(basename "$0") -G genome.fa -T fq -x /path/to/CS-IAAS_v1.1.fa -p 16 --index-bam
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -G|--genome) genome="$2"; shift 2 ;;
    -T|--transdb) transdb="$2"; shift 2 ;;
    -t|--tiegff) tiegff="$2"; shift 2 ;;
    -x|--index) index_prefix="$2"; shift 2 ;;
    -p|--threads) threads="$2"; shift 2 ;;
    -d|--debug) debug=1; shift ;;
    --index-bam) index_bam=1; shift ;;
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

if [[ -z "$index_prefix" && -n "${HISAT2_INDEX:-}" ]]; then
  index_prefix="$HISAT2_INDEX"
fi

sam_threads=$(( threads > 4 ? 4 : threads ))
sort_threads=$(( threads > 4 ? 4 : threads ))
sort_mem="${SAMTOOLS_SORT_MEM:-1G}"
mapq_min="${MAPQ_MIN:-10}"
unique_only="${UNIQUE_ONLY:-1}"
hisat2_extra="${HISAT2_EXTRA:-}"
if [[ "$unique_only" -eq 1 ]]; then
  unique_awk='BEGIN{FS=OFS="\t"} /^@/ || $0 ~ /NH:i:1($|\t)/'
else
  unique_awk='BEGIN{FS=OFS="\t"} /^@/ || 1'
fi

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

index_ok() {
  local p="$1"
  [[ -f "${p}.1.ht2l" || -f "${p}.1.ht2" ]]
}

# Write per-sample HISAT2 overall rate under 08.align_stats/ (from 02/*.hisat2.log).
write_hisat2_rate() {
  local sample="$1" logf="$2"
  local rate_tsv="08.align_stats/${sample}_hisat2_rate.tsv"
  local rate=""
  if [[ -s "$logf" ]]; then
    rate=$(grep -oE '[0-9.]+% overall alignment rate' "$logf" | head -1 | grep -oE '[0-9.]+' || true)
  fi
  if [[ -z "$rate" ]]; then
    echo "WARN: no overall alignment rate in $logf" >&2
    rate="NA"
  fi
  printf 'sample\toverall_alignment_rate_pct\tsource_log\n' > "$rate_tsv"
  printf '%s\t%s\t%s\n' "$sample" "$rate" "$logf" >> "$rate_tsv"
  echo "  hisat2 overall alignment rate: ${rate}%"
}

# Rebuild project-level rate table from per-sample TSVs.
refresh_hisat2_rate_table() {
  local out="08.align_stats/hisat2_alignment_rates.tsv"
  {
    printf 'sample\toverall_alignment_rate_pct\tsource_log\n'
    local f
    for f in 08.align_stats/*_hisat2_rate.tsv; do
      [[ -f "$f" ]] || continue
      # skip header of each per-sample file
      tail -n +2 "$f"
    done | sort -u
  } > "${out}.tmp"
  mv -f "${out}.tmp" "$out"
}

# CSI only (BAI fails when any ref length > 2^29, e.g. wheat). Never delete BAM on index failure.
maybe_index_bam_csi() {
  local bam="$1"
  [[ "$index_bam" -eq 1 ]] || return 0
  if [[ -s "${bam}.csi" ]]; then
    echo "  CSI 已存在，跳过: ${bam}.csi"
    return 0
  fi
  rm -f "${bam}.bai"
  if samtools index -c -@"${sam_threads}" "$bam"; then
    echo "  CSI 索引完成: ${bam}.csi"
  else
    echo "WARN: samtools index -c 失败（BAM 仍可用于 featureCounts）: $bam" >&2
  fi
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

# Keep numbered dirs continuous: 01 → 02 → 03 … (08 used for stats; 04–07 by StringTie)
for dir in 01.ref 02.hisat2-mapping 03.stringtie_sort 08.align_stats; do
  mkdir -p "$dir"
done

find 02.hisat2-mapping 03.stringtie_sort -name '*.tmp' -type f -delete 2>/dev/null || true

genome=$(realpath "$genome")
transdb=$(realpath "$transdb")
ref_dir=$(realpath "01.ref")
genome_basename=$(basename "$genome")
local_index="${ref_dir}/${genome_basename}"

echo "=== MAP-ONLY v${VERSION} ==="
echo "genome=$genome"
echo "fastq=$transdb"
echo "R1_GLOB=$R1_GLOB  R1_SUFFIX=$R1_SUFFIX  R2_SUFFIX=$R2_SUFFIX"
echo "threads=$threads sort_mem=$sort_mem force=$force debug=$debug index_bam=$index_bam"
echo "mode=$([[ "$debug" -eq 1 ]] && echo multi-write-debug || echo onewrite-pipe)"
echo "cwd=$(pwd)"
n_r1=$(grep -c . <<<"$trans_files" || true)
echo "R1 samples: $n_r1"

echo -e "\n=== 步骤1: 参考索引（先验再用，禁止无故重建）==="
hisat_x=""
if [[ -n "$index_prefix" ]]; then
  index_prefix=$(realpath -s "$index_prefix")
  if index_ok "$index_prefix"; then
    hisat_x="$index_prefix"
    echo "使用外部索引: $hisat_x"
    # Keep 01.ref populated for continuity: symlink ht2(l) if missing locally
    ln -sfn "$genome" "${ref_dir}/${genome_basename}"
    for ext in 1.ht2 2.ht2 3.ht2 4.ht2 5.ht2 6.ht2 7.ht2 8.ht2 \
               1.ht2l 2.ht2l 3.ht2l 4.ht2l 5.ht2l 6.ht2l 7.ht2l 8.ht2l; do
      src="${index_prefix}.${ext}"
      dst="${local_index}.${ext}"
      if [[ -f "$src" && ! -e "$dst" ]]; then
        ln -sfn "$src" "$dst"
      fi
    done
  else
    echo "ERROR: -x/--index 无效（缺少 ${index_prefix}.1.ht2 或 .1.ht2l）" >&2
    exit 1
  fi
elif index_ok "$local_index"; then
  hisat_x="$local_index"
  echo "索引已存在，跳过: $hisat_x"
else
  echo "01.ref 无索引，开始 hisat2-build（基因组索引，不含 --ss/--exon）"
  ln -sfn "$genome" "${ref_dir}/${genome_basename}"
  hisat2-build -p "$threads" "${ref_dir}/${genome_basename}" "$local_index"
  hisat_x="$local_index"
fi
index_ok "$hisat_x" || { echo "ERROR: HISAT2 index not usable: $hisat_x" >&2; exit 1; }

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
  hisat_log="02.hisat2-mapping/${sample_name}.hisat2.log"
  sorted_bam="03.stringtie_sort/${sample_name}_sorted.bam"
  stats_raw="08.align_stats/${sample_name}_raw.flagstat.txt"
  stats_unique="08.align_stats/${sample_name}_unique.flagstat.txt"

  if [[ "$force" -eq 0 ]] && bam_ok "$sorted_bam"; then
    echo "已完成，跳过: $sample_name"
    # Backfill rate table from existing log if missing
    if [[ ! -s "08.align_stats/${sample_name}_hisat2_rate.tsv" && -s "$hisat_log" ]]; then
      write_hisat2_rate "$sample_name" "$hisat_log"
    fi
    maybe_index_bam_csi "$sorted_bam"
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
  id2=$(fq_mate_path "$id") || exit 1
  [[ -f "$id2" ]] || { echo "ERROR: missing mate $id2" >&2; exit 1; }

  if [[ "$debug" -eq 1 ]]; then
    # Legacy multi-write path (02 intermediates kept)
    if bam_ok "$raw_bam"; then
      echo "  [debug] 原始BAM已存在，跳过比对"
    else
      rm -f "$raw_bam"
      if [[ -f "$legacy_sam" ]]; then
        samtools view -bS -@"${sam_threads}" "$legacy_sam" > "${raw_bam}.tmp"
        atomic_mv "${raw_bam}.tmp" "$raw_bam"
      else
        echo "  [debug] hisat2 → raw BAM [-p ${threads}, --dta --mm]"
        hisat2 -x "$hisat_x" -p "$threads" $hisat2_extra \
          --dta --mm \
          --fr -1 "$id" -2 "$id2" \
          2> "$hisat_log" \
          | samtools view -@"${sam_threads}" -bS - > "${raw_bam}.tmp"
        atomic_mv "${raw_bam}.tmp" "$raw_bam"
      fi
    fi
    [[ -f "$stats_raw" && -s "$stats_raw" ]] || samtools flagstat -@"${sam_threads}" "$raw_bam" > "$stats_raw"

    if bam_ok "$unique_bam"; then
      echo "  [debug] unique BAM 已存在，跳过筛选"
    else
      rm -f "$unique_bam"
      echo "  [debug] 筛选 unique → BAM"
      samtools view -h -q "$mapq_min" -@"${sam_threads}" "$raw_bam" | \
        awk "$unique_awk" | \
        samtools view -bS -@"${sam_threads}" - > "${unique_bam}.tmp"
      atomic_mv "${unique_bam}.tmp" "$unique_bam"
    fi
    [[ -f "$stats_unique" && -s "$stats_unique" ]] || samtools flagstat -@"${sam_threads}" "$unique_bam" > "$stats_unique"

    echo "  [debug] 排序 → $sorted_bam"
    samtools sort -m "${sort_mem}" -@"${sort_threads}" -o "${sorted_bam}.tmp" "$unique_bam"
    atomic_mv "${sorted_bam}.tmp" "$sorted_bam"
  else
    # Default: one write. Resume unfinished legacy intermediates if present.
    if bam_ok "$raw_bam" || bam_ok "$unique_bam" || [[ -f "$legacy_sam" ]]; then
      echo "  检测到旧式中间文件，走兼容收尾（写完 sorted 后删除中间 BAM）"
      if ! bam_ok "$raw_bam"; then
        rm -f "$raw_bam"
        if [[ -f "$legacy_sam" ]]; then
          samtools view -bS -@"${sam_threads}" "$legacy_sam" > "${raw_bam}.tmp"
          atomic_mv "${raw_bam}.tmp" "$raw_bam"
        else
          hisat2 -x "$hisat_x" -p "$threads" $hisat2_extra --dta --mm --fr -1 "$id" -2 "$id2" \
            2> "$hisat_log" | samtools view -@"${sam_threads}" -bS - > "${raw_bam}.tmp"
          atomic_mv "${raw_bam}.tmp" "$raw_bam"
        fi
      fi
      if ! bam_ok "$unique_bam"; then
        rm -f "$unique_bam"
        samtools view -h -q "$mapq_min" -@"${sam_threads}" "$raw_bam" | \
          awk "$unique_awk" | \
          samtools view -bS -@"${sam_threads}" - > "${unique_bam}.tmp"
        atomic_mv "${unique_bam}.tmp" "$unique_bam"
      fi
      samtools sort -m "${sort_mem}" -@"${sort_threads}" -o "${sorted_bam}.tmp" "$unique_bam"
      atomic_mv "${sorted_bam}.tmp" "$sorted_bam"
      rm -f "$raw_bam" "$unique_bam" "$legacy_sam"
    else
      echo "  hisat2 | unique | sort → $sorted_bam [-p ${threads}, --dta --mm]"
      # shellcheck disable=SC2094
      hisat2 -x "$hisat_x" -p "$threads" $hisat2_extra \
        --dta --mm \
        --fr -1 "$id" -2 "$id2" \
        2> "$hisat_log" \
        | samtools view -h -q "$mapq_min" -@"${sam_threads}" - \
        | awk "$unique_awk" \
        | samtools sort -m "${sort_mem}" -@"${sort_threads}" -o "${sorted_bam}.tmp" -
      pipe_rc=( "${PIPESTATUS[@]}" )
      # hisat2=0, view=1, awk=2, sort=3
      if [[ "${pipe_rc[0]}" -ne 0 || "${pipe_rc[3]}" -ne 0 ]]; then
        rm -f "${sorted_bam}.tmp"
        echo "ERROR: pipe failed for $sample_name (hisat2=${pipe_rc[0]} sort=${pipe_rc[3]}); see $hisat_log" >&2
        exit 1
      fi
      atomic_mv "${sorted_bam}.tmp" "$sorted_bam"
    fi
  fi

  bam_ok "$sorted_bam" || { echo "ERROR: sorted BAM invalid after write: $sorted_bam" >&2; exit 1; }
  # Final BAM is already unique-filtered; keep *_unique.flagstat name for QC compatibility
  samtools flagstat -@"${sam_threads}" "$sorted_bam" > "$stats_unique"
  write_hisat2_rate "$sample_name" "$hisat_log"
  maybe_index_bam_csi "$sorted_bam"

  echo "  完成: $sample_name"
  n_done=$((n_done + 1))
done <<< "$trans_files"

refresh_hisat2_rate_table
echo -e "\n=== MAP-ONLY 完成 ==="
echo "本任务新完成: $n_done  已存在跳过: $n_skip  清理坏文件: $n_fail_clean"
echo "BAM: 03.stringtie_sort/  日志: 02.hisat2-mapping/  统计: 08.align_stats/"
echo "比对率表: 08.align_stats/hisat2_alignment_rates.tsv"
echo "（未跑定量；续跑同一命令即可跳过已完成样本）"
