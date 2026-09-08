#!/bin/bash
# Unified RNA-seq pipeline: one mapping, then StringTie and/or featureCounts.
# All stages support resume (skip finished work; --force to redo).
#
#   HISAT2 (--dta --mm) → sorted BAM
#     ├─ StringTie -e → prepDE / getTPM / getFPKM   (04–07.*)
#     └─ featureCounts → count / FPKM / TPM         (09.featureCounts/)
#
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERSION="3.1.0-dual-resume"

threads=12
debug=0
force=0
mode="both"
genome=""
transdb=""
tiegff=""
annot_gtf=""
skip_map=0

map_script="${PIPELINE_DIR}/map_only.sh"
post_st="${PIPELINE_DIR}/post_stringtie.sh"
post_fc="${PIPELINE_DIR}/post_featurecounts.sh"

usage() {
  cat <<EOF
RNA-seq dual pipeline v${VERSION} (断点续传)
Usage: $(basename "$0") -G <genome.fa> -T <fastq_dir> [options]

必需:
  -G, --genome       参考基因组 fasta
  -T, --transdb      PE fastq 目录（*clean_1.fastq.gz）；quant-only 时可省略

模式 (-m):
  map | stringtie | featurecounts | both | quant-only   (默认: both)

注释:
  -t, --tiegff       StringTie 用 GFF/GFF3
  -a, --annot        featureCounts 用 GTF

其它:
  -p, --threads      线程数 (默认: ${threads})
  -d, --debug        比对保留中间 BAM
  --force            强制重跑（向下传递给 map/post）
  -h, --help

断点续传: 直接重新提交同一命令即可；已完成样本/步骤会自动跳过。
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -G|--genome) genome="$2"; shift 2 ;;
    -T|--transdb) transdb="$2"; shift 2 ;;
    -t|--tiegff) tiegff="$2"; shift 2 ;;
    -a|--annot) annot_gtf="$2"; shift 2 ;;
    -m|--mode) mode="$2"; shift 2 ;;
    -p|--threads) threads="$2"; shift 2 ;;
    -d|--debug) debug=1; shift ;;
    --force) force=1; shift ;;
    -h|--help) usage ;;
    *) echo "未知选项: $1" >&2; usage ;;
  esac
done

case "$mode" in
  map|stringtie|featurecounts|both|quant-only) ;;
  *) echo "ERROR: invalid -m $mode" >&2; usage ;;
esac

need_st=0
need_fc=0
case "$mode" in
  stringtie) need_st=1 ;;
  featurecounts) need_fc=1 ;;
  both) need_st=1; need_fc=1 ;;
  quant-only)
    [[ -n "$tiegff" ]] && need_st=1
    [[ -n "$annot_gtf" ]] && need_fc=1
    [[ "$need_st" -eq 1 || "$need_fc" -eq 1 ]] || {
      echo "ERROR: quant-only 至少提供 -t 或 -a" >&2; exit 1
    }
    skip_map=1
    ;;
  map) ;;
esac

[[ "$mode" != "quant-only" ]] && [[ -z "$genome" || -z "$transdb" ]] && usage
[[ "$mode" == "quant-only" ]] || { [[ -f "$genome" && -d "$transdb" ]]; }
[[ "$need_st" -eq 1 && ( -z "$tiegff" || ! -f "$tiegff" ) ]] && {
  echo "ERROR: StringTie 路线需要有效 -t GFF" >&2; exit 1
}
[[ "$need_fc" -eq 1 && ( -z "$annot_gtf" || ! -f "$annot_gtf" ) ]] && {
  echo "ERROR: featureCounts 路线需要有效 -a GTF" >&2; exit 1
}

echo "========================================"
echo "RNA-seq dual v${VERSION}"
echo "mode=$mode  threads=$threads  force=$force  cwd=$(pwd)"
echo "st=$need_st  fc=$need_fc  skip_map=$skip_map"
echo "Start: $(date)"
echo "========================================"

if [[ "$skip_map" -eq 0 ]]; then
  map_args=(-G "$genome" -T "$transdb" -p "$threads")
  [[ -n "$tiegff" ]] && map_args+=(-t "$tiegff")
  [[ "$debug" -eq 1 ]] && map_args+=(-d)
  [[ "$force" -eq 1 ]] && map_args+=(--force)
  bash "$map_script" "${map_args[@]}"
fi

if [[ "$mode" == "map" ]]; then
  echo "mode=map：到此结束。BAM 在 03.stringtie_sort/"
  exit 0
fi

if [[ "$need_st" -eq 1 ]]; then
  echo -e "\n########## StringTie 定量 ##########"
  if [[ "$force" -eq 1 ]]; then
    bash "$post_st" -t "$tiegff" -p "$threads" --force
  else
    bash "$post_st" -t "$tiegff" -p "$threads"
  fi
fi

if [[ "$need_fc" -eq 1 ]]; then
  echo -e "\n########## featureCounts 定量 ##########"
  if [[ "$force" -eq 1 ]]; then
    bash "$post_fc" -a "$annot_gtf" -p "$threads" --force
  else
    bash "$post_fc" -a "$annot_gtf" -p "$threads"
  fi
fi

echo -e "\n========================================"
echo "全部完成: $(date)"
echo "StringTie: 04.stringtie_out/ 05.count_out/ 06.TPM_out/ 07.FPKM_out/"
echo "featureCounts: 09.featureCounts/"
echo "========================================"
