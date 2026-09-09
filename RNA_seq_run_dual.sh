#!/bin/bash
# Unified RNA-seq pipeline: one mapping, then StringTie and/or featureCounts.
# All stages support resume (skip finished work; --force to redo).
#
#   HISAT2 (--dta --mm) → sorted BAM
#     ├─ StringTie -e → prepDE / getTPM / getFPKM   (04–07.*)
#     └─ featureCounts → count / FPKM / TPM         (09.featureCounts/)
#     └─ QC tables + MultiQC                        (10.qc/)
#
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=pipeline_common.sh
source "${PIPELINE_DIR}/pipeline_common.sh"

VERSION="3.2.0-dual-compat"

threads=12
debug=0
force=0
mode="both"
genome=""
transdb=""
tiegff=""
annot_gtf=""
skip_map=0
run_qc=1

# featureCounts passthrough
fc_strand=0
fc_mapq=10
fc_feature="exon"
fc_attr="gene_id"
fc_skip_preflight=0
fc_min_assigned=50

map_script="${PIPELINE_DIR}/map_only.sh"
post_st="${PIPELINE_DIR}/post_stringtie.sh"
post_fc="${PIPELINE_DIR}/post_featurecounts.sh"
qc_script="${PIPELINE_DIR}/run_qc.sh"

usage() {
  cat <<EOF
RNA-seq dual pipeline v${VERSION}（断点续传 / 兼容增强）

一次比对，可选 StringTie 与/或 featureCounts；结束后默认跑 QC 汇总（含 MultiQC，若已安装）。

Usage:
  $(basename "$0") -G <genome.fa> -T <fastq_dir> [options]
  $(basename "$0") -m quant-only -t <gff> -a <gtf> [fc options]

──────────────── 必需（非 quant-only）────────────────
  -G, --genome PATH      参考基因组 FASTA
  -T, --transdb DIR      PE FASTQ 目录

──────────────── 模式 -m / --mode ────────────────
  map            只比对
  stringtie      比对 + StringTie（需 -t）
  featurecounts  比对 + featureCounts（需 -a）
  both           比对 + 两条定量（默认；需 -t 与 -a）
  quant-only     跳过比对；有 -t 跑 StringTie，有 -a 跑 featureCounts

──────────────── 注释 ────────────────
  -t, --tiegff PATH      StringTie 用 GFF/GFF3
  -a, --annot PATH       featureCounts 用 GTF

──────────────── FASTQ 命名（map 阶段）────────────────
  --r1-suffix STR        默认: ${R1_SUFFIX}
  --r2-suffix STR        默认: ${R2_SUFFIX}
  --r1-glob GLOB         默认: ${R1_GLOB}（实验室 clean 数据）
  例: --r1-suffix _1.fq.gz --r2-suffix _2.fq.gz --r1-glob '*_1.fq.gz'

──────────────── featureCounts 透传 ────────────────
  -s, --strand 0|1|2     链特异（默认: ${fc_strand}）
  -Q, --mapq N           最低 MAPQ（默认: ${fc_mapq}）
  --fc-feature TYPE      特征类型（默认: ${fc_feature}）
  --fc-attr ATTR         基因属性（默认: ${fc_attr}）
  --skip-preflight       跳过 BAM/GTF 染色体与属性检查
  --min-assigned PCT     Assigned% 告警阈值（默认: ${fc_min_assigned}）

──────────────── 其它 ────────────────
  -p, --threads N        线程（默认: ${threads}）
  -d, --debug            比对保留中间 BAM
  --force                强制重跑（向下传递）
  --no-qc                结束后不跑 10.qc / MultiQC
  --qc-only              仅对当前目录已有结果跑 QC（需在输出目录执行）
  -h, --help             本帮助

断点续传: 同一命令重新提交即可；已完成样本/步骤会跳过。

输出目录约定:
  01.ref/  03.stringtie_sort/  04–07.*  08.align_stats/
  09.featureCounts/  10.qc/

示例:
  # 双路径 + 默认 clean FASTQ 命名
  bash $(basename "$0") -G genome.fa -T clean_fastq -t genes.gff3 -a genes.gtf -m both -p 32

  # BAM 已齐，只跑 featureCounts（反向链文库）
  bash $(basename "$0") -m quant-only -a genes.gtf -s 2 -p 32

  # 自定义 FASTQ 后缀
  bash $(basename "$0") -G genome.fa -T fq -a genes.gtf -m featurecounts \\
    --r1-suffix _R1.fastq.gz --r2-suffix _R2.fastq.gz --r1-glob '*_R1.fastq.gz'

  # 只刷新 QC
  cd /path/to/rnaseq_quant_XXX && bash $(basename "$0") --qc-only --title "my project"
EOF
  exit 1
}

qc_title=""
qc_only=0

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
    --r1-suffix) R1_SUFFIX="$2"; shift 2 ;;
    --r2-suffix) R2_SUFFIX="$2"; shift 2 ;;
    --r1-glob) R1_GLOB="$2"; shift 2 ;;
    -s|--strand) fc_strand="$2"; shift 2 ;;
    -Q|--mapq) fc_mapq="$2"; shift 2 ;;
    --fc-feature) fc_feature="$2"; shift 2 ;;
    --fc-attr) fc_attr="$2"; shift 2 ;;
    --skip-preflight) fc_skip_preflight=1; shift ;;
    --min-assigned) fc_min_assigned="$2"; shift 2 ;;
    --no-qc) run_qc=0; shift ;;
    --qc-only) qc_only=1; shift ;;
    --title) qc_title="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "未知选项: $1" >&2; usage ;;
  esac
done

if [[ "$qc_only" -eq 1 ]]; then
  echo "=== QC-only ==="
  qc_args=(-o 10.qc)
  [[ -n "$qc_title" ]] && qc_args+=(--title "$qc_title")
  bash "$qc_script" "${qc_args[@]}"
  exit 0
fi

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
echo "st=$need_st  fc=$need_fc  skip_map=$skip_map  qc=$run_qc"
echo "FQ: R1_GLOB=$R1_GLOB R1_SUFFIX=$R1_SUFFIX R2_SUFFIX=$R2_SUFFIX"
if [[ "$need_fc" -eq 1 ]]; then
  echo "FC: strand=$fc_strand mapq=$fc_mapq feature=$fc_feature attr=$fc_attr"
fi
echo "Start: $(date)"
echo "========================================"

if [[ "$skip_map" -eq 0 ]]; then
  map_args=(-G "$genome" -T "$transdb" -p "$threads"
            --r1-suffix "$R1_SUFFIX" --r2-suffix "$R2_SUFFIX" --r1-glob "$R1_GLOB")
  [[ -n "$tiegff" ]] && map_args+=(-t "$tiegff")
  [[ "$debug" -eq 1 ]] && map_args+=(-d)
  [[ "$force" -eq 1 ]] && map_args+=(--force)
  bash "$map_script" "${map_args[@]}"
fi

if [[ "$mode" == "map" ]]; then
  echo "mode=map：到此结束。BAM 在 03.stringtie_sort/"
  if [[ "$run_qc" -eq 1 ]]; then
    bash "$qc_script" -o 10.qc --title "${qc_title:-RNA-seq map-only QC}" || true
  fi
  exit 0
fi

if [[ "$need_st" -eq 1 ]]; then
  echo -e "\n########## StringTie 定量 ##########"
  st_args=(-t "$tiegff" -p "$threads")
  [[ "$force" -eq 1 ]] && st_args+=(--force)
  bash "$post_st" "${st_args[@]}"
fi

if [[ "$need_fc" -eq 1 ]]; then
  echo -e "\n########## featureCounts 定量 ##########"
  fc_args=(-a "$annot_gtf" -p "$threads"
           -s "$fc_strand" -Q "$fc_mapq"
           -t "$fc_feature" -g "$fc_attr"
           --min-assigned "$fc_min_assigned")
  [[ "$force" -eq 1 ]] && fc_args+=(--force)
  [[ "$fc_skip_preflight" -eq 1 ]] && fc_args+=(--skip-preflight)
  bash "$post_fc" "${fc_args[@]}"
fi

if [[ "$run_qc" -eq 1 ]]; then
  echo -e "\n########## QC / MultiQC ##########"
  qc_args=(-o 10.qc --min-assigned "$fc_min_assigned")
  [[ -n "$qc_title" ]] && qc_args+=(--title "$qc_title")
  bash "$qc_script" "${qc_args[@]}" || echo "WARN: QC step reported issues (see above)"
fi

echo -e "\n========================================"
echo "全部完成: $(date)"
echo "StringTie: 04.stringtie_out/ 05.count_out/ 06.TPM_out/ 07.FPKM_out/"
echo "featureCounts: 09.featureCounts/"
echo "QC: 10.qc/"
echo "========================================"
