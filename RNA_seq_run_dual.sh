#!/bin/bash
# Unified RNA-seq pipeline: one mapping, then StringTie and/or featureCounts.
# All stages support resume (skip finished work; --force to redo).
#
#   HISAT2 (--dta --mm) | unique | sort → 03.*_sorted.bam (default one-write)
#     ├─ 02.*.hisat2.log
#     ├─ StringTie -e → prepDE / getTPM / getFPKM   (04–07.*)
#     ├─ featureCounts → count / FPKM / TPM         (09.featureCounts/)
#     └─ QC tables + MultiQC                        (10.qc/)
#
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=pipeline_common.sh
source "${PIPELINE_DIR}/pipeline_common.sh"

VERSION="5.1.0"

threads=12
debug=0
force=0
mode="both"
genome=""
transdb=""
tiegff=""
annot_gtf=""
index_prefix=""
index_bam=0
skip_map=0
run_qc=1
min_hisat2=0

# integrated fastp QC (optional)
use_fastp=0
fastp_dir=""
fastp_report_dir=""
fastp_threads=8
fastp_min_qual=20
fastp_min_len=40
fastp_compression=4
fastp_detect_adapter=1
fastp_remove_raw=0
fastp_script="${PIPELINE_DIR}/fastp_qc.sh"

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

──────────────── 内置 fastp（可选）────────────────
  --fastp                输入 -T 为 raw FASTQ，先跑 fastp（固定实验室参数）再比对
  --fastp-dir DIR        clean FASTQ 输出目录（默认: <cwd>/00.clean_fastq）
  --fastp-report-dir DIR fastp 报告目录（默认: <cwd>/10.qc/fastp）
  --fastp-threads N      fastp 线程（默认: 8）
  --fastp-min-qual N     合格碱基 Phred（默认: 20）
  --fastp-min-len N      最短 read 长度（默认: 40）
  --fastp-compression N  gzip 级别（默认: 4）
  --no-detect-adapter    关闭 PE 自动接头检测
  --fastp-remove-raw     fastp 成功后删除 raw FASTQ（默认保留）

──────────────── featureCounts 透传 ────────────────
  -s, --strand 0|1|2     链特异（默认: ${fc_strand}）
  -Q, --mapq N           最低 MAPQ（默认: ${fc_mapq}）
  --fc-feature TYPE      特征类型（默认: ${fc_feature}）
  --fc-attr ATTR         基因属性（默认: ${fc_attr}）
  --skip-preflight       跳过 BAM/GTF 染色体与属性检查
  --min-assigned PCT     Assigned% 告警阈值（默认: ${fc_min_assigned}）

──────────────── 其它 ────────────────
  -p, --threads N        线程（默认: ${threads}）
  -d, --debug            比对多步落盘并保留 02 中间 BAM（默认一次落盘）
  -x, --index PREFIX     复用已有 HISAT2 索引前缀（传给 map_only）
  --index-bam            对 sorted BAM 建 CSI 索引（小麦推荐）
  --force                强制重跑（向下传递）
  --no-qc                结束后不跑 10.qc / MultiQC
  --qc-only              仅对当前目录已有结果跑 QC（需在输出目录执行）
  --min-hisat2 PCT       QC 时 hisat2 overall% 低于此值告警（默认 0=关闭）
  -h, --help             本帮助

断点续传: 同一命令重新提交即可；已完成样本/步骤会跳过。

输出目录约定（编号连续）:
  01.ref/  02.hisat2-mapping/  03.stringtie_sort/  04–07.*
  08.align_stats/  09.featureCounts/  10.qc/

示例:
  # 双路径 + 默认 clean FASTQ 命名
  bash $(basename "$0") -G genome.fa -T clean_fastq -t genes.gff3 -a genes.gtf -m both -p 32

  # BAM 已齐，只跑 featureCounts（反向链文库）
  bash $(basename "$0") -m quant-only -a genes.gtf -s 2 -p 32

  # 从 raw FASTQ 开始，内置 fastp
  bash $(basename "$0") -G genome.fa -T raw_fastq --fastp -a genes.gtf -m featurecounts -p 32

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
    -x|--index) index_prefix="$2"; shift 2 ;;
    --index-bam) index_bam=1; shift ;;
    --force) force=1; shift ;;
    --min-hisat2) min_hisat2="$2"; shift 2 ;;
    --r1-suffix) R1_SUFFIX="$2"; shift 2 ;;
    --r2-suffix) R2_SUFFIX="$2"; shift 2 ;;
    --r1-glob) R1_GLOB="$2"; shift 2 ;;
    --fastp) use_fastp=1; shift ;;
    --fastp-dir) fastp_dir="$2"; shift 2 ;;
    --fastp-report-dir) fastp_report_dir="$2"; shift 2 ;;
    --fastp-threads) fastp_threads="$2"; shift 2 ;;
    --fastp-min-qual) fastp_min_qual="$2"; shift 2 ;;
    --fastp-min-len) fastp_min_len="$2"; shift 2 ;;
    --fastp-compression) fastp_compression="$2"; shift 2 ;;
    --no-detect-adapter) fastp_detect_adapter=0; shift ;;
    --fastp-remove-raw) fastp_remove_raw=1; shift ;;
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
  qc_args=(-o 10.qc --min-assigned "$fc_min_assigned" --min-hisat2 "$min_hisat2")
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
[[ "$use_fastp" -eq 1 ]] && echo "fastp: ON  threads=$fastp_threads  min_qual=$fastp_min_qual  min_len=$fastp_min_len"
if [[ "$need_fc" -eq 1 ]]; then
  echo "FC: strand=$fc_strand mapq=$fc_mapq feature=$fc_feature attr=$fc_attr"
fi
echo "Start: $(date)"
echo "========================================"

if [[ "$skip_map" -eq 0 ]]; then
  n_clean=$(find "$transdb" -maxdepth 1 -name '*clean_1.fastq.gz' | wc -l)
  n_raw=$(find "$transdb" -maxdepth 1 -name '*_1.fastq.gz' ! -name '*clean_1.fastq.gz' | wc -l)
  if [[ "$use_fastp" -eq 1 ]]; then
    echo "INPUT mode: raw FASTQ → fastp → clean（找到 ${n_raw} 个 raw R1）"
    if [[ "$n_raw" -eq 0 && "$n_clean" -gt 0 ]]; then
      echo "WARN: 已加 --fastp，但目录里只有 clean FASTQ；会对 clean 再 trim 一次。" >&2
    fi
  else
    echo "INPUT mode: clean FASTQ（fastp 关闭；找到 ${n_clean} 个 *clean_1.fastq.gz）"
    if [[ "$n_clean" -eq 0 && "$n_raw" -gt 0 ]]; then
      echo "WARN: 未加 --fastp，但目录里只有 raw FASTQ；如果这是 raw 数据，请加 --fastp。" >&2
    fi
  fi
fi

if [[ "$use_fastp" -eq 1 && "$skip_map" -eq 0 ]]; then
  if [[ ! -x "$fastp_script" && ! -f "$fastp_script" ]]; then
    echo "ERROR: missing fastp_qc.sh: $fastp_script" >&2; exit 1
  fi
  fastp_clean="${fastp_dir:-$(pwd)/00.clean_fastq}"
  fastp_rep="${fastp_report_dir:-$(pwd)/10.qc/fastp}"
  fp_args=(-T "$transdb" -o "$fastp_clean" -j "$fastp_rep" -p "$fastp_threads"
           --min-qual "$fastp_min_qual" --min-len "$fastp_min_len"
           --compression "$fastp_compression")
  [[ "$fastp_detect_adapter" -eq 0 ]] && fp_args+=(--no-detect-adapter)
  [[ "$fastp_remove_raw" -eq 1 ]] && fp_args+=(--remove-raw)
  [[ "$force" -eq 1 ]] && fp_args+=(--force)
  echo -e "\n########## fastp QC ##########"
  bash "$fastp_script" "${fp_args[@]}"
  transdb="$fastp_clean"
  echo "fastp done → transdb=$transdb"
fi

if [[ "$skip_map" -eq 0 ]]; then
  map_args=(-G "$genome" -T "$transdb" -p "$threads"
            --r1-suffix "$R1_SUFFIX" --r2-suffix "$R2_SUFFIX" --r1-glob "$R1_GLOB")
  [[ -n "$tiegff" ]] && map_args+=(-t "$tiegff")
  [[ -n "$index_prefix" ]] && map_args+=(-x "$index_prefix")
  [[ "$debug" -eq 1 ]] && map_args+=(-d)
  [[ "$index_bam" -eq 1 ]] && map_args+=(--index-bam)
  [[ "$force" -eq 1 ]] && map_args+=(--force)
  bash "$map_script" "${map_args[@]}"
fi

if [[ "$mode" == "map" ]]; then
  echo "mode=map：到此结束。BAM=03.stringtie_sort/ 日志=02.hisat2-mapping/"
  if [[ "$run_qc" -eq 1 ]]; then
    bash "$qc_script" -o 10.qc --min-assigned "$fc_min_assigned" --min-hisat2 "$min_hisat2" \
      --title "${qc_title:-RNA-seq map-only QC}" || true
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
  qc_args=(-o 10.qc --min-assigned "$fc_min_assigned" --min-hisat2 "$min_hisat2")
  [[ -n "$qc_title" ]] && qc_args+=(--title "$qc_title")
  bash "$qc_script" "${qc_args[@]}" || echo "WARN: QC step reported issues (see above)"
fi

echo -e "\n========================================"
echo "全部完成: $(date)"
echo "StringTie: 04.stringtie_out/ 05.count_out/ 06.TPM_out/ 07.FPKM_out/"
echo "featureCounts: 09.featureCounts/"
echo "QC: 10.qc/"
echo "========================================"
