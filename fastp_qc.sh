#!/bin/bash
# fastp QC stage: raw FASTQ → clean FASTQ with the lab-standard parameters.
# PE: <sample>_clean_1/2.fastq.gz
# SE: <sample>_clean.fastq.gz
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=pipeline_common.sh
source "${PIPELINE_DIR}/pipeline_common.sh"

VERSION="1.2.0"
FASTP="${FASTP:-fastp}"

rawdir=""
outdir="00.clean_fastq"
reportdir="10.qc/fastp"
threads=8
min_qual=20
min_len=40
compression=4
detect_adapter=1
remove_raw=0
force=0
FP_R1_GLOB="*_1.fastq.gz"
FP_R2_SUFFIX="_2.fastq.gz"
FP_SE_GLOB="*.fastq.gz"

usage() {
  cat <<USAGE
fastp QC stage v${VERSION}

Usage:
  $(basename "$0") -T <raw_fastq_dir> [options]

必需:
  -T, --transdb DIR      raw FASTQ 目录（PE 或 SE）

可选:
  -o, --outdir DIR       clean FASTQ 输出目录（默认: ${outdir}）
  -j, --report-dir DIR   fastp 报告目录（默认: ${reportdir}）
  -p, --threads N        fastp 线程数（默认: ${threads}）
  --min-qual N           合格碱基 Phred 阈值（默认: ${min_qual}）
  --min-len N            过滤后最短 read 长度（默认: ${min_len}）
  --compression N        gzip 压缩级别（默认: ${compression}）
  --no-detect-adapter    关闭 PE 自动接头检测
  --remove-raw           fastp 成功后删除 raw FASTQ（默认保留）
  --r1-glob GLOB         raw R1 的 glob（默认: ${FP_R1_GLOB}）
  --r2-suffix STR        raw R2 后缀（默认: ${FP_R2_SUFFIX}）
  --se-glob GLOB         单端 raw FASTQ 的 glob（默认: ${FP_SE_GLOB}）
  --force                已存在 clean 文件也重跑
  -h, --help             本帮助

固定参数（与实验室 252 样本项目一致）:
  --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 40 --compression 4
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -T|--transdb) rawdir="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    -j|--report-dir) reportdir="$2"; shift 2 ;;
    -p|--threads) threads="$2"; shift 2 ;;
    --min-qual) min_qual="$2"; shift 2 ;;
    --min-len) min_len="$2"; shift 2 ;;
    --compression) compression="$2"; shift 2 ;;
    --no-detect-adapter) detect_adapter=0; shift ;;
    --remove-raw) remove_raw=1; shift ;;
    --r1-glob) FP_R1_GLOB="$2"; shift 2 ;;
    --r2-suffix) FP_R2_SUFFIX="$2"; shift 2 ;;
    --se-glob) FP_SE_GLOB="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage ;;
    *) echo "未知选项: $1" >&2; usage ;;
  esac
done

[[ -n "$rawdir" && -d "$rawdir" ]] || { echo "ERROR: -T 目录不存在: $rawdir" >&2; exit 1; }
command -v "$FASTP" >/dev/null || { echo "ERROR: fastp not found: $FASTP" >&2; exit 1; }

mkdir -p "$outdir" "$reportdir" "${reportdir}/unpaired"

echo "========================================"
echo "fastp QC v${VERSION}"
echo "raw=$rawdir  clean=$outdir  reports=$reportdir"
echo "threads=$threads  min_qual=$min_qual  min_len=$min_len  compression=$compression"
echo "detect_adapter=$detect_adapter  remove_raw=$remove_raw  force=$force"
echo "Start: $(date)"
echo "========================================"

n=0; n_skip=0; n_fail=0

while IFS= read -r r1; do
  [[ -n "$r1" ]] || continue
  base=$(basename "$r1")

  if [[ "$base" == *"${FP_R1_GLOB#\*}" ]]; then
    sample="${base%${FP_R1_GLOB#\*}}"
    r2="${rawdir}/${sample}${FP_R2_SUFFIX}"
    out1="${outdir}/${sample}_clean_1.fastq.gz"
    out2="${outdir}/${sample}_clean_2.fastq.gz"
    is_pe=1
  elif [[ "$base" == *.fastq.gz ]]; then
    sample="${base%.fastq.gz}"
    r2=""
    out1="${outdir}/${sample}_clean.fastq.gz"
    out2=""
    is_pe=0
  else
    continue
  fi

  unpaired="${reportdir}/unpaired/${sample}.unpaired.fq.gz"
  json="${reportdir}/${sample}.fastp.json"
  html="${reportdir}/${sample}.fastp.html"
  log="${reportdir}/${sample}.fastp.log"

  if [[ "$is_pe" -eq 1 ]]; then
    if [[ "$force" -eq 0 && -s "$out1" && -s "$out2" ]]; then
      n_skip=$((n_skip + 1)); echo "[SKIP] ${sample} (PE)"; continue
    fi
    if [[ ! -f "$r2" ]]; then
      echo "[FAIL] ${sample}: R2 missing ($r2)" >&2
      n_fail=$((n_fail + 1)); continue
    fi
    fp_args=(--in1 "$r1" --in2 "$r2" --out1 "$out1" --out2 "$out2"
             --unpaired1 "$unpaired" --unpaired2 "$unpaired"
             --thread "$threads" --compression "$compression"
             --qualified_quality_phred "$min_qual" --length_required "$min_len"
             --json "$json" --html "$html" --report_title "$sample")
    [[ "$detect_adapter" -eq 1 ]] && fp_args+=(--detect_adapter_for_pe)
  else
    if [[ "$force" -eq 0 && -s "$out1" ]]; then
      n_skip=$((n_skip + 1)); echo "[SKIP] ${sample} (SE)"; continue
    fi
    fp_args=(--in1 "$r1" --out1 "$out1"
             --thread "$threads" --compression "$compression"
             --qualified_quality_phred "$min_qual" --length_required "$min_len"
             --json "$json" --html "$html" --report_title "$sample")
  fi

  echo "[START] ${sample}"
  if "$FASTP" "${fp_args[@]}" > "$log" 2>&1; then
    n=$((n + 1)); echo "[DONE] ${sample}"
    if [[ "$remove_raw" -eq 1 ]]; then
      rm -f "$r1"
      [[ "$is_pe" -eq 1 ]] && rm -f "$r2"
    fi
  else
    n_fail=$((n_fail + 1)); echo "[FAIL] ${sample}: see $log" >&2
  fi
done < <({
  find "$rawdir" -maxdepth 1 -name "$FP_R1_GLOB" -print
  find "$rawdir" -maxdepth 1 -name "$FP_SE_GLOB" -print | grep -vE '_1\.fastq\.gz$|_2\.fastq\.gz$' || true
} | sort -u)

echo "========================================"
echo "fastp finished: done=$n skip=$n_skip fail=$n_fail"
echo "End: $(date)"
echo "========================================"
[[ "$n_fail" -eq 0 ]]
