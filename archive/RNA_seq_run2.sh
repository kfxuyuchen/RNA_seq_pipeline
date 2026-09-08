#!/bin/bash
set -euo pipefail

VERSION="2.0.0"
SCRIPT_NAME=$(basename "$0")

threads=20
debug=0
genome=""
transdb=""
tiegff=""
prepDE_script="/public/home/xuyuchen/local/RNA_seq_pipeline/prepDE.py"
getTPM_script="/public/home/xuyuchen/local/RNA_seq_pipeline/getTPM.py"
getFPKM_script="/public/home/xuyuchen/local/RNA_seq_pipeline/getFPKM.py"

usage() {
  echo "RNA-seq数据分析流程脚本 v${VERSION}"
  echo "Usage: ${SCRIPT_NAME} -G <参考基因组文件> -T <转录组数据文件夹> -t <输入gff文件> [选项]"
  echo "必需参数:"
  echo "  -G, --genome      参考基因组文件（fasta格式）"
  echo "  -T, --transdb     转录组数据文件夹（包含*clean_1.fastq.gz文件）"
  echo "  -t, --tiegff      基因注释GFF文件"
  echo "可选参数:"
  echo "  -p, --threads     线程数 (默认: ${threads})"
  echo "  -d, --debug       调试模式（保留中间BAM/SAM）"
  echo "  -h, --help        显示帮助信息"
  echo "  --prepDE          prepDE.py脚本路径 (默认: ${prepDE_script})"
  echo "  --getTPM          getTPM.py脚本路径 (默认: ${getTPM_script})"
  echo "  --getFPKM         getFPKM.py脚本路径 (默认: ${getFPKM_script})"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -G|--genome) genome="$2"; shift 2 ;;
    -T|--transdb) transdb="$2"; shift 2 ;;
    -t|--tiegff) tiegff="$2"; shift 2 ;;
    -p|--threads) threads="$2"; shift 2 ;;
    -d|--debug) debug=1; shift ;;
    --prepDE) prepDE_script="$2"; shift 2 ;;
    --getTPM) getTPM_script="$2"; shift 2 ;;
    --getFPKM) getFPKM_script="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "错误：未知选项 $1" >&2; usage ;;
  esac
done

if [[ -z "$genome" || -z "$transdb" || -z "$tiegff" ]]; then
  echo "错误：缺少必需参数！" >&2
  usage
fi

check_dependency() {
  if ! command -v "$1" &> /dev/null; then
    echo "错误：依赖工具 $1 未安装或未在PATH中！" >&2
    exit 1
  fi
}

check_file() {
  if [[ ! -f "$1" ]]; then
    echo "错误：文件 $1 不存在！" >&2
    exit 1
  fi
}

check_dir() {
  if [[ ! -d "$1" ]]; then
    echo "错误：目录 $1 不存在！" >&2
    exit 1
  fi
}

sam_threads=$(( threads > 4 ? 4 : threads ))
sort_threads=$(( threads > 4 ? 4 : threads ))
sort_mem="${SAMTOOLS_SORT_MEM:-2G}"
st_parallel=$(( threads / 4 ))
[[ "$st_parallel" -lt 1 ]] && st_parallel=1

echo "=== 环境检查 ==="
check_dependency "hisat2"
check_dependency "hisat2-build"
check_dependency "samtools"
check_dependency "stringtie"
check_dependency "python3"
check_file "$prepDE_script"
check_file "$getTPM_script"
check_file "$getFPKM_script"

check_file "$genome"
check_dir "$transdb"
check_file "$tiegff"

trans_files=$(find "$transdb" -maxdepth 1 -name "*clean_1.fastq.gz" -print | sort)
if [[ -z "$trans_files" ]]; then
  echo "错误：转录组目录 $transdb 中未找到*clean_1.fastq.gz文件！" >&2
  exit 1
fi

echo -e "\n=== 创建工作目录 ==="
workdirs=("01.ref" "02.hisat2-mapping" "03.stringtie_sort"
          "04.stringtie_out" "05.count_out" "06.TPM_out" "07.FPKM_out" "08.align_stats")
for dir in "${workdirs[@]}"; do
  mkdir -p "$dir"
  echo "创建目录: $dir"
done

genome=$(realpath "$genome")
transdb=$(realpath "$transdb")
tiegff=$(realpath "$tiegff")
ref_dir=$(realpath "01.ref")
genome_basename=$(basename "$genome")

# 步骤1: 参考基因组索引构建
echo -e "\n=== 步骤1/8: 参考基因组索引构建 ==="
if [[ ! -f "${ref_dir}/${genome_basename}.1.ht2l" ]]; then
  echo "创建参考基因组软链接..."
  ln -sf "$genome" "${ref_dir}/"
  echo "开始构建hisat2索引..."
  hisat2-build -p "$threads" "${ref_dir}/${genome_basename}" "${ref_dir}/${genome_basename}"
  echo "索引构建完成: ${ref_dir}/${genome_basename}.*.ht2l"
else
  echo "索引文件已存在，跳过构建"
fi

# 步骤2: 逐样本比对→筛选→排序（不落盘SAM，兼容run1遗留SAM）
echo -e "\n=== 步骤2/8: 比对、筛选与排序（pipe BAM） ==="
while IFS= read -r id; do
  [[ -n "$id" ]] || continue
  sample_name=$(basename "$id" "_1.fastq.gz")
  legacy_sam="02.hisat2-mapping/${sample_name}.sam"
  raw_bam="02.hisat2-mapping/${sample_name}_raw.bam"
  unique_bam="02.hisat2-mapping/${sample_name}_unique.bam"
  sorted_bam="03.stringtie_sort/${sample_name}_sorted.bam"
  stats_raw="08.align_stats/${sample_name}_raw.flagstat.txt"
  stats_unique="08.align_stats/${sample_name}_unique.flagstat.txt"

  if [[ -f "$sorted_bam" ]]; then
    echo "已排序，跳过: $sample_name"
    continue
  fi

  echo "处理样本: $sample_name"

  if [[ ! -f "$raw_bam" ]]; then
    if [[ -f "$legacy_sam" ]]; then
      echo "  从run1遗留SAM转BAM: $legacy_sam"
      samtools view -bS -@"${sam_threads}" "$legacy_sam" > "$raw_bam"
    else
      id2="$(dirname "$id")/$(basename "$id" "_1.fastq.gz")_2.fastq.gz"
      if [[ ! -f "$id2" ]]; then
        echo "错误：未找到配对文件 $id2" >&2
        exit 1
      fi
      echo "  hisat2 → samtools view (pipe BAM)"
      hisat2 -x "${ref_dir}/${genome_basename}" -p "$threads" \
        --fr -1 "$id" -2 "$id2" \
        | samtools view -@"${sam_threads}" -bS - > "$raw_bam"
    fi
  else
    echo "  原始BAM已存在，跳过比对"
  fi

  if [[ ! -f "$stats_raw" ]]; then
    echo "  生成原始比对统计: $stats_raw"
    samtools flagstat -@"${sam_threads}" "$raw_bam" > "$stats_raw"
  fi

  if [[ ! -f "$unique_bam" ]]; then
    echo "  筛选unique mapping → BAM"
    samtools view -h -q 10 -@"${sam_threads}" "$raw_bam" | \
      grep -E "^@|NH:i:1$|NH:i:1[^0-9]" | \
      samtools view -bS -@"${sam_threads}" - > "$unique_bam"
  else
    echo "  unique BAM已存在，跳过筛选"
  fi

  if [[ ! -f "$stats_unique" ]]; then
    echo "  生成筛选后比对统计: $stats_unique"
    samtools flagstat -@"${sam_threads}" "$unique_bam" > "$stats_unique"
  fi

  if [[ "$debug" -eq 0 ]]; then
    rm -f "$raw_bam"
    rm -f "$legacy_sam"
  fi

  echo "  排序BAM: $sorted_bam (sort -m ${sort_mem} -@${sort_threads})"
  samtools sort -m "${sort_mem}" -@"${sort_threads}" -o "$sorted_bam" "$unique_bam"

  if [[ "$debug" -eq 0 ]]; then
    rm -f "$unique_bam"
  fi

  echo "  完成: $sample_name"
done <<< "$trans_files"

# 步骤3: Stringtie定量
echo -e "\n=== 步骤3/8: Stringtie表达量计算 ==="
find 03.stringtie_sort -name "*_sorted.bam" -print0 | xargs -0 -I {} -P "${st_parallel}" sh -c '
    id="{}"
    sample_name=$(basename "$id" "_sorted.bam")
    output_dir="04.stringtie_out/${sample_name}"
    mkdir -p "$output_dir"

    if [[ ! -f "${output_dir}/${sample_name}.assemble.gtf" ]]; then
      echo "定量样本: $sample_name"
      stringtie -e -A "${output_dir}/${sample_name}_gene_abund.tab" \
        -B -p 8 -G "$1" -o "${output_dir}/${sample_name}.assemble.gtf" "$id"
    else
      echo "定量结果已存在，跳过: $sample_name"
    fi
  ' sh "$tiegff"

# 步骤4: 生成count矩阵
echo -e "\n=== 步骤4/8: 生成count矩阵 ==="
if [[ ! -f "05.count_out/gene_count.csv" ]]; then
  find 04.stringtie_out -name "*.assemble.gtf" -print0 | while IFS= read -r -d '' gtf; do
    sample_name=$(basename "$(dirname "$gtf")")
    echo -e "${sample_name}\t${gtf}"
  done > sample_lst

  python2 "$prepDE_script" -i sample_lst \
    -g "05.count_out/gene_count.csv" \
    -t "05.count_out/transcript_count.csv"

  if [[ "$debug" -eq 0 ]]; then
    rm -f sample_lst
  fi
else
  echo "count矩阵已存在，跳过生成"
fi

# 步骤5: 生成TPM矩阵
echo -e "\n=== 步骤5/8: 生成TPM矩阵 ==="
if [[ ! -f "06.TPM_out/gene_TPM.csv" ]]; then
  if [[ ! -f "sample_lst" ]]; then
    find 04.stringtie_out -name "*.assemble.gtf" -print0 | while IFS= read -r -d '' gtf; do
      sample_name=$(basename "$(dirname "$gtf")")
      echo -e "${sample_name}\t${gtf}"
    done > sample_lst
  fi

  python2 "$getTPM_script" -i sample_lst \
    -g "06.TPM_out/gene_TPM.csv" \
    -t "06.TPM_out/transcript_TPM.csv"
else
  echo "TPM矩阵已存在，跳过生成"
fi

# 步骤6: 生成FPKM矩阵
echo -e "\n=== 步骤6/8: 生成FPKM矩阵 ==="
if [[ ! -f "07.FPKM_out/gene_FPKM.csv" ]]; then
  if [[ ! -f "sample_lst" ]]; then
    find 04.stringtie_out -name "*.assemble.gtf" -print0 | while IFS= read -r -d '' gtf; do
      sample_name=$(basename "$(dirname "$gtf")")
      echo -e "${sample_name}\t${gtf}"
    done > sample_lst
  fi

  python2 "$getFPKM_script" -i sample_lst \
    -g "07.FPKM_out/gene_FPKM.csv" \
    -t "07.FPKM_out/transcript_FPKM.csv"

  if [[ "$debug" -eq 0 ]]; then
    rm -f sample_lst
  fi
else
  echo "FPKM矩阵已存在，跳过生成"
fi

# 步骤7: 汇总统计信息
echo -e "\n=== 步骤7/8: 汇总比对统计信息 ==="
echo "比对统计信息已保存至: 08.align_stats/"
echo "  - 原始比对统计: *raw.flagstat.txt"
echo "  - 筛选后比对统计: *unique.flagstat.txt"

echo -e "\n=== 步骤8/8: 分析流程完成 ==="
echo "结果文件:"
echo "count矩阵: 05.count_out/"
echo "TPM矩阵: 06.TPM_out/"
echo "FPKM矩阵: 07.FPKM_out/"
echo "比对统计: 08.align_stats/"
echo "运行日志: $(date +'%Y%m%d_%H%M%S').log"
