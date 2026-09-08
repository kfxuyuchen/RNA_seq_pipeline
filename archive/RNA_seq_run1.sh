#!/bin/bash
set -euo pipefail  # 严格模式：未定义变量/命令失败/管道错误时退出

# 版本信息
VERSION="1.0.1"  # 更新版本号
SCRIPT_NAME=$(basename "$0")

# 默认参数
threads=20
debug=0
genome=""
transdb=""
tiegff=""
prepDE_script="/public/home/xuyuchen/local/RNA_seq_pipeline/prepDE.py"
getTPM_script="/public/home/xuyuchen/local/RNA_seq_pipeline/getTPM.py"
getFPKM_script="/public/home/xuyuchen/local/RNA_seq_pipeline/getFPKM.py"

# 帮助信息
usage() {
  echo "RNA-seq数据分析流程脚本 v${VERSION}"
  echo "Usage: ${SCRIPT_NAME} -G <参考基因组文件> -T <转录组数据文件夹> -t <输入gff文件> [选项]"
  echo "必需参数:"
  echo "  -G, --genome      参考基因组文件（fasta格式）"
  echo "  -T, --transdb     转录组数据文件夹（包含*clean_1.fastq.gz文件）"
  echo "  -t, --tiegff      基因注释GFF文件"
  echo "可选参数:"
  echo "  -p, --threads     线程数 (默认: ${threads})"
  echo "  -d, --debug       调试模式（保留中间文件）"
  echo "  -h, --help        显示帮助信息"
  echo "  --prepDE          prepDE.py脚本路径 (默认: ${prepDE_script})"
  echo "  --getTPM          getTPM.py脚本路径 (默认: ${getTPM_script})"
  echo "  --getFPKM         getFPKM.py脚本路径 (默认: ${getFPKM_script})"
  exit 1
}

# 参数解析
while [[ $# -gt 0 ]]; do
  case "$1" in
    -G|--genome)
      genome="$2"
      shift 2
      ;;
    -T|--transdb)
      transdb="$2"
      shift 2
      ;;
    -t|--tiegff)
      tiegff="$2"
      shift 2
      ;;
    -p|--threads)
      threads="$2"
      shift 2
      ;;
    -d|--debug)
      debug=1
      shift
      ;;
    --prepDE)
      prepDE_script="$2"
      shift 2
      ;;
    --getTPM)
      getTPM_script="$2"
      shift 2
      ;;
    --getFPKM)
      getFPKM_script="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "错误：未知选项 $1" >&2
      usage
      ;;
  esac
done

# 必需参数校验
if [[ -z "$genome" || -z "$transdb" || -z "$tiegff" ]]; then
  echo "错误：缺少必需参数！" >&2
  usage
fi

# 环境检查函数
check_dependency() {
  if ! command -v "$1" &> /dev/null; then
    echo "错误：依赖工具 $1 未安装或未在PATH中！" >&2
    exit 1
  fi
}

# 文件检查函数
check_file() {
  if [[ ! -f "$1" ]]; then
    echo "错误：文件 $1 不存在！" >&2
    exit 1
  fi
}

# 目录检查函数
check_dir() {
  if [[ ! -d "$1" ]]; then
    echo "错误：目录 $1 不存在！" >&2
    exit 1
  fi
}

# 依赖工具检查
echo "=== 环境检查 ==="
check_dependency "hisat2"
check_dependency "hisat2-build"
check_dependency "samtools"
check_dependency "stringtie"
check_dependency "python3"
check_file "$prepDE_script"
check_file "$getTPM_script"
check_file "$getFPKM_script"

# 输入文件检查
check_file "$genome"
check_dir "$transdb"
check_file "$tiegff"

# 转录组数据检查
trans_files=$(find "$transdb" -maxdepth 1 -name "*clean_1.fastq.gz" -print)
if [[ -z "$trans_files" ]]; then
  echo "错误：转录组目录 $transdb 中未找到*clean_1.fastq.gz文件！" >&2
  exit 1
fi

# 创建工作目录（新增统计目录）
echo -e "\n=== 创建工作目录 ==="
workdirs=("01.ref" "02.hisat2-mapping" "03.stringtie_sort" 
          "04.stringtie_out" "05.count_out" "06.TPM_out" "07.FPKM_out" "08.align_stats")  # 新增统计目录
for dir in "${workdirs[@]}"; do
  mkdir -p "$dir"
  echo "创建目录: $dir"
done

# 路径标准化
genome=$(realpath "$genome")
transdb=$(realpath "$transdb")
tiegff=$(realpath "$tiegff")
ref_dir=$(realpath "01.ref")
genome_basename=$(basename "$genome")

# 步骤1: 参考基因组索引构建
echo -e "\n=== 步骤1/10: 参考基因组索引构建 ==="  # 更新步骤计数
if [[ ! -f "${ref_dir}/${genome_basename}.1.ht2l" ]]; then
  echo "创建参考基因组软链接..."
  ln -sf "$genome" "${ref_dir}/"
  
  echo "开始构建hisat2索引..."
  hisat2-build -p "$threads" "${ref_dir}/${genome_basename}" "${ref_dir}/${genome_basename}"
  echo "索引构建完成: ${ref_dir}/${genome_basename}.*.ht2l"
else
  echo "索引文件已存在，跳过构建"
fi

# 步骤2: 转录组数据比对（新增原始比对统计）
echo -e "\n=== 步骤2/10: 转录组数据比对 ==="  # 更新步骤计数
echo "$trans_files" | while IFS= read -r id; do
  sample_name=$(basename "$id" "_1.fastq.gz")  
  output_sam="02.hisat2-mapping/${sample_name}.sam"
  stats_raw="08.align_stats/${sample_name}_raw.flagstat.txt"  # 原始比对统计文件
  
  if [[ ! -f "$output_sam" ]]; then
    echo "处理样本: $sample_name"
    id2=$(dirname "$id")/$(basename "$id" "_1.fastq.gz")_2.fastq.gz
    
    if [[ ! -f "$id2" ]]; then
      echo "错误：未找到配对文件 $id2" >&2
      exit 1
    fi
    
    # --mm: memory-map index to reduce RSS (important for large genomes e.g. wheat)
    hisat2 -x "${ref_dir}/${genome_basename}" -p "$threads" \
      --mm \
      --fr -1 "$id" -2 "$id2" -S "$output_sam"
    echo "生成SAM文件: $output_sam"
    
    # 新增：生成原始比对统计信息
    echo "生成原始比对统计: $stats_raw"
    samtools flagstat -@4 "$output_sam" > "$stats_raw"
  else
    echo "SAM文件已存在，跳过比对: $output_sam"
    # 补全缺失的统计文件
    if [[ ! -f "$stats_raw" ]]; then
      echo "补全原始比对统计: $stats_raw"
      samtools flagstat -@4 "$output_sam" > "$stats_raw"
    fi
  fi
done

# 步骤3: 筛选unique mapping（核心修复）
echo -e "\n=== 步骤3/10: 筛选unique mapping ==="  # 更新步骤计数
find 02.hisat2-mapping -name "*_clean.sam" -print0 |xargs -0 -I {} -P "$((threads / 4))" sh -c '
    id="{}"
    sample_name=$(basename "$id" ".sam")
    output_sam="${id%.sam}_unique.sam"
    
    if [[ ! -f "$output_sam" ]]; then
      echo "筛选样本: $sample_name"
      # 修复MAPQ阈值参数（使用固定值10），移除错误的$2传递
      samtools view -h -q 10 -@ 4 "$id" | \
        grep -E "^@|NH:i:1$|NH:i:1[^0-9]" > "$output_sam"
    else
      echo "已筛选，跳过: $output_sam"
    fi
  '

# 步骤4: SAM转BAM（新增筛选后比对统计）
echo -e "\n=== 步骤4/10: SAM转BAM ==="  # 更新步骤计数
find 02.hisat2-mapping -name "*_unique.sam" -print0 | xargs -0 -I {} -P "$((threads / 4))" sh -c '
    id="{}"
    output_bam="${id%.sam}.bam"
    sample_name=$(basename "$output_bam" "_unique.bam")
    stats_unique="08.align_stats/${sample_name}_unique.flagstat.txt"  # 筛选后比对统计文件
    
    if [[ ! -f "$output_bam" ]]; then
      echo "转换样本: $sample_name"
      samtools view -bS -@ 4 "$id" > "$output_bam"
      
      # 新增：生成筛选后比对统计信息
      echo "生成筛选后比对统计: $stats_unique"
      samtools flagstat -@4 "$output_bam" > "$stats_unique"
      
      # 非调试模式删除SAM文件
      if [[ "$1" -eq 0 ]]; then
        rm -f "$id"
      fi
    else
      echo "BAM文件已存在，跳过: $output_bam"
      # 补全缺失的统计文件
      if [[ ! -f "$stats_unique" ]]; then
        echo "补全筛选后比对统计: $stats_unique"
        samtools flagstat -@4 "$output_bam" > "$stats_unique"
      fi
    fi
  ' _ "$debug"

# 步骤5: BAM排序
echo -e "\n=== 步骤5/10: BAM文件排序 ==="  # 更新步骤计数
find 02.hisat2-mapping -name "*_unique.bam" -print0 | xargs -0 -I {} -P "$((threads / 4))" sh -c '
    id="{}"
    sample_name=$(basename "$id" "_unique.bam")
    output_bam="03.stringtie_sort/${sample_name}_sorted.bam"
    
    if [[ ! -f "$output_bam" ]]; then  
      echo "排序样本: $sample_name"
      samtools sort -@ 8 -o "$output_bam" "$id"
      
      # 非调试模式删除原始BAM
      if [[ "$1" -eq 0 ]]; then
        rm -f "$id"
      fi
    else  
      echo "排序BAM已存在，跳过: $output_bam"; 
    fi
  ' _ "$debug"

# 步骤6: Stringtie定量
echo -e "\n=== 步骤6/10: Stringtie表达量计算 ===$tiegff"  # 更新步骤计数
find 03.stringtie_sort -name "*_sorted.bam" -print0 | xargs -0 -I {} -P "$((threads / 8))" sh -c '
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

# 步骤7: 生成count矩阵
echo -e "\n=== 步骤7/10: 生成count矩阵 ==="  # 更新步骤计数
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

# 步骤8: 生成TPM矩阵
echo -e "\n=== 步骤8/10: 生成TPM矩阵 ==="  # 更新步骤计数
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

# 步骤9: 生成FPKM矩阵
echo -e "\n=== 步骤9/10: 生成FPKM矩阵 ==="  # 更新步骤计数
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

# 步骤10: 汇总统计信息（新增）
echo -e "\n=== 步骤10/10: 汇总比对统计信息 ==="  # 新增步骤
echo "比对统计信息已保存至: 08.align_stats/"
echo "  - 原始比对统计: *raw.flagstat.txt"
echo "  - 筛选后比对统计: *unique.flagstat.txt"

echo -e "\n=== 分析流程完成 ==="
echo "结果文件:"
echo "count矩阵: 05.count_out/"
echo "TPM矩阵: 06.TPM_out/"
echo "FPKM矩阵: 07.FPKM_out/"
echo "比对统计: 08.align_stats/"  # 新增统计结果输出
echo "运行日志: $(date +'%Y%m%d_%H%M%S').log"
