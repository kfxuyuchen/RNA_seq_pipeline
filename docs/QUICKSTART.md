# 快速开始

## 1. 环境

```bash
conda activate RNA-seq
# 需要: hisat2, samtools, stringtie, featureCounts, gffread, python3
# 可选: fastp（PATH 中或 FASTP=/path/to/fastp）、multiqc
```

## 2. 从 clean reads 开始（默认）

```bash
bash RNA_seq_run_dual.sh \
  -G genome.fa -T clean_fastq \
  -t annotation.gff3 -a annotation.gtf \
  -m both -p 16 -x /path/to/hisat2_prefix --index-bam
```

- `clean_fastq/` 内为 `<sample>_clean_1.fastq.gz` / `<sample>_clean_2.fastq.gz`；
- `-m both` = StringTie + featureCounts；只要基因层面可用 `-m featurecounts`。

## 3. 从 raw reads 开始

```bash
bash RNA_seq_run_dual.sh \
  -G genome.fa -T raw_fastq --fastp \
  -a annotation.gtf -m featurecounts \
  -p 16 -x /path/to/hisat2_prefix
```

fastp 参数固定为实验室标准，输出到 `00.clean_fastq/`。

## 4. 已有 BAM，只补定量

```bash
cd /path/to/rnaseq_quant
bash RNA_seq_run_dual.sh -m quant-only -t annotation.gff3 -a annotation.gtf -p 32
```

## 5. 输出

```
01.ref/                 HISAT2 索引（或软链）
02.hisat2-mapping/      *.hisat2.log
03.stringtie_sort/      *_sorted.bam（MAPQ≥10 + NH:i:1）
04.stringtie_out/       StringTie GTF
05.count_out/ 06.TPM_out/ 07.FPKM_out/   StringTie 矩阵
08.align_stats/         flagstat / 比对率
09.featureCounts/       gene_counts.txt / gene_TPM.tsv / gene_FPKM.tsv
10.qc/                  QC 汇总 / MultiQC
```

## 6. 续跑与冒烟

- 同一条命令重新提交即可跳过已完成步骤；
- `--force` 强制重跑；
- `bash smoke_test.sh` 应输出 `SMOKE_PASS`。
