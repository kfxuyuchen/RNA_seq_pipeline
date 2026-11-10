# RNA-seq Pipeline

**Version 5.4.1** · 一次比对，StringTie / featureCounts 双路径定量，带断点续传、QC 与并行映射。

权威目录：`/public/home/xuyuchen/local/RNA_seq_pipeline/`  
仓库：<https://github.com/kfxuyuchen/RNA_seq_pipeline>

---

## 设计原则

| 原则 | 做法 |
|---|---|
| 单一职责 | 只负责 `FASTQ → BAM → 表达矩阵`；下载、批次编排、矩阵合并放在外部 |
| 默认稳定 | 默认参数与历史项目一致；新增能力不改变旧行为 |
| 可复现 | 版本、参数、日志、QC 全部落盘 |
| 可续跑 | 已完成步骤自动跳过，`--force` 重跑 |
| 先验后用 | `smoke_test.sh` 冒烟通过再全量 |
| 少写盘 | 默认一次落盘到 `03.stringtie_sort/*_sorted.bam` |

## 快速开始

```bash
conda activate RNA-seq

# clean reads（默认）
bash RNA_seq_run_dual.sh -G genome.fa -T clean_fastq \
  -t annotation.gff3 -a annotation.gtf -m both -p 16 -x <hisat2_prefix> --index-bam

# raw reads（内置 fastp）
bash RNA_seq_run_dual.sh -G genome.fa -T raw_fastq --fastp \
  -a annotation.gtf -m featurecounts -p 16 -x <hisat2_prefix>

# 已有 BAM，只补定量
cd /path/to/rnaseq_quant && bash RNA_seq_run_dual.sh -m quant-only -a annotation.gtf -p 32
```

详细用法：[`docs/QUICKSTART.md`](docs/QUICKSTART.md)

## 输出目录

```
01.ref/  02.hisat2-mapping/  03.stringtie_sort/
04.stringtie_out/  05.count_out/  06.TPM_out/  07.FPKM_out/
08.align_stats/  09.featureCounts/  10.qc/
```

主分析矩阵用 `09.featureCounts/`；StringTie 结果为可选/历史兼容，见 [`docs/LEGACY.md`](docs/LEGACY.md)。

## 参数与覆盖

常用参数见 [`docs/PARAMETERS.md`](docs/PARAMETERS.md)。默认值即合同，按需从外部覆盖：

```bash
--mapq-min 20            # 比对最低 MAPQ
--keep-multimap          # 保留多重比对
--hisat2-extra "..."     # 追加 HISAT2 参数
--samtools-threads 8     # samtools 线程
```

## 并行映射

一样本一 task，适合大样本量：

```bash
sbatch --array=0-N%8 --export=ALL orchestration/submit_map_array_amd.slurm
```

见 [`docs/ORCHESTRATION.md`](docs/ORCHESTRATION.md)。SLURM 模板位于 `orchestration/`，`templates/` 保留兼容软链。

## 验收

```bash
bash smoke_test.sh          # 期望 SMOKE_PASS
```

## 组件

| 文件 | 作用 |
|---|---|
| `RNA_seq_run_dual.sh` | 主入口 |
| `map_only.sh` | 比对（唯一过滤 + 排序） |
| `fastp_qc.sh` | raw → clean |
| `post_stringtie.sh` / `post_featurecounts.sh` | 定量 |
| `run_qc.sh` / `summarize_qc.py` | QC 汇总 |
| `counts_to_fpkm_tpm.py` | featureCounts 路线的 TPM/FPKM |
| `orchestration/` | SLURM 模板与提交脚本 |
| `docs/` | 快速开始、参数、legacy、编排 |
