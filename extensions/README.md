# 项目级扩展脚本

这些脚本不修改 pipeline 本体，只调用其固定参数口径，用于：
- 下载完成即跑单样本
- 每样本独立目录
- 自动合并矩阵

## 文件

| 文件 | 作用 |
|---|---|
| `quant_one.sh` | 单样本：raw FASTQ → fastp → HISAT2（唯一过滤）→ featureCounts → TPM/FPKM |
| `merge_results.py` | 合并所有样本的 count/TPM/FPKM |
| `status.py` | 对照样本列表报告完成度 |
| `submit_quant_one.slurm` | SLURM 模板 |

## 用法

```bash
export THREADS=8
bash quant_one.sh <study> <run>

# 全部完成后
python3 merge_results.py /path/to/quant_root
python3 status.py /path/to/quant_root
```

## 可配置环境变量

`REF`、`GTF`、`RAWROOT`、`OUTROOT`、`THREADS`、`KEEP_BAM`、`KEEP_CLEAN`、`HISAT2_MM`、`SAMPLE_ID`。

默认参数与 `v5.6.x` pipeline 一致（fastp Q20/len40、MAPQ≥10、NH:i:1、featureCounts 同参数）。
