# 参数参考

## 常用

| 参数 | 说明 | 默认 |
|---|---|---|
| `-G` | 参考基因组 FASTA | 必填（非 quant-only） |
| `-T` | FASTQ 目录 | 必填（非 quant-only） |
| `-t` | StringTie 用 GFF/GFF3 | 无 |
| `-a` | featureCounts 用 GTF | 无 |
| `-m` | `map` / `stringtie` / `featurecounts` / `both` / `quant-only` | `both` |
| `-p` | 线程数 | 12 |
| `-x` | 复用 HISAT2 索引前缀 | 无 |
| `--index-bam` | 对 BAM 建 CSI（小麦推荐） | 关 |
| `--force` | 强制重跑 | 关 |

## 输入模式

| 参数 | 说明 |
|---|---|
| `--fastp` | 输入 raw FASTQ，先跑 fastp |
| `--fastp-dir` | clean 输出目录（默认 `00.clean_fastq`） |
| `--fastp-report-dir` | fastp 报告（默认 `10.qc/fastp`） |
| `--fastp-threads` | fastp 线程（默认 8） |
| `--fastp-min-qual` | 合格 Phred（默认 20） |
| `--fastp-min-len` | 最短 read（默认 40） |
| `--fastp-compression` | gzip 级别（默认 4） |
| `--no-detect-adapter` | 关闭 PE 自动接头检测 |
| `--fastp-remove-raw` | QC 后删除 raw（默认保留） |

## 比对过滤

| 参数 | 说明 | 默认 |
|---|---|---|
| `--mapq-min` | 唯一比对最低 MAPQ | 10 |
| `--keep-multimap` | 保留多重比对 | 关（只留 NH:i:1） |
| `--hisat2-extra` | 追加 HISAT2 参数 | 无 |
| `--samtools-threads` | view/sort 线程 | `min(-p,4)` |

环境变量：`MAPQ_MIN`、`UNIQUE_ONLY`、`HISAT2_EXTRA`、`SAMTOOLS_THREADS`、`SAMTOOLS_SORT_MEM`。

## FASTQ 命名

| 参数 | 说明 | 默认 |
|---|---|---|
| `--r1-suffix` | R1 后缀 | `_1.fastq.gz` |
| `--r2-suffix` | R2 后缀 | `_2.fastq.gz` |
| `--r1-glob` | R1 glob | `*clean_1.fastq.gz` |

## featureCounts 透传

| 参数 | 说明 | 默认 |
|---|---|---|
| `-s` | 链特异 | 0 |
| `-Q` | 最低 MAPQ | 10 |
| `--fc-feature` | 特征类型 | exon |
| `--fc-attr` | 基因属性 | gene_id |
| `--skip-preflight` | 跳过 BAM/GTF 检查 | 关 |
| `--min-assigned` | Assigned% 告警阈值 | 50 |

## QC

| 参数 | 说明 |
|---|---|
| `--no-qc` | 不跑 10.qc |
| `--qc-only` | 只刷新 QC |
| `--title` | MultiQC 标题 |
| `--min-hisat2` | 比对率告警阈值 |

## 内存估算

- `samtools sort` 峰值 ≈ `SAMTOOLS_THREADS × SAMTOOLS_SORT_MEM`；
- 小麦 large index 常驻内存约 25–30 GB（`--mm` 时主要占 page cache）。
