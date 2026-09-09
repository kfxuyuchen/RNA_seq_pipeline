# RNA-seq Pipeline（双路径）

统一入口：**比对一次 → StringTie 和/或 featureCounts → QC**。

路径：`/public/home/xuyuchen/local/RNA_seq_pipeline/`

版本：`RNA_seq_run_dual.sh` **v3.2.0**（FASTQ 命名可配 / FC 参数透传 / preflight / QC+MultiQC）

---

## 目录（唯一权威脚本）

| 文件 | 作用 |
|------|------|
| `RNA_seq_run_dual.sh` | **主入口**（推荐）；`bash … -h` 看完整帮助 |
| `map_only.sh` | HISAT2 比对（`--dta --mm`）→ sorted BAM |
| `post_stringtie.sh` | StringTie `-e` → count/TPM/FPKM 矩阵 |
| `post_featurecounts.sh` | featureCounts → count/FPKM/TPM（含 preflight） |
| `preflight_check.py` | BAM↔注释染色体 / gene 属性开跑前检查 |
| `counts_to_fpkm_tpm.py` | featureCounts → FPKM/TPM（按 gene_id 排序） |
| `sort_expr_matrix.py` | 按首列 ID 排序矩阵行 |
| `summarize_qc.py` / `run_qc.sh` | flagstat + featureCounts 汇总；可选 MultiQC |
| `pipeline_common.sh` | FASTQ 命名等公共函数 |
| `prepDE.py` / `getTPM.py` / `getFPKM.py` | StringTie GTF 汇总（系统 **python2**） |
| `RNA.yaml` / `RNA-seq.nobuilds.yaml` | conda 环境 |
| `archive/` | 旧版单路径脚本（已停用） |

> 项目目录里的同名 submit 脚本应为 **薄包装**，逻辑只在本目录维护。

---

## 流程示意

```
clean PE FASTQ（命名可配）
    │
    ▼
HISAT2 (--dta --mm) → unique → sort
    │
    ▼
03.stringtie_sort/*_sorted.bam
    │
    ├──────────────┬──────────────┐
    ▼              ▼              ▼
StringTie -e   featureCounts     QC
prepDE/TPM/FPKM  counts→FPKM/TPM  flagstat汇总
(python2)        (python3)        + MultiQC
    ▼              ▼              ▼
04–07.*        09.featureCounts/  10.qc/
```

---

## 环境

```bash
conda activate RNA-seq
# MultiQC（可选，用于 HTML 报告）:
#   conda install -c bioconda multiqc
#   # 或: pip install multiqc
```

| 工具 | 来源 |
|------|------|
| hisat2, stringtie, samtools, featureCounts, gffread, python3 | conda `RNA-seq` |
| multiqc（可选） | bioconda / pip |
| python2（仅 StringTie 汇总） | 系统 `/usr/bin/python2` |

---

## 用法

在**输出目录**下运行：

```bash
cd /path/to/rnaseq_quant_XXX
PIPE=/public/home/xuyuchen/local/RNA_seq_pipeline

# 完整帮助
bash $PIPE/RNA_seq_run_dual.sh -h

# 双路径（默认 *clean_1.fastq.gz）
bash $PIPE/RNA_seq_run_dual.sh \
  -G genome.fa -T /path/to/clean_fastq \
  -t genes.gff3 -a genes.gtf \
  -m both -p 32

# BAM 已齐：只跑 featureCounts（透传链特异等）
bash $PIPE/RNA_seq_run_dual.sh -m quant-only -a genes.gtf -s 2 -Q 10 -p 32

# 自定义 FASTQ 后缀
bash $PIPE/RNA_seq_run_dual.sh -G genome.fa -T fq_dir -a genes.gtf -m featurecounts \
  --r1-suffix _R1.fastq.gz --r2-suffix _R2.fastq.gz --r1-glob '*_R1.fastq.gz' -p 16

# 只刷新 QC / MultiQC
bash $PIPE/RNA_seq_run_dual.sh --qc-only --title "SRP585076 T2T"
```

### `-m` 模式

| 模式 | 行为 |
|------|------|
| `map` | 只比对 |
| `stringtie` | 比对 + StringTie（需 `-t`） |
| `featurecounts` | 比对 + featureCounts（需 `-a`） |
| `both` | 比对 + 两条（需 `-t` 与 `-a`） |
| `quant-only` | 跳过比对；按给定 `-t`/`-a` 定量 |

### FASTQ 命名

| 参数 | 默认 | 含义 |
|------|------|------|
| `--r1-glob` | `*clean_1.fastq.gz` | 在 `-T` 下查找 R1 |
| `--r1-suffix` | `_1.fastq.gz` | 从文件名剥出 sample id |
| `--r2-suffix` | `_2.fastq.gz` | 配对 R2 |

### featureCounts 透传（经 dual → post）

| 参数 | 默认 | 含义 |
|------|------|------|
| `-s / --strand` | `0` | 0/1/2 |
| `-Q / --mapq` | `10` | 最低 MAPQ |
| `--fc-feature` | `exon` | `-t` 特征 |
| `--fc-attr` | `gene_id` | `-g` 属性 |
| `--skip-preflight` | off | 跳过染色体/属性检查 |
| `--min-assigned` | `50` | Assigned% 告警阈值 |
| `--no-qc` | off | 结束后不跑 10.qc |

也可分步：

```bash
bash $PIPE/map_only.sh -G genome.fa -T fastq_dir -p 12
bash $PIPE/post_stringtie.sh -t genes.gff3 -p 16
bash $PIPE/post_featurecounts.sh -a genes.gtf -s 0 -p 16
bash $PIPE/run_qc.sh --title "my run"
```

---

## 注释文件

| 用途 | 格式 |
|------|------|
| StringTie（`-t`） | GFF/GFF3 |
| featureCounts（`-a`） | GTF（或 `gffread -T` 转换的 `*.gffread.gtf`） |

```bash
gffread annotation.gff3 -T -o annotation.gffread.gtf
```

`post_featurecounts` 默认会跑 `preflight_check.py`：检查特征/属性是否存在，以及 BAM 与注释染色体名是否有交集。

---

## 输出目录约定

| 路径 | 内容 |
|------|------|
| `01.ref/` | HISAT2 索引 |
| `03.stringtie_sort/` | `*_sorted.bam` |
| `04.stringtie_out/` | 每样本 StringTie GTF |
| `05.count_out/` `06.TPM_out/` `07.FPKM_out/` | StringTie 矩阵（行按 ID 排序） |
| `08.align_stats/` | flagstat |
| `09.featureCounts/` | count / FPKM / TPM |
| `10.qc/` | `qc_overview.tsv`、`*_summary.tsv`、`multiqc_report.html` |

---

## 断点续传

| 阶段 | 跳过条件 |
|------|----------|
| 比对 | 已有合法 `*_sorted.bam` |
| StringTie | 已有非空 `*.assemble.gtf`；矩阵 stamp 未变 |
| featureCounts | 矩阵存在且 BAM/注释/参数 fingerprint 未变 |

```bash
bash RNA_seq_run_dual.sh ... --force
```

---

## 已废弃

`archive/RNA_seq_run1.sh`、`archive/RNA_seq_run2.sh` → 请用 `RNA_seq_run_dual.sh`。
