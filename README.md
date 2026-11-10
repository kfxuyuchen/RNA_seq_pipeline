# RNA-seq Pipeline

**Version 4.0.0** · 一次比对，StringTie / featureCounts 双路径定量，带断点续传与 QC。

权威脚本目录：`/public/home/xuyuchen/local/RNA_seq_pipeline/`  
仓库：<https://github.com/kfxuyuchen/RNA_seq_pipeline>

---

## 设计原则

| 原则 | 做法 |
|------|------|
| 目录可读 | 固定 `01`–`10` 连续编号，跨项目一致 |
| 默认省 I/O | 比对管道 **一次落盘** 到 sorted BAM；不写 raw/unique 中间文件 |
| 先验再用 | 已有 HISAT2 索引用 `-x` / `HISAT2_INDEX`；禁止无故删除旧索引 |
| 小麦友好 | 可选 `--index-bam` → **CSI**（勿用默认 BAI；超长染色体会失败） |
| 可验收 | `smoke_test.sh` 小数据冒烟；QC 表含 hisat2 overall% 与 Assigned% |
| 集群习惯 | Slurm 模板 **不写 `--mem`**，用分区核存比换算 `-c` |

---

## 流程一览

```
clean PE FASTQ
      │
      ▼
 HISAT2 (--dta --mm) | unique | sort     ← 默认只写一份 BAM
      │                    │
      ▼                    ▼
 02.hisat2-mapping/   03.stringtie_sort/*_sorted.bam
   *.hisat2.log         (+ optional *.csi)
      │
      ├──────────────┬──────────────┐
      ▼              ▼              ▼
 StringTie -e   featureCounts     QC
 04–07.*        09.featureCounts/  10.qc/
```

`--debug` 时恢复旧式 `02` 中间 BAM（raw → unique → sort），便于排错。

---


## 并行映射（v5.4.0+）

样本多时不要顺序跑，用 SLURM array：一个 task 一个样本。

```bash
# 1) 样本名单（sample 名，不含 _clean 后缀）
cat > samples.txt <<LIST
sampleA
sampleB
sampleC
LIST

# 2) 提交 array，最多 8 个并发
export OUTDIR=/path/to/rnaseq_quant
export FASTQ_DIR=/path/to/clean_fastq
export GENOME=/path/to/genome.fa
export INDEX=/path/to/hisat2_prefix
export SAMPLE_LIST=/path/to/samples.txt
export THREADS=8
export SAMTOOLS_THREADS=4
export SAMTOOLS_SORT_MEM=1G
sbatch --array=0-2%8 --export=ALL templates/submit_map_array_amd.slurm

# 3) 全部映射完成后，一次性汇总
cd "$OUTDIR"
bash RNA_seq_run_dual.sh -m quant-only -t annotation.gff3 -a annotation.gtf -p 32
```

`--mm` 下多个 task 共享同一份索引 page cache，内存不会随并发数线性增长。

## 参数覆盖（v5.2.0+）

默认口径与 252 样本项目一致；需要改变时不必改脚本，直接外部传参：

```bash
# 保留多重比对
bash RNA_seq_run_dual.sh ... --keep-multimap

# 改最低 MAPQ
bash RNA_seq_run_dual.sh ... --mapq-min 20

# 给 hisat2 追加参数
bash RNA_seq_run_dual.sh ... --hisat2-extra "--max-intronlen 5000"
```

`samtools sort` 的内存估算：**线程数 × `-m`**。例如 `--samtools-threads 8` + `SAMTOOLS_SORT_MEM=2G` 约 16 GB。
默认仍是 `min(-p, 4)` × `1G`，保持历史内存行为。

也可直接设环境变量调用 `map_only.sh`：

```bash
MAPQ_MIN=20 UNIQUE_ONLY=0 HISAT2_EXTRA="--max-intronlen 5000" \
  bash map_only.sh -G genome.fa -T clean_fastq -p 24
```

其它常用覆盖：

| 需求 | 入口 |
|---|---|
| 线程 | `-p N` |
| 链特异 | `-s 0/1/2` |
| featureCounts MAPQ | `-Q N` |
| 特征类型/属性 | `--fc-feature` / `--fc-attr` |
| samtools view/sort 线程 | `--samtools-threads N` 或环境变量 `SAMTOOLS_THREADS` |
| samtools sort 内存 | 环境变量 `SAMTOOLS_SORT_MEM` |
| fastp 参数 | `--fastp-min-qual` / `--fastp-min-len` / `--fastp-compression` / `--no-detect-adapter` |

## 两种输入模式（v5.1.0+）

| 模式 | 输入 | 命令 | 说明 |
|---|---|---|---|
| clean | `*clean_1.fastq.gz` | 不加 `--fastp` | 与旧版行为一致；适合已经 QC 过的数据 |
| raw | `*_1.fastq.gz` | 加 `--fastp` | 先 fastp 再比对；同一条命令完成 QC+定量 |

启动时会打印当前模式，例如：

```
INPUT mode: raw FASTQ → fastp → clean（找到 6 个 raw R1）
```

如果目录内容与模式不匹配，会给出 `WARN`：

- 没加 `--fastp` 但只有 raw 文件 → 提示加 `--fastp`；
- 加了 `--fastp` 但只有 clean 文件 → 提示会二次 trim。

## 内置 fastp（v5.0.0+）

如果手上只有 raw FASTQ，可以跳过外置 QC，直接用 `--fastp`：

```bash
bash RNA_seq_run_dual.sh \
  -G genome.fa -T raw_fastq --fastp \
  -a genes.gtf -m featurecounts \
  -p 32 -x /path/to/hisat2_prefix
```

- fastp 固定参数与实验室 252 样本项目一致：
  `--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 40 --compression 4`
- clean 输出：`00.clean_fastq/<sample>_clean_1/2.fastq.gz`
- 报告：`10.qc/fastp/<sample>.fastp.{json,html,log}`
- 默认保留 raw FASTQ；需要删除时加 `--fastp-remove-raw`
- 不加 `--fastp` 时，行为与旧版完全一致（输入为 clean FASTQ）

单独运行：

```bash
bash fastp_qc.sh -T raw_fastq -o 00.clean_fastq -j 10.qc/fastp -p 8
```

依赖：`fastp`（PATH 中可用，或设置 `FASTP=/path/to/fastp`）。

## 输出目录（01–10 连续）

| 目录 | 内容 |
|------|------|
| `01.ref/` | HISAT2 索引，或指向已有索引的软链 |
| `02.hisat2-mapping/` | `*.hisat2.log`；仅 `--debug` 时保留 `*_raw.bam` / `*_unique.bam` |
| `03.stringtie_sort/` | `*_sorted.bam`（定量共用）；`--index-bam` 时另有 `*.csi` |
| `04.stringtie_out/` | 每样本 StringTie GTF |
| `05.count_out/` `06.TPM_out/` `07.FPKM_out/` | StringTie 矩阵（行按 ID 排序） |
| `08.align_stats/` | `*_unique.flagstat.txt`、`*_hisat2_rate.tsv`、`hisat2_alignment_rates.tsv` |
| `09.featureCounts/` | `gene_counts.txt`、count / FPKM / TPM |
| `10.qc/` | `qc_overview.tsv`、`hisat2_alignment_rates.tsv`、MultiQC（可选） |

---

## 快速开始

```bash
conda activate RNA-seq
PIPE=/public/home/xuyuchen/local/RNA_seq_pipeline
cd /path/to/rnaseq_quant_XXX

# 官方冒烟（玩具数据，应打印 SMOKE_PASS）
bash $PIPE/smoke_test.sh

# 双路径（默认 FASTQ：*clean_1.fastq.gz / *clean_2.fastq.gz）
bash $PIPE/RNA_seq_run_dual.sh \
  -G genome.fa -T /path/to/clean_fastq \
  -t genes.gff3 -a genes.gtf \
  -m both -p 16 \
  -x /path/to/existing_hisat2_prefix \
  --index-bam \
  --min-hisat2 60

# BAM 已齐：只跑 featureCounts
bash $PIPE/RNA_seq_run_dual.sh -m quant-only -a genes.gtf -s 0 -Q 10 -p 16

# 只刷新 QC
bash $PIPE/RNA_seq_run_dual.sh --qc-only --title "my project"
```

### 模式 `-m`

| 模式 | 行为 |
|------|------|
| `map` | 只比对 |
| `stringtie` | 比对 + StringTie（需 `-t`） |
| `featurecounts` | 比对 + featureCounts（需 `-a`） |
| `both` | 比对 + 两条定量（默认） |
| `quant-only` | 跳过比对；按 `-t` / `-a` 定量 |

### 常用参数

| 参数 | 含义 |
|------|------|
| `-x` / `HISAT2_INDEX` | 复用已有 HISAT2 索引前缀 |
| `--index-bam` | 对 sorted BAM 建 **CSI** 索引 |
| `-d` | 多步落盘并保留 `02` 中间 BAM |
| `--r1-glob` / `--r1-suffix` / `--r2-suffix` | FASTQ 命名 |
| `-s` `-Q` `--fc-feature` `--fc-attr` | featureCounts 透传 |
| `--min-assigned` / `--min-hisat2` | QC 告警阈值 |
| `--force` | 强制重跑 |

分步调用：

```bash
bash $PIPE/map_only.sh -G genome.fa -T fastq_dir -p 12 -x INDEX_PREFIX --index-bam
bash $PIPE/post_stringtie.sh -t genes.gff3 -p 16
bash $PIPE/post_featurecounts.sh -a genes.gtf -s 0 -Q 10 -p 16
bash $PIPE/run_qc.sh --title "SRP… T2T" --min-hisat2 60
```

---

## Slurm（AMD 分区模板）

模板在 `templates/`，**禁止 `#SBATCH --mem`**，用 `-c` 换内存（`AMD_9A14` ≈ 3600 MB/CPU）。

```bash
export OUTDIR=/path/to/rnaseq_quant_XXX
export GENOME=/path/to/genome.fa
export FASTQ_DIR=/path/to/clean_fastq
export HISAT2_INDEX=/path/to/CS-IAAS_v1.1.fa   # 可选，强烈推荐复用
export GTF=/path/to/genes.gffread.gtf

mkdir -p "$OUTDIR/logs"
sbatch --export=ALL $PIPE/templates/submit_map_amd.slurm
# 比对完成后：
sbatch --export=ALL $PIPE/templates/submit_featurecounts_amd.slurm
```

---

## 环境

```bash
conda activate RNA-seq
# 可选 MultiQC：
#   conda install -c bioconda multiqc
```

| 工具 | 用途 |
|------|------|
| hisat2, samtools, stringtie, featureCounts, gffread, python3 | 主流程 |
| multiqc | `10.qc` HTML（可选） |
| python2 | 仅 StringTie 汇总（`prepDE` / getTPM / getFPKM） |

环境文件：`RNA.yaml`、`RNA-seq.nobuilds.yaml`。

---

## 注释

| 用途 | 格式 |
|------|------|
| StringTie（`-t`） | GFF/GFF3 |
| featureCounts（`-a`） | GTF（可用 `gffread -T` 转换） |

```bash
gffread annotation.gff3 -T -o annotation.gffread.gtf
```

`post_featurecounts` 默认跑 `preflight_check.py`（特征/属性、BAM↔注释染色体交集）。

featureCounts 默认与实验室 252 样本 T2T 基线一致：

`-t exon -g gene_id -Q 10 -s 0 -p --countReadPairs --primary -C`

---

## 断点续传与验收

| 阶段 | 跳过条件 |
|------|----------|
| 比对 | 已有合法 `03.stringtie_sort/*_sorted.bam`（可补写比对率 / CSI） |
| StringTie | 已有非空 assemble GTF；矩阵 stamp 未变 |
| featureCounts | 矩阵存在且 BAM/注释/参数 fingerprint 未变 |

强制重跑：加 `--force`。

冒烟：

```bash
bash smoke_test.sh
# 期望：SMOKE_PASS
```

---

## 脚本清单

| 文件 | 作用 |
|------|------|
| `RNA_seq_run_dual.sh` | **主入口** v4.0.0 |
| `map_only.sh` | 比对（默认一次落盘） |
| `post_stringtie.sh` / `post_featurecounts.sh` | 定量 |
| `run_qc.sh` / `summarize_qc.py` | QC 表 + MultiQC |
| `preflight_check.py` | FC 开跑前检查 |
| `smoke_test.sh` | 官方小数据冒烟 |
| `templates/*.slurm` | AMD 分区提交模板 |
| `pipeline_common.sh` | FASTQ 命名公共函数 |
| `archive/` | 旧版单路径脚本（停用） |

变更记录见 [`CHANGELOG.md`](CHANGELOG.md)。版本号见 [`VERSION`](VERSION)。

---

## 常见问题

**Q: 为什么默认不写 `02` 里的 BAM？**  
A: 定量只需要 sorted BAM；多写 raw/unique 再删会浪费 I/O。排错时用 `-d`。

**Q: 小麦 `samtools index` 报 BAI / out of range？**  
A: 用 `--index-bam`（CSI）。featureCounts **不需要** BAM 索引即可计数。

**Q: 索引要不要重建？**  
A: 先复用。`-x` 指向已有前缀；`01.ref` 已有 `.ht2`/`.ht2l` 会自动跳过 build。不要擅自 `rm` 旧索引。

**Q: overall 比对率在哪？**  
A: `02/*.hisat2.log`，以及 `08.align_stats/hisat2_alignment_rates.tsv` → `10.qc/qc_overview.tsv` 的 `hisat2_overall_pct`。
