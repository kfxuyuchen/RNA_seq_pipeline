# RNA-seq Pipeline（双路径）

统一入口：**比对一次 → StringTie 和/或 featureCounts**。

路径：`/public/home/xuyuchen/local/RNA_seq_pipeline/`

---

## 目录（唯一权威脚本）

| 文件 | 作用 |
|------|------|
| `RNA_seq_run_dual.sh` | **主入口**（推荐） |
| `map_only.sh` | HISAT2 比对（`--dta --mm`）→ sorted BAM |
| `post_stringtie.sh` | StringTie `-e` → count/TPM/FPKM 矩阵 |
| `post_featurecounts.sh` | featureCounts → count/FPKM/TPM |
| `counts_to_fpkm_tpm.py` | 由 featureCounts 矩阵计算 FPKM/TPM（Python 3；输出按 gene_id 排序） |
| `sort_expr_matrix.py` | 按首列 ID 对矩阵行排序（StringTie 路径在 post 末尾调用） |
| `prepDE.py` / `getTPM.py` / `getFPKM.py` | StringTie GTF 汇总矩阵（需系统 **python2**） |
| `RNA.yaml` | conda 环境（含 hisat2 / stringtie / subread / gffread） |
| `RNA-seq.nobuilds.yaml` | 跨机器宽松版环境文件 |
| `archive/` | 旧版 `RNA_seq_run1/2.sh`（已停用，勿再引用） |

> 项目目录 `wheat_T2T/scripts/` 里同名分析脚本仅为 **薄包装**，真正逻辑只在本目录维护。

---

## 流程示意

```
clean PE FASTQ
    │
    ▼
HISAT2 (--dta --mm) → unique → sort
    │
    ▼
03.stringtie_sort/*_sorted.bam
    │
    ├──────────────┬──────────────┐
    ▼              ▼              │
StringTie -e   featureCounts      │
    ▼              ▼              │
prepDE/TPM/FPKM  counts_to_fpkm_tpm
(python2)        (python3)
    ▼              ▼
04–07.*        09.featureCounts/
```

---

## 环境

```bash
conda activate RNA-seq
# 重建 / 迁移:
#   mamba env create -f RNA.yaml
#   # 或: mamba env create -n RNA-seq -f RNA-seq.nobuilds.yaml
```

| 工具 | 来源 |
|------|------|
| hisat2, stringtie, samtools, featureCounts, gffread, python3 | conda `RNA-seq` |
| python2（仅 StringTie 汇总脚本） | 系统 `/usr/bin/python2` |

自检：

```bash
hisat2 --version; stringtie --version; samtools --version
featureCounts -v; gffread --version; python3 --version; python2 --version
```

---

## 用法

在**输出目录**下运行（会写 `01.ref` / `03.stringtie_sort` / `04.*` / `09.*` 等）：

```bash
cd /path/to/rnaseq_quant_XXX
PIPE=/public/home/xuyuchen/local/RNA_seq_pipeline

# 两条定量都做（默认 -m both）
bash $PIPE/RNA_seq_run_dual.sh \
  -G genome.fa \
  -T /path/to/clean_for_quant \
  -t genes.gff3 \
  -a genes.gtf \
  -m both -p 12

# 只要 featureCounts
bash $PIPE/RNA_seq_run_dual.sh -G genome.fa -T fastq_dir -a genes.gtf -m featurecounts -p 12

# BAM 已齐，只做定量（给 -t 和/或 -a）
bash $PIPE/RNA_seq_run_dual.sh -m quant-only -t genes.gff3 -a genes.gtf -p 32
```

### `-m` 模式

| 模式 | 行为 |
|------|------|
| `map` | 只比对 |
| `stringtie` | 比对 + StringTie（需 `-t`） |
| `featurecounts` | 比对 + featureCounts（需 `-a`） |
| `both` | 比对 + 两条（需 `-t` 与 `-a`） |
| `quant-only` | 跳过比对；有 `-t` 跑 StringTie，有 `-a` 跑 featureCounts，两个都给则两条都跑 |

也可分步调用：

```bash
bash $PIPE/map_only.sh -G genome.fa -T fastq_dir -p 12
bash $PIPE/post_stringtie.sh -t genes.gff3 -p 16
bash $PIPE/post_featurecounts.sh -a genes.gtf -p 16
```

---

## 注释文件

| 用途 | 格式 |
|------|------|
| StringTie（`-t`） | GFF/GFF3 |
| featureCounts（`-a`） | GTF（官方或 `gffread -T` 转换的 `*.gffread.gtf`） |

GFF3 → GTF 示例：

```bash
gffread annotation.gff3 -T -o annotation.gffread.gtf
```

---

## 输出目录约定

| 路径 | 内容 |
|------|------|
| `01.ref/` | HISAT2 索引 |
| `03.stringtie_sort/` | `*_sorted.bam`（两条定量共用） |
| `04.stringtie_out/` | 每样本 StringTie GTF |
| `05.count_out/` `06.TPM_out/` `07.FPKM_out/` | StringTie 汇总矩阵 |
| `08.align_stats/` | flagstat |
| `09.featureCounts/` | `gene_count.tsv` / `gene_FPKM.tsv` / `gene_TPM.tsv` |

---

## 断点续传

所有阶段默认支持续跑，**同一命令重新提交即可**：

| 阶段 | 跳过条件 | 中断保护 |
|------|----------|----------|
| 比对 `map_only.sh` | 已有合法 `*_sorted.bam` | 写 `*.tmp` 再 `mv`；启动时清理残留 tmp；坏 BAM 删除重做 |
| StringTie | 已有非空 `*.assemble.gtf` | 输出先写 tmp；矩阵用 stamp，GTF 集合不变则跳过重算 |
| featureCounts | 最终矩阵存在且 BAM/注释 fingerprint 未变 | 输出先写 tmp；变了才重跑 |

强制全部重做：加 `--force`（会传给 map / post）。

```bash
bash RNA_seq_run_dual.sh ... --force
```


集群提交脚本在：

`/scratch/share/taoyongfu/xuyuchen/data/Wheat_RNA_Seq/wheat_T2T/scripts/`

其中 `submit_*.slurm` 负责资源申请；分析逻辑一律调用本 pipeline 目录，避免两套脚本漂移。

---

## 已废弃

`archive/RNA_seq_run1.sh`、`archive/RNA_seq_run2.sh` 为旧单路径全流程，**请改用 `RNA_seq_run_dual.sh`**。
