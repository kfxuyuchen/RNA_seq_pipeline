# Changelog

## 5.6.1 — 2026-11-10

### Fixes
- `post_featurecounts.sh` 自动识别单端 BAM：SE 时不再传 `-p/--countReadPairs`。
  修复 SE 样本定量时报错 `No paired-end reads were detected in paired-end read library`。
- 混合 PE/SE BAM 的目录会明确报错，提示分开定量。

## 5.6.0 — 2026-11-10

### Highlights
- **单端（SE）支持**：`map_only.sh` 自动识别 SE clean FASTQ（默认 `*_clean.fastq.gz`），使用 `hisat2 -U`，沿用同一套 MAPQ/NH:i:1 过滤与排序。
- `fastp_qc.sh` 支持 SE：`<sample>.fastq.gz → <sample>_clean.fastq.gz`。
- `RNA_seq_run_dual.sh --se-glob` 透传 SE glob。
- PE 行为与命名完全不变；SE BAM 命名为 `<sample>_clean_sorted.bam`，与 PE 列名风格一致。

### Component versions
- `RNA_seq_run_dual.sh` **5.6.0**
- `map_only.sh` **2.7.0**
- `fastp_qc.sh` **1.2.0**

## 5.5.0 — 2026-11-10

### Highlights
- **可选关闭 `--mm`**：新增环境变量 `HISAT2_MM`（默认 1）与 `RNA_seq_run_dual.sh --no-mmap`。
  GPFS 上 `--mm` 会导致 hisat2 大量页错误、CPU 利用率仅约 1 核；关掉后索引常驻内存，速度明显提升。
- 默认行为不变（仍使用 `--mm`）。

### Component versions
- `RNA_seq_run_dual.sh` **5.5.0**
- `map_only.sh` **2.6.0**

## 5.4.2 — 2026-11-10

### Cleanup（减法）
- 删除 `archive/`（旧版脚本，历史保留在 git）。
- 删除 `templates/` 兼容软链；SLURM 模板统一使用 `orchestration/`。
- 删除 `__pycache__/`（已在 `.gitignore`）。
- README 同步移除兼容路径说明。
- 不改变任何脚本行为与默认参数。

## 5.4.1 — 2026-11-10

### Refactor（不改变任何默认行为）
- **编排移出核心**：SLURM 模板与提交脚本移至 `orchestration/`；`templates/` 保留兼容软链。
- **文档分层**：`docs/QUICKSTART.md`、`docs/PARAMETERS.md`、`docs/LEGACY.md`、`docs/ORCHESTRATION.md`；README 精简为总览 + 快速开始。
- **标注 legacy**：`getTPM.py` / `getFPKM.py` / `prepDE.py` / `post_stringtie.sh` 头部注明仅用于可选 StringTie 路线。
- 未修改任何脚本参数、默认值或输出格式。

## 5.4.0 — 2026-11-10

### Highlights
- **按样本子集并行映射**：`map_only.sh` 新增 `--only S1,S2` 与 `--sample-list FILE`。
- 新增 `templates/submit_map_array_amd.slurm` 与 `submit_map_array.sh`：
  一个 SLURM array task 只映射一个样本，多任务共享 `--mm` page cache。
- 适合“下载分批、先到先跑”的场景：映射并行，最后一次性 `post_featurecounts` 汇总。

### Component versions
- `map_only.sh` **2.5.0**
- `RNA_seq_run_dual.sh` **5.3.0**

## 5.3.0 — 2026-11-10

### Highlights
- **samtools 线程可外部配置**：新增环境变量 `SAMTOOLS_THREADS` 与 `RNA_seq_run_dual.sh --samtools-threads N`。
- 默认仍为 `min(-p, 4)`，保持历史行为与内存占用；只有显式传参才放大。
- 文档补充内存估算：`samtools sort` 内存 ≈ `SAMTOOLS_THREADS × SAMTOOLS_SORT_MEM`。

### Component versions
- `RNA_seq_run_dual.sh` **5.3.0**
- `map_only.sh` **2.4.3**

## 5.2.0 — 2026-11-10

### Highlights
- **比对过滤参数外部可配**：`map_only.sh` 新增环境变量 `MAPQ_MIN`（默认 10）、`UNIQUE_ONLY`（默认 1）、`HISAT2_EXTRA`。
- `RNA_seq_run_dual.sh` 新增同名命令行入口：`--mapq-min N`、`--keep-multimap`、`--hisat2-extra "ARGS"`。
- 默认值保持不变，与 252 样本项目口径完全一致，仅提供可选覆盖。

### Component versions
- `RNA_seq_run_dual.sh` **5.2.0**
- `map_only.sh` **2.4.2**
- `fastp_qc.sh` **1.0.0**

## 5.1.0 — 2026-11-10

### Highlights
- **双输入模式显式化**：`RNA_seq_run_dual.sh` 启动时打印 `INPUT mode: clean` 或 `INPUT mode: raw → fastp`。
- **输入自检**：未加 `--fastp` 却只找到 raw FASTQ，或加了 `--fastp` 却只找到 clean FASTQ 时给出告警，避免拿错目录。
- 兼容两套工作习惯：raw reads 入口（`--fastp`）与 clean reads 入口（默认）共用同一套比对/定量逻辑。

### Component versions
- `RNA_seq_run_dual.sh` **5.1.0**
- `fastp_qc.sh` **1.0.0**
- `map_only.sh` **2.4.1**

## 5.0.1 — 2026-11-10

### Fixes
- `map_only.sh` 的 `-x/--index` 不再对索引前缀做 symlink 解析：改用 `realpath -s`。
  修复当参考 FASTA 是指向真实文件的软链接、索引放在软链接旁边时，`.1.ht2(l)` 找不到的问题。
- `map_only.sh` **2.4.1**。

## 5.0.0 — 2026-11-10

### Highlights
- **内置 fastp 步骤**：新增 `fastp_qc.sh` v1.0.0，`RNA_seq_run_dual.sh --fastp` 可直接以 raw FASTQ 为输入，先 QC 再比对。
- 固定实验室参数：`--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 40 --compression 4`，与 252 样本项目一致。
- 输出命名 `<sample>_clean_1/2.fastq.gz`，可被 `map_only.sh` 默认 `*clean_1.fastq.gz` 直接消费。
- fastp 报告写入 `10.qc/fastp/`（json/html/log/unpaired），默认**不删除** raw FASTQ（`--fastp-remove-raw` 可显式开启）。
- 兼容原有 clean-input 用法：不加 `--fastp` 时行为完全不变。

### Component versions
- `RNA_seq_run_dual.sh` **5.0.0**
- `fastp_qc.sh` **1.0.0**
- `map_only.sh` **2.4.0**

## 4.0.0 — 2026-10-29

### Highlights
- **Default one-write mapping**: `hisat2 | unique(MAPQ≥10, NH:i:1) | sort` → only `03.stringtie_sort/*_sorted.bam`.
- **Continuous dirs kept**: `01.ref` → `02.hisat2-mapping` (logs) → `03` … → `10.qc`.
- **HISAT2 overall rate persisted**: `08.align_stats/*_hisat2_rate.tsv` + `hisat2_alignment_rates.tsv`; QC overview includes `hisat2_overall_pct`.
- **Optional CSI BAM index** (`--index-bam`): always `samtools index -c` (wheat-safe; BAI fails when chr > 2^29). Index failure does not invalidate BAM.
- **Reuse external HISAT2 index**: `-x` / `HISAT2_INDEX`; never delete existing indexes by default.
- **Official smoke**: `bash smoke_test.sh` (≤5 min toy data).
- **Slurm templates**: `templates/submit_map_amd.slurm`, `templates/submit_featurecounts_amd.slurm` (no `#SBATCH --mem`).

### Component versions
- `RNA_seq_run_dual.sh` **4.0.0**
- `map_only.sh` **2.4.0**
- `run_qc.sh` / `summarize_qc.py` **1.1.0**

## 3.2.0 — 2026-09-10

- Configurable FASTQ naming, featureCounts passthrough, QC + MultiQC.
