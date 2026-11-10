# 编排（SLURM）

Pipeline 只负责单样本/目录的比对与定量；调度、并发、批次属于编排层。

## 单作业（样本少）

```bash
sbatch --export=ALL orchestration/submit_map_amd.slurm
```

需要环境变量：`OUTDIR`、`GENOME`、`FASTQ_DIR`（可选 `HISAT2_INDEX`、`GTF`）。

## 一样本一 task（样本多，推荐）

```bash
export OUTDIR=/path/to/quant
export FASTQ_DIR=/path/to/clean_fastq
export GENOME=/path/to/genome.fa
export INDEX=/path/to/hisat2_prefix
export SAMPLE_LIST=/path/to/samples.txt   # 一行一个 sample 名
export THREADS=8
export SAMTOOLS_THREADS=4
export SAMTOOLS_SORT_MEM=1G

sbatch --array=0-N%8 --export=ALL orchestration/submit_map_array_amd.slurm
```

每个 task 调用：

```bash
map_only.sh -G genome.fa -T clean_fastq -p 8 -x index --only <sample>
```

- `--mm` 下同节点任务共享索引 page cache；
- 4–8 并发、每任务 8 线程通常最划算；
- 所有 task 写入同一个 `03.stringtie_sort/`。

## 汇总定量

映射全部完成后，**只跑一次**定量：

```bash
cd "$OUTDIR"
bash RNA_seq_run_dual.sh -m quant-only -t annotation.gff3 -a annotation.gtf -p 32
```

得到单张合并矩阵（与历史项目一致）。

## 依赖关系

```bash
JID=$(sbatch --parsable --array=0-N%8 --export=ALL orchestration/submit_map_array_amd.slurm)
sbatch --dependency=afterok:$JID --export=ALL orchestration/submit_featurecounts_amd.slurm
```

用 `afterok` 而不是 `afterany`：后者在数组部分 task 结束时就可能触发。

## 加样本

已提交的 array 不能追加 task。对新增样本再提交一个 array，输出到同一个 `OUTDIR`；`map_only.sh` 会跳过已有 BAM。
