# Legacy：StringTie 路线

## 状态

StringTie 路线（`04.stringtie_out/`、`05.count_out/`、`06.TPM_out/`、`07.FPKM_out/`）保留用于：

- 转录本/isoform 层面分析；
- 与历史项目（252 样本）对齐；
- 交叉验证 featureCounts 结果。

新项目的**主矩阵建议只用 featureCounts（`09.featureCounts/`）**，原因：

- featureCounts 基于已知注释、binomial 模型稳定；
- StringTie 依赖组装与注释，跨批次更容易漂移；
- StringTie 路线的 `getTPM.py` / `getFPKM.py` / `prepDE.py` 依赖 python2。

## 三个 python2 脚本

| 脚本 | 作用 |
|---|---|
| `prepDE.py` | 从 StringTie GTF 汇总 count 矩阵 |
| `getTPM.py` | 计算 TPM |
| `getFPKM.py` | 计算 FPKM |

由 `post_stringtie.sh` 调用，顺序固定。环境需要 `/usr/bin/python2`。

## 如何只跑 featureCounts

```bash
bash RNA_seq_run_dual.sh -m featurecounts -a annotation.gtf -p 32
```

或在已有 BAM 上：

```bash
bash RNA_seq_run_dual.sh -m quant-only -a annotation.gtf -p 32
```

这样不会调用 StringTie 和三个 python2 脚本。

## 何时才需要 both

```bash
bash RNA_seq_run_dual.sh -m both -t annotation.gff3 -a annotation.gtf -p 32
```

仅在需要 isoform/转录本结果，或要与旧项目完全对齐时使用。
