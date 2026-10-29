# Changelog

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
