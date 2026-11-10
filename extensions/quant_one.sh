#!/bin/bash
# 自定义单样本定量：raw FASTQ -> fastp -> HISAT2 -> featureCounts -> TPM/FPKM
# 参考 RNA_seq_pipeline 的参数口径，但输出更自由，可按需组合。
#
# 用法: quant_one.sh <study> <run>
# 环境变量（可选）:
#   THREADS=8 OUTROOT=... RAWROOT=... KEEP_BAM=1 KEEP_CLEAN=1
#   REF=... GTF=... SAMPLE_ID=...
set -euo pipefail

STUDY="$1"; RUN="$2"
P=/public/home/xuyuchen/data/RNA-seq_public/01.public1
REF="${REF:-/public/home/xuyuchen/01.scratch_home/data/Wheat_RNA_Seq/wheat_T2T/rnaseq_quant_T2T/01.ref/CS-IAAS_v1.1.fa}"
GTF="${GTF:-/public/home/xuyuchen/data/ref_genome/Wheat_CS_ref/3CS_IAAS_T2T/raw_data/CS-IAAS_v1.1_HC.gffread.gtf}"
RAWROOT="${RAWROOT:-$P/00.data/fastq}"
OUTROOT="${OUTROOT:-$P/15.quant}"
THREADS="${THREADS:-8}"
SAMPLE_ID="${SAMPLE_ID:-$RUN}"
KEEP_BAM="${KEEP_BAM:-0}"
KEEP_CLEAN="${KEEP_CLEAN:-0}"
MKFIFO="${HISAT2_MM:-0}"   # 默认不 mmap

FASTP=/public/home/xuyuchen/local/bin/fastp
ENV=/public/home/xuyuchen/miniconda3/envs/RNA-seq/bin
OUT="$OUTROOT/$SAMPLE_ID"
RAW="$RAWROOT/$STUDY"
mkdir -p "$OUT"
source /public/home/xuyuchen/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

r1="$RAW/${RUN}_1.fastq.gz"; [ -f "$r1" ] || r1="$RAW/${RUN}_r1.fq.gz"
r2="$RAW/${RUN}_2.fastq.gz"; [ -f "$r2" ] || r2="$RAW/${RUN}_r2.fq.gz"
c1="$OUT/${SAMPLE_ID}_clean_1.fastq.gz"; c2="$OUT/${SAMPLE_ID}_clean_2.fastq.gz"
cse="$OUT/${SAMPLE_ID}_clean.fastq.gz"
bam="$OUT/${SAMPLE_ID}.sorted.bam"
cnt="$OUT/${SAMPLE_ID}.counts.tsv"

[ -s "$OUT/gene_TPM.tsv" ] && { echo "[SKIP] $SAMPLE_ID"; exit 0; }

if [ -f "$r2" ]; then
  if [ ! -s "$c1" ] || [ ! -s "$c2" ]; then
    "$FASTP" --in1 "$r1" --in2 "$r2" --out1 "$c1" --out2 "$c2" \
      --unpaired1 "$OUT/unpaired.fq.gz" --unpaired2 "$OUT/unpaired.fq.gz" \
      --thread 8 --compression 4 --detect_adapter_for_pe \
      --qualified_quality_phred 20 --length_required 40 \
      --json "$OUT/fastp.json" --html "$OUT/fastp.html" --report_title "$SAMPLE_ID" \
      > "$OUT/fastp.log" 2>&1
  fi
  mmflag=""; [ "$MKFIFO" = "1" ] && mmflag="--mm"
  "$ENV/hisat2" -x "$REF" -p "$THREADS" --dta $mmflag --fr -1 "$c1" -2 "$c2" \
    2> "$OUT/hisat2.log" \
    | "$ENV/samtools" view -h -q 10 -@4 - \
    | awk 'BEGIN{FS=OFS="\t"} /^@/ || $0 ~ /NH:i:1($|\t)/' \
    | "$ENV/samtools" sort -m 1G -@4 -o "$bam" -
  FC_LAYOUT=(-p --countReadPairs)
else
  if [ ! -s "$cse" ]; then
    "$FASTP" --in1 "$r1" --out1 "$cse" --thread 8 --compression 4 \
      --qualified_quality_phred 20 --length_required 40 \
      --json "$OUT/fastp.json" --html "$OUT/fastp.html" --report_title "$SAMPLE_ID" \
      > "$OUT/fastp.log" 2>&1
  fi
  mmflag=""; [ "$MKFIFO" = "1" ] && mmflag="--mm"
  "$ENV/hisat2" -x "$REF" -p "$THREADS" --dta $mmflag -U "$cse" \
    2> "$OUT/hisat2.log" \
    | "$ENV/samtools" view -h -q 10 -@4 - \
    | awk 'BEGIN{FS=OFS="\t"} /^@/ || $0 ~ /NH:i:1($|\t)/' \
    | "$ENV/samtools" sort -m 1G -@4 -o "$bam" -
  FC_LAYOUT=()
fi

"$ENV/featureCounts" -T "$THREADS" -a "$GTF" -o "$cnt" -t exon -g gene_id \
  -Q 10 -s 0 "${FC_LAYOUT[@]}" --primary -C "$bam" > "$OUT/featureCounts.log" 2>&1

/public/home/xuyuchen/miniconda3/bin/python - "$cnt" "$OUT" <<'PY'
import sys, pandas as pd, numpy as np
cnt, out = sys.argv[1], sys.argv[2]
df = pd.read_csv(cnt, sep="\t", comment="#")
df = df[["Geneid","Length", df.columns[-1]]]
df.columns = ["gene_id","length","count"]
rate = df["count"]/df["length"]
df["TPM"] = rate/rate.sum()*1e6 if rate.sum()>0 else 0
lib = df["count"].sum()/1e6
df["FPKM"] = df["count"]/(df["length"]/1000*lib) if lib>0 else 0
df[["gene_id","length","count"]].to_csv(f"{out}/gene_counts.tsv", sep="\t", index=False)
df[["gene_id","TPM"]].to_csv(f"{out}/gene_TPM.tsv", sep="\t", index=False)
df[["gene_id","FPKM"]].to_csv(f"{out}/gene_FPKM.tsv", sep="\t", index=False)
print(f"genes={len(df)} counted={int((df['count']>0).sum())}")
PY

[ "$KEEP_BAM" -eq 0 ] && rm -f "$bam" "$bam.csi"
[ "$KEEP_CLEAN" -eq 0 ] && rm -f "$c1" "$c2" "$cse" "$OUT/unpaired.fq.gz"
echo "[DONE] $SAMPLE_ID"
