#!/usr/bin/env python3
"""Merge per-sample quant dirs into gene x sample matrices (count/TPM/FPKM)."""
import sys, glob, os
import pandas as pd

root = sys.argv[1] if len(sys.argv) > 1 else "/public/home/xuyuchen/data/RNA-seq_public/01.public1/15.quant"
outdir = os.path.join(root, "merged"); os.makedirs(outdir, exist_ok=True)

SPECS = [("count", "gene_counts.tsv", "gene_counts.tsv"),
         ("TPM", "gene_TPM.tsv", "gene_TPM.tsv"),
         ("FPKM", "gene_FPKM.tsv", "gene_FPKM.tsv")]

for col, infile, outfile in SPECS:
    mats = {}
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(d): continue
        f = os.path.join(d, infile)
        if not os.path.exists(f): continue
        s = os.path.basename(d)
        df = pd.read_csv(f, sep="\t")
        mats[s] = df.set_index("gene_id")[col]
    if not mats:
        print(f"skip {col}: no samples")
        continue
    m = pd.DataFrame(mats).sort_index()
    m.index.name = "gene_id"
    if col == "count":
        m.to_csv(os.path.join(outdir, outfile), sep="\t", float_format="%d")
    else:
        m.to_csv(os.path.join(outdir, outfile), sep="\t", float_format="%.6f")
    print(col, m.shape)
