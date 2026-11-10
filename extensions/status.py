#!/usr/bin/env python3
"""Report per-sample quant status against the selection list."""
import sys, os, glob
import pandas as pd
P="/public/home/xuyuchen/data/RNA-seq_public/01.public1"
sel=pd.read_csv("/scratch/share/taoyongfu/xuyuchen/data/Wheat_RNA_Seq/wheat_T2T/rnaseq_quant_T2T/13.external_wheat/catalogue/final_selection_sizes.tsv",sep="\t",dtype=str)
sel=sel[sel.study!="CRA016389"].drop_duplicates("run")
root=sys.argv[1] if len(sys.argv)>1 else f"{P}/15.quant"
done=[]; other=[]
for _,r in sel.iterrows():
    p=f"{root}/{r.run}/gene_TPM.tsv"
    (done if os.path.exists(p) else other).append(r.run)
print(f"done={len(done)} pending={len(other)} total={len(sel)}")
if other: print("pending sample:", " ".join(other[:20]), "..." if len(other)>20 else "")
