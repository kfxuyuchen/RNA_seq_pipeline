#!/usr/bin/env python3
"""Compute gene-level FPKM and TPM from a featureCounts matrix.

Input: featureCounts output (Geneid, ..., Length, sample columns).
Output: gene_count.tsv, gene_FPKM.tsv, gene_TPM.tsv (and optional CSV copies).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_featurecounts(path: Path):
    genes = []
    lengths = []
    samples = []
    counts = []  # list of lists, rows=genes
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if row[0] == "Geneid":
                # Geneid Chr Start End Strand Length <samples...>
                if "Length" not in row:
                    raise SystemExit(f"ERROR: no Length column in {path}")
                li = row.index("Length")
                samples = row[li + 1 :]
                # strip path/bam suffixes from sample headers for readability
                samples = [Path(s).name.replace("_sorted.bam", "").replace(".bam", "") for s in samples]
                continue
            genes.append(row[0])
            lengths.append(int(float(row[li])))
            counts.append([int(float(x)) for x in row[li + 1 :]])
    if not samples:
        raise SystemExit(f"ERROR: failed to parse header in {path}")
    return genes, lengths, samples, counts


def write_matrix(path: Path, genes, samples, matrix, fmt="%.6f"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id"] + samples)
        for g, row in zip(genes, matrix):
            if fmt == "int":
                w.writerow([g] + [str(int(x)) for x in row])
            else:
                w.writerow([g] + [fmt % x for x in row])


def to_fpkm_tpm(counts, lengths):
    """counts: list[list[int]] genes x samples; lengths: list[int] bp."""
    n_genes = len(counts)
    n_samples = len(counts[0]) if n_genes else 0
    fpkm = [[0.0] * n_samples for _ in range(n_genes)]
    tpm = [[0.0] * n_samples for _ in range(n_genes)]

    for j in range(n_samples):
        lib = sum(counts[i][j] for i in range(n_genes))
        rpk = [0.0] * n_genes
        for i in range(n_genes):
            L = lengths[i]
            c = counts[i][j]
            if L > 0:
                rpk[i] = c / (L / 1000.0)
                if lib > 0:
                    fpkm[i][j] = c / ((L / 1000.0) * (lib / 1e6))
            else:
                rpk[i] = 0.0
                fpkm[i][j] = 0.0
        rpk_sum = sum(rpk)
        if rpk_sum > 0:
            for i in range(n_genes):
                tpm[i][j] = rpk[i] / rpk_sum * 1e6
    return fpkm, tpm


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True, help="featureCounts counts file")
    ap.add_argument("-o", "--outdir", required=True, help="output directory")
    ap.add_argument("--prefix", default="gene", help="output file prefix (default: gene)")
    args = ap.parse_args()

    inp = Path(args.input)
    out = Path(args.outdir)
    genes, lengths, samples, counts = read_featurecounts(inp)
    fpkm, tpm = to_fpkm_tpm(counts, lengths)

    write_matrix(out / f"{args.prefix}_count.tsv", genes, samples, counts, fmt="int")
    write_matrix(out / f"{args.prefix}_FPKM.tsv", genes, samples, fpkm)
    write_matrix(out / f"{args.prefix}_TPM.tsv", genes, samples, tpm)

    # also write length table for reference
    with (out / f"{args.prefix}_length.tsv").open("w") as f:
        f.write("gene_id\tlength_bp\n")
        for g, L in zip(genes, lengths):
            f.write(f"{g}\t{L}\n")

    print(f"genes={len(genes)} samples={len(samples)}")
    print(f"wrote {out / (args.prefix + '_count.tsv')}")
    print(f"wrote {out / (args.prefix + '_FPKM.tsv')}")
    print(f"wrote {out / (args.prefix + '_TPM.tsv')}")
    print(f"wrote {out / (args.prefix + '_length.tsv')}")


if __name__ == "__main__":
    main()
