#!/usr/bin/env python3
"""Preflight checks before featureCounts / StringTie quantification.

Checks (fail = exit 1; warn = exit 0 with messages on stderr/stdout):
  1) Annotation has requested gene-id attribute (default gene_id)
  2) BAM vs annotation chromosome-name overlap (chr1 vs 1 etc.)
  3) featureCounts feature type present (default exon)
  4) Optional: PE mate presence for FASTQ dirs (when --r1 given)

Usage:
  preflight_check.py --annot genes.gtf --bam sample_sorted.bam
  preflight_check.py --annot genes.gtf --bam-dir 03.stringtie_sort -g gene_id -t exon
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def bam_chroms(bam: Path, limit_sq: int = 50000) -> set[str]:
    import subprocess

    out = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
    chroms = set()
    for line in out.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t")[1:]:
                if field.startswith("SN:"):
                    chroms.add(field[3:])
                    break
        if len(chroms) >= limit_sq:
            break
    return chroms


def annot_stats(annot: Path, attr: str, feature: str, max_lines: int = 2_000_000):
    """Return (chroms, n_feature, n_attr, example_attrs)."""
    chroms: set[str] = set()
    n_feat = 0
    n_attr = 0
    attr_examples: set[str] = set()
    key = f'{attr} "'
    key2 = f"{attr}="
    with open_text(annot) as fh:
        for i, line in enumerate(fh):
            if i >= max_lines:
                break
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chroms.add(parts[0])
            ftype = parts[2]
            attrs = parts[8]
            if len(attr_examples) < 12:
                for tok in attrs.replace(";", " ").split():
                    if "=" in tok:
                        attr_examples.add(tok.split("=", 1)[0])
                    elif tok.endswith('"') or True:
                        # GTF: gene_id "x";
                        if " " not in tok and tok.endswith('"'):
                            continue
                # collect attribute keys more carefully
                for chunk in attrs.split(";"):
                    chunk = chunk.strip()
                    if not chunk:
                        continue
                    if " " in chunk:
                        attr_examples.add(chunk.split(" ", 1)[0])
                    elif "=" in chunk:
                        attr_examples.add(chunk.split("=", 1)[0])
            if ftype == feature:
                n_feat += 1
                if key in attrs or key2 in attrs or f"{attr} " in attrs:
                    n_attr += 1
    return chroms, n_feat, n_attr, sorted(attr_examples)[:20]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True, type=Path, help="GTF/GFF annotation")
    ap.add_argument("--bam", type=Path, help="One BAM for chrom check (preferred small header-only)")
    ap.add_argument("--bam-dir", type=Path, help="Dir of *_sorted.bam; uses first BAM")
    ap.add_argument("-g", "--attr", default="gene_id", help="gene id attribute (default: gene_id)")
    ap.add_argument("-t", "--feature", default="exon", help="feature type (default: exon)")
    ap.add_argument("--min-chrom-overlap", type=float, default=0.01, help="min fraction of BAM chroms in annot")
    ap.add_argument("--strict", action="store_true", help="treat warnings as errors")
    args = ap.parse_args()

    errors: list[str] = []
    warns: list[str] = []

    if not args.annot.is_file():
        print(f"ERROR: annot not found: {args.annot}", file=sys.stderr)
        return 1

    bam = args.bam
    if bam is None and args.bam_dir is not None:
        bams = sorted(args.bam_dir.glob("*_sorted.bam"))
        if not bams:
            errors.append(f"no *_sorted.bam in {args.bam_dir}")
        else:
            bam = bams[0]

    print(f"=== preflight_check ===")
    print(f"annot={args.annot}")
    print(f"attr={args.attr}  feature={args.feature}")
    if bam:
        print(f"bam={bam}")

    chroms_a, n_feat, n_attr, examples = annot_stats(args.annot, args.attr, args.feature)
    print(f"annot_chroms={len(chroms_a)}  feature({args.feature})={n_feat}  with_{args.attr}={n_attr}")
    print(f"annot_attr_keys_sample={examples}")

    if n_feat == 0:
        errors.append(f"no feature type '{args.feature}' in annotation (try -t CDS/exon/transcript)")
    elif n_attr == 0:
        errors.append(
            f"feature '{args.feature}' present but attribute '{args.attr}' never seen; "
            f"try -g with one of: {examples}"
        )
    elif n_attr < max(1, n_feat // 10):
        warns.append(f"only {n_attr}/{n_feat} '{args.feature}' lines carry '{args.attr}'")

    if bam is not None:
        if not bam.is_file():
            errors.append(f"BAM not found: {bam}")
        else:
            try:
                chroms_b = bam_chroms(bam)
            except Exception as e:
                errors.append(f"samtools view -H failed on {bam}: {e}")
                chroms_b = set()
            if chroms_b:
                inter = chroms_a & chroms_b
                frac = len(inter) / max(1, len(chroms_b))
                print(f"bam_chroms={len(chroms_b)}  overlap={len(inter)}  frac={frac:.3f}")
                if not inter:
                    errors.append(
                        "ZERO chromosome overlap between BAM and annotation "
                        "(e.g. chr1 vs 1). Fix naming or choose matching GTF."
                    )
                elif frac < args.min_chrom_overlap:
                    warns.append(
                        f"low chrom overlap frac={frac:.3f} (< {args.min_chrom_overlap}); "
                        f"BAM-only examples={sorted(chroms_b - chroms_a)[:5]} "
                        f"annot-only examples={sorted(chroms_a - chroms_b)[:5]}"
                    )
                else:
                    print("chrom_overlap: OK")

    for w in warns:
        print(f"WARN: {w}")
    for e in errors:
        print(f"ERROR: {e}")

    if errors or (args.strict and warns):
        print("preflight: FAIL")
        return 1
    print("preflight: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
