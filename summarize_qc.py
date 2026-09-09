#!/usr/bin/env python3
"""Summarize RNA-seq pipeline QC into tables under 10.qc/.

Inputs (auto-detected under cwd or --outdir parent):
  08.align_stats/*_{raw,unique}.flagstat.txt
  09.featureCounts/gene_counts.txt.summary

Outputs:
  10.qc/align_flagstat_summary.tsv
  10.qc/featurecounts_summary.tsv
  10.qc/qc_overview.tsv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_flagstat(path: Path) -> dict:
    text = path.read_text()
    d: dict = {"file": str(path)}
    # samtools flagstat classic lines
    m = re.search(r"^(\d+) \+ \d+ in total", text, re.M)
    if m:
        d["total"] = int(m.group(1))
    m = re.search(r"^(\d+) \+ \d+ mapped \(([\d.]+)", text, re.M)
    if m:
        d["mapped"] = int(m.group(1))
        d["mapped_pct"] = float(m.group(2))
    m = re.search(r"^(\d+) \+ \d+ properly paired \(([\d.]+)", text, re.M)
    if m:
        d["properly_paired"] = int(m.group(1))
        d["properly_paired_pct"] = float(m.group(2))
    m = re.search(r"^(\d+) \+ \d+ singletons \(([\d.]+)", text, re.M)
    if m:
        d["singletons"] = int(m.group(1))
        d["singletons_pct"] = float(m.group(2))
    return d


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="10.qc", help="QC output dir (default: 10.qc)")
    ap.add_argument("--align-dir", default="08.align_stats")
    ap.add_argument("--fc-summary", default="09.featureCounts/gene_counts.txt.summary")
    ap.add_argument("--min-assigned-pct", type=float, default=50.0, help="warn if Assigned%% below this")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- align flagstat ----
    align_rows = []
    align_dir = Path(args.align_dir)
    if align_dir.is_dir():
        for p in sorted(align_dir.glob("*.flagstat.txt")):
            name = p.name.replace(".flagstat.txt", "")
            kind = "unique" if name.endswith("_unique") else ("raw" if name.endswith("_raw") else "other")
            sample = name
            for suf in ("_unique", "_raw"):
                if sample.endswith(suf):
                    sample = sample[: -len(suf)]
                    break
            d = parse_flagstat(p)
            align_rows.append(
                {
                    "sample": sample,
                    "stage": kind,
                    "total": d.get("total", ""),
                    "mapped": d.get("mapped", ""),
                    "mapped_pct": d.get("mapped_pct", ""),
                    "properly_paired": d.get("properly_paired", ""),
                    "properly_paired_pct": d.get("properly_paired_pct", ""),
                    "singletons": d.get("singletons", ""),
                    "singletons_pct": d.get("singletons_pct", ""),
                    "source": p.name,
                }
            )
    align_cols = [
        "sample",
        "stage",
        "total",
        "mapped",
        "mapped_pct",
        "properly_paired",
        "properly_paired_pct",
        "singletons",
        "singletons_pct",
        "source",
    ]
    write_tsv(out / "align_flagstat_summary.tsv", align_rows, align_cols)
    print(f"wrote {out / 'align_flagstat_summary.tsv'} rows={len(align_rows)}")

    # ---- featureCounts summary ----
    fc_rows = []
    fc_path = Path(args.fc_summary)
    if fc_path.is_file():
        lines = fc_path.read_text().splitlines()
        header = lines[0].split("\t")
        samples = []
        for h in header[1:]:
            s = Path(h).name
            s = s.replace("_sorted.bam", "").replace(".bam", "")
            samples.append(s)
        status_map: dict[str, list[str]] = {}
        for line in lines[1:]:
            parts = line.split("\t")
            status_map[parts[0]] = parts[1:]
        for i, sample in enumerate(samples):
            assigned = float(status_map.get("Assigned", ["0"] * len(samples))[i])
            nofeat = float(status_map.get("Unassigned_NoFeatures", ["0"] * len(samples))[i])
            ambig = float(status_map.get("Unassigned_Ambiguity", ["0"] * len(samples))[i])
            multi = float(status_map.get("Unassigned_MultiMapping", ["0"] * len(samples))[i])
            total = assigned + nofeat + ambig + multi
            # also sum all numeric statuses for denom
            all_sum = 0.0
            for vals in status_map.values():
                try:
                    all_sum += float(vals[i])
                except ValueError:
                    pass
            denom = all_sum if all_sum > 0 else total
            pct = 100.0 * assigned / denom if denom > 0 else 0.0
            warn = "LOW_ASSIGNED" if pct < args.min_assigned_pct else "OK"
            fc_rows.append(
                {
                    "sample": sample,
                    "Assigned": int(assigned),
                    "Unassigned_NoFeatures": int(nofeat),
                    "Unassigned_Ambiguity": int(ambig),
                    "Unassigned_MultiMapping": int(multi),
                    "Assigned_pct": f"{pct:.2f}",
                    "status": warn,
                }
            )
        fc_cols = [
            "sample",
            "Assigned",
            "Unassigned_NoFeatures",
            "Unassigned_Ambiguity",
            "Unassigned_MultiMapping",
            "Assigned_pct",
            "status",
        ]
        write_tsv(out / "featurecounts_summary.tsv", fc_rows, fc_cols)
        print(f"wrote {out / 'featurecounts_summary.tsv'} rows={len(fc_rows)}")
        low = [r for r in fc_rows if r["status"] != "OK"]
        if low:
            print(f"WARN: {len(low)} sample(s) Assigned_pct < {args.min_assigned_pct}")
            for r in low:
                print(f"  {r['sample']}: {r['Assigned_pct']}%")
    else:
        print(f"NOTE: no featureCounts summary at {fc_path}")

    # ---- overview join on unique align + fc ----
    unique = {r["sample"]: r for r in align_rows if r["stage"] == "unique"}
    overview = []
    samples = sorted(set(unique) | {r["sample"] for r in fc_rows})
    for s in samples:
        u = unique.get(s, {})
        f = next((r for r in fc_rows if r["sample"] == s), {})
        overview.append(
            {
                "sample": s,
                "unique_mapped_pct": u.get("mapped_pct", ""),
                "unique_properly_paired_pct": u.get("properly_paired_pct", ""),
                "fc_Assigned_pct": f.get("Assigned_pct", ""),
                "fc_status": f.get("status", ""),
            }
        )
    write_tsv(
        out / "qc_overview.tsv",
        overview,
        ["sample", "unique_mapped_pct", "unique_properly_paired_pct", "fc_Assigned_pct", "fc_status"],
    )
    print(f"wrote {out / 'qc_overview.tsv'} rows={len(overview)}")


if __name__ == "__main__":
    main()
