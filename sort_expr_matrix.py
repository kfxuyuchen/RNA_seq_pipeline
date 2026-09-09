#!/usr/bin/env python3
"""Sort expression matrices by first-column ID (gene/transcript), keep header.

Sample columns are left as-is (pipeline already sorts sample names when building
sample_lst / BAM list). Safe to re-run (idempotent).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def detect_delim(header: str) -> str:
    if header.count("\t") >= header.count(","):
        return "\t"
    return ","


def sort_file(path: Path) -> tuple[int, bool]:
    text = path.read_text()
    lines = text.splitlines()
    if not lines:
        raise SystemExit(f"ERROR: empty matrix {path}")
    header = lines[0]
    delim = detect_delim(header)
    body = lines[1:]
    key = lambda L: L.split(delim, 1)[0].strip('"')
    body_sorted = sorted(body, key=key)
    changed = body_sorted != body
    tmp = path.with_suffix(path.suffix + ".sorttmp")
    tmp.write_text(header + "\n" + "\n".join(body_sorted) + ("\n" if body_sorted else ""))
    tmp.replace(path)
    return len(body_sorted), changed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("matrices", nargs="+", type=Path, help="CSV/TSV matrix files")
    args = ap.parse_args()
    for p in args.matrices:
        if not p.is_file():
            print(f"ERROR: missing {p}", file=sys.stderr)
            sys.exit(1)
        n, changed = sort_file(p)
        state = "rewrote" if changed else "already-sorted"
        print(f"sort_expr_matrix: {state} rows={n} {p}")


if __name__ == "__main__":
    main()
