#!/usr/bin/env python3
# peptide_design/cli.py
"""
Read a CSV or XLSX in the peptide-design format and emit HELM + validation status.

Usage:
    python -m peptide_design.cli input.csv
    python -m peptide_design.cli input.xlsx
    python -m peptide_design.cli input.csv --stdout
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

from peptide_design.core import parse_row
from peptide_design.generator import generate_helm
from peptide_design.validator import validate_helm
from scripts.table_io import read_table


def process(rows: list[dict]) -> list[dict]:
    out = []
    for row in rows:
        compound = parse_row(row)
        helm = generate_helm(compound)
        ok, msg = validate_helm(helm)
        result = dict(row)
        result["HELM"] = helm
        result["HELM_status"] = "OK" if ok else f"ERROR: {msg}"
        out.append(result)
    return out


def _write_csv(path: Path, orig_headers: list[str], rows: list[dict]) -> None:
    extra = ["HELM", "HELM_status"]
    fieldnames = list(orig_headers) + [h for h in extra if h not in orig_headers]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description="Peptide design -> HELM")
    p.add_argument("input", help="CSV or XLSX input file")
    p.add_argument("--stdout", action="store_true", help="Print HELM strings to stdout")
    args = p.parse_args(argv)

    path = Path(args.input)
    headers, rows = read_table(path)

    results = process(rows)

    if args.stdout:
        for r in results:
            print(f"{r.get('Name','?')}\t{r['HELM']}\t{r['HELM_status']}")
        return

    out_path = path.with_stem(path.stem + "_helm").with_suffix(".csv")
    _write_csv(out_path, headers, results)
    print(f"Written {len(results)} rows -> {out_path}")
    errors = [r for r in results if r["HELM_status"] != "OK"]
    if errors:
        print(f"  {len(errors)} validation errors:")
        for r in errors:
            print(f"    {r.get('Name','?')}: {r['HELM_status']}")


if __name__ == "__main__":
    main()
