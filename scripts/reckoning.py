"""Reckoning — fragmentation rule validation harness.

Loads data/ground_truth.json and scores the current smiles_to_helm engine
against every curated SMILES↔HELM pair. Reports coverage per compound and
flags failures for rule authoring.

Usage
-----
    python scripts/reckoning.py                    # score all pairs, print table
    python scripts/reckoning.py --out results.csv  # also write CSV
    python scripts/reckoning.py --id cyclosporin_a # single compound

Output columns
--------------
    id, name, n_residues, n_known, n_unknown, coverage,
    expected_coverage, status, new_symbols, notes
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import logging
logging.disable(logging.WARNING)

from monomer_db.monomer_db import MonomerDB
from scripts.smiles_to_helm import fragment_smiles, smiles_to_helm


_GT_PATH = ROOT / 'data' / 'ground_truth.json'
_RULES_PATH = ROOT / 'data' / 'fragmentation_rules.json'


# ── Core scoring ───────────────────────────────────────────────────────────────

def score_pair(gt: dict, db: MonomerDB) -> dict:
    """Run fragmentation on one ground-truth entry and return a result dict."""
    smiles = gt['smiles']
    cyclic = gt['cyclic']
    known_helm = gt.get('helm')
    expected_cov = gt.get('expected_coverage')
    limitation = gt.get('known_limitation')

    helm_out, new_mons = smiles_to_helm(smiles, cyclic=cyclic, db=db)

    if helm_out is None:
        return {
            'id':               gt['id'],
            'name':             gt['name'],
            'n_residues':       gt['n_residues'],
            'n_known':          0,
            'n_unknown':        0,
            'coverage':         0.0,
            'expected_coverage': expected_cov,
            'status':           'fragmentation_failed',
            'new_symbols':      '',
            'known_limitation': limitation or '',
            'helm_out':         '',
        }

    frags = fragment_smiles(smiles, db)
    n_total = len(frags)
    n_unknown = len(new_mons)
    n_known = n_total - n_unknown
    coverage = round(n_known / n_total, 3) if n_total else 0.0

    if limitation:
        status = f'known_limitation:{limitation}'
    elif expected_cov is not None and coverage >= expected_cov:
        status = 'pass'
    elif expected_cov is not None:
        status = f'partial ({coverage:.0%} vs expected {expected_cov:.0%})'
    else:
        status = 'unvalidated'

    new_syms = ', '.join(m['symbol'] for m in new_mons)

    return {
        'id':                gt['id'],
        'name':              gt['name'],
        'n_residues':        gt['n_residues'],
        'n_known':           n_known,
        'n_unknown':         n_unknown,
        'coverage':          coverage,
        'expected_coverage': expected_cov,
        'status':            status,
        'new_symbols':       new_syms,
        'known_limitation':  limitation or '',
        'helm_out':          helm_out,
    }


# ── New-compound mode ──────────────────────────────────────────────────────────

def score_new_compound(smiles: str, name: str, cyclic: bool = False,
                       db: MonomerDB | None = None,
                       pending_log: Path | None = None) -> dict:
    """Score a novel compound with no pre-validated HELM.

    Returns the generated HELM string and any pending monomers.
    Optionally writes UNK entries to pending_log for chemist review.
    """
    if db is None:
        db = MonomerDB()
    helm, new_mons = smiles_to_helm(smiles, cyclic=cyclic, db=db,
                                    new_monomer_log=pending_log)
    frags = fragment_smiles(smiles, db)
    n_total = len(frags)
    n_unknown = len(new_mons)
    n_known = n_total - n_unknown
    coverage = round(n_known / n_total, 3) if n_total else 0.0
    return {
        'name':        name,
        'helm':        helm,
        'n_residues':  n_total,
        'n_known':     n_known,
        'n_unknown':   n_unknown,
        'coverage':    coverage,
        'new_monomers': new_mons,
    }


# ── Rule-set health ────────────────────────────────────────────────────────────

def check_rule_health(rules_path: Path = _RULES_PATH) -> None:
    """Print a summary of active vs planned rules from the rule library."""
    with open(rules_path) as f:
        lib = json.load(f)
    rules = lib['rules']
    active   = [r for r in rules if r['status'] == 'active']
    planned  = [r for r in rules if r['status'] == 'planned']
    deprecated = [r for r in rules if r['status'] == 'deprecated']
    print(f'\nRule library v{lib["version"]}')
    print(f'  Active:     {len(active)}  ({", ".join(r["id"] for r in active)})')
    print(f'  Planned:    {len(planned)}  ({", ".join(r["id"] for r in planned)})')
    print(f'  Deprecated: {len(deprecated)}')
    print()
    for r in planned:
        vcount = len(r.get('validated_on', []))
        print(f'  [{r["layer"]}.{r["priority"]}] {r["id"]} — {r["description"][:60]}')
        print(f'           validated_on: {vcount} pairs | needed for: {r["notes"][:60]}')


# ── Batch runner ───────────────────────────────────────────────────────────────

def run_reckoning(
    gt_path: Path = _GT_PATH,
    filter_id: str | None = None,
    out_csv: Path | None = None,
    verbose: bool = True,
) -> list[dict]:
    with open(gt_path) as f:
        ground_truth = json.load(f)

    if filter_id:
        ground_truth = [g for g in ground_truth if g['id'] == filter_id]
        if not ground_truth:
            print(f'No entry with id={filter_id!r}')
            return []

    db = MonomerDB()
    results = []

    for gt in ground_truth:
        r = score_pair(gt, db)
        results.append(r)

    if verbose:
        _print_table(results)

    if out_csv:
        _write_csv(results, out_csv)
        print(f'\nWritten → {out_csv}')

    return results


def _print_table(results: list[dict]) -> None:
    print()
    header = f'{"Name":<22} {"Res":>4} {"Known":>5} {"UNK":>4} {"Cov":>6}  Status'
    print(header)
    print('─' * len(header))
    for r in results:
        cov_str = f'{r["coverage"]:.0%}'
        unk_note = f'  [{r["new_symbols"]}]' if r['new_symbols'] else ''
        print(f'{r["name"]:<22} {r["n_residues"]:>4} {r["n_known"]:>5} '
              f'{r["n_unknown"]:>4} {cov_str:>6}  {r["status"]}{unk_note}')
    print()
    n_pass = sum(1 for r in results if r['status'] == 'pass')
    n_lim  = sum(1 for r in results if 'known_limitation' in r['status'])
    n_fail = len(results) - n_pass - n_lim
    print(f'Summary: {n_pass} pass  {n_lim} known-limitation  {n_fail} fail  '
          f'(of {len(results)} total)')


def _write_csv(results: list[dict], path: Path) -> None:
    fields = ['id', 'name', 'n_residues', 'n_known', 'n_unknown', 'coverage',
              'expected_coverage', 'status', 'new_symbols', 'known_limitation']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(results)


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reckoning — fragmentation validation harness')
    parser.add_argument('--id',      help='Score a single ground-truth entry by ID')
    parser.add_argument('--out',     help='Write results to CSV file')
    parser.add_argument('--rules',   action='store_true', help='Show rule-set health report')
    parser.add_argument('--new',     help='Score a novel SMILES string (not in ground truth)')
    parser.add_argument('--name',    default='novel_compound', help='Name for --new compound')
    parser.add_argument('--cyclic',  action='store_true', help='Treat --new compound as cyclic')
    args = parser.parse_args()

    if args.rules:
        check_rule_health()

    if args.new:
        r = score_new_compound(args.new, args.name, cyclic=args.cyclic,
                               pending_log=ROOT / 'data' / 'pending_monomers.json')
        print(f'\n{r["name"]}')
        print(f'  HELM:      {r["helm"]}')
        print(f'  Residues:  {r["n_residues"]}  ({r["n_known"]} known, {r["n_unknown"]} UNK)')
        print(f'  Coverage:  {r["coverage"]:.0%}')
        if r['new_monomers']:
            print(f'  Pending:   {[m["symbol"] for m in r["new_monomers"]]}')
            print(f'             → data/pending_monomers.json')
    elif not args.rules:
        run_reckoning(
            filter_id=args.id,
            out_csv=Path(args.out) if args.out else None,
        )
