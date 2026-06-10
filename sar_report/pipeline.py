"""Build the aligned row set from (name, HELM) pairs."""
from __future__ import annotations

from scripts.helm_parser import HELMParser
from scripts.helm_alignment import HELMAlignment
from monomer_db.monomer_db import MonomerDB
from sar_report.mw import calc_mw


def build_data(pairs: list[tuple[str, str]], db: MonomerDB) -> tuple:
    ref_name, ref_helm = next(
        ((n, h) for n, h in pairs if 'parent' in n.lower()),
        pairs[0]
    )
    ref_obj = HELMParser.parse(ref_helm)
    ref_mw  = calc_mw(ref_obj)
    ref_row = {
        'id': ref_name, 'helm': ref_helm, 'obj': ref_obj,
        'mw': ref_mw, 'score': 1.0, 'is_ref': True,
    }

    aligner = HELMAlignment(ref_obj)
    rows = [ref_row]

    for name, helm in pairs:
        if name.lower() == ref_name.lower():
            continue
        try:
            obj = HELMParser.parse(helm)
            if len(obj.positions()) != len(ref_obj.positions()):
                rows.append({'id': name, 'helm': helm, 'obj': obj,
                             'mw': calc_mw(obj), 'score': None, 'is_ref': False})
                continue
            result = aligner.align(obj)
            rows.append({
                'id':    name,
                'helm':  helm,
                'obj':   result.rotated,
                'mw':    calc_mw(obj),
                'score': result.best_score,
                'is_ref': False,
            })
        except Exception as e:
            rows.append({'id': name, 'helm': helm, 'obj': None,
                         'mw': None, 'score': None, 'is_ref': False, 'err': str(e)})

    return ref_obj, rows


def read_numbers_helms(path: str) -> list[tuple[str, str]]:
    from scripts.table_io import read_table
    headers, rows = read_table(path)
    pairs = []
    for row in rows:
        vals = [row.get(h, '') for h in headers]
        name = vals[0] if vals else ''
        helm = vals[1] if len(vals) > 1 else ''
        if helm and 'PEPTIDE' in helm:
            pairs.append((name, helm))
    return pairs
