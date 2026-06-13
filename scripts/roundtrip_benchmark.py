"""Roundtrip benchmark: HELM ↔ SMILES.

Benchmark peptides
------------------
Gramicidin S      — head-to-tail cyclic decapeptide (confirmed PubChem CID 73357)
Cyclosporin A     — head-to-tail cyclic 11-mer with N-methyl AAs (confirmed PubChem CID 5284373)
Somatostatin-14   — cyclic disulfide (self-consistency: HELM→SMILES→HELM)
Octreotide        — somatostatin analog, Cys2-Cys7 disulfide, D-Thr-ol C-terminus
Pramlintide       — amylin analog 37-mer with Cys2-Cys7 disulfide

For each peptide:
  HELM→SMILES : parse HELM, assemble SMILES, compute Morgan Tanimoto vs
                reference SMILES (PubChem where confirmed, else self-generated).
  SMILES→HELM : fragment reference SMILES, match monomers by canonical SMILES
                (symbol aliases meG/Sar are normalised via DB lookup).

Results written to data/roundtrip_benchmark.csv.
"""
from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

logging.disable(logging.WARNING)

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.helm_parser import HELMParser
from scripts.smiles_to_helm import fragment_smiles, smiles_to_helm
from monomer_db.monomer_db import MonomerDB

# ── Reference data ─────────────────────────────────────────────────────────────
# PubChem SMILES marked with source. 'self' = generated from HELM (no external ref).

BENCHMARKS: list[dict] = [
    {
        'name': 'Gramicidin S',
        'helm': (
            'PEPTIDE1{V.[Orn].L.[dF].P.V.[Orn].L.[dF].P}'
            '$PEPTIDE1,PEPTIDE1,1:R1-10:R2$$V2.0'
        ),
        # PubChem CID 73357 — confirmed Tanimoto=1.0
        'ref_smiles': (
            'CC(C)[C@H]1C(=O)N[C@H](CCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H]'
            '(Cc2ccccc2)C(=O)N2CCC[C@H]2C(=O)N[C@@H](C(C)C)C(=O)N[C@H]'
            '(CCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](Cc2ccccc2)C(=O)N2CCC'
            '[C@H]2C(=O)N1'
        ),
        'ref_source': 'PubChem CID 73357',
        'linear': False,
        'note': 'head-to-tail cyclic decapeptide',
    },
    {
        'name': 'Cyclosporin A',
        'helm': (
            'PEPTIDE1{[MeBmt].[Abu].[Sar].[MeLeu].V.[MeLeu].A.[dA].[MeLeu].[MeLeu].[MeVal]}'
            '$PEPTIDE1,PEPTIDE1,1:R1-11:R2$$V2.0'
        ),
        # PubChem CID 5284373 — confirmed Tanimoto=1.0
        'ref_smiles': (
            'C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)'
            'N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)'
            'C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)'
            '[C@@H](CC(C)C)C(=O)N(C)[C@@H](C(C)C)C(=O)N1C'
        ),
        'ref_source': 'PubChem CID 5284373',
        'linear': False,
        'note': 'cyclic 11-mer, N-methyl AAs, no disulfide',
    },
    {
        'name': 'Somatostatin-14',
        'helm': (
            'PEPTIDE1{A.G.C.K.N.F.F.W.K.T.F.T.S.C}'
            '$PEPTIDE1,PEPTIDE1,3:R3-14:R3$$V2.0'
        ),
        'ref_smiles': None,   # generated self-consistently below
        'ref_source': 'self',
        'linear': False,
        'note': 'cyclic disulfide Cys3-Cys14',
    },
    {
        'name': 'Octreotide',
        'helm': (
            'PEPTIDE1{[dF].C.F.[dW].K.T.C.[D-Thr-ol]}'
            '$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$V2.0'
        ),
        'ref_smiles': None,   # generated self-consistently below
        'ref_source': 'self',
        'linear': False,
        'note': 'Cys2-Cys7 disulfide, D-Thr-ol C-terminus',
    },
    {
        'name': 'Pramlintide',
        'helm': (
            'PEPTIDE1{K.C.N.T.A.T.C.A.T.Q.R.L.A.N.F.L.V.H.S.S.N.N.F.G.P.I.L.P.P.T.N.V.G.S.N.T.Y}'
            '$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$V2.0'
        ),
        'ref_smiles': None,
        'ref_source': 'self',
        'linear': True,
        'note': 'amylin analog 37-mer, Cys2-Cys7 disulfide',
    },
]

# ── Helpers ────────────────────────────────────────────────────────────────────

def _fp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def _tanimoto(smi_a: str, smi_b: str) -> float | None:
    fa, fb = _fp(smi_a), _fp(smi_b)
    if fa is None or fb is None:
        return None
    return DataStructs.TanimotoSimilarity(fa, fb)


def _helm_to_smiles(helm: str, cap_c: str = 'amide') -> str | None:
    try:
        obj = HELMParser.parse(helm)
        return obj.to_smiles(cap_c=cap_c)
    except Exception as exc:
        print(f'  HELM→SMILES error: {exc}')
        return None


def _helm_canonical_seq(helm: str, db: MonomerDB) -> list[str]:
    """Return sequence as list of canonical SMILES per position (normalises symbol aliases)."""
    obj = HELMParser.parse(helm)
    chain = obj.get_chain()
    if chain is None:
        return []
    result = []
    for m in chain['monomers']:
        entry = m.get('entry')
        if entry:
            mol = Chem.MolFromSmiles(entry['smiles'].replace('[H:', '[*:').replace('[OH:', '[*:'))
            result.append(Chem.MolToSmiles(mol) if mol else m['symbol'])
        else:
            result.append(m['symbol'])
    return result


def _compare_sequences(helm_known: str, helm_recovered: str, db: MonomerDB) -> str:
    """Compare two HELM strings by canonical monomer SMILES (ignores symbol aliases)."""
    try:
        seq_k = _helm_canonical_seq(helm_known, db)
        seq_r = _helm_canonical_seq(helm_recovered, db)
        if seq_k == seq_r:
            return 'exact'
        if len(seq_k) == len(seq_r):
            n_match = sum(a == b for a, b in zip(seq_k, seq_r))
            return f'partial_{n_match}/{len(seq_k)}'
        return 'mismatch'
    except Exception:
        return 'error'


# ── Main benchmark ─────────────────────────────────────────────────────────────

def run_benchmark(out_csv: Path | None = None) -> list[dict]:
    db = MonomerDB()
    rows: list[dict] = []

    # Fill in self-generated reference SMILES
    for bp in BENCHMARKS:
        if bp['ref_smiles'] is None:
            bp['ref_smiles'] = _helm_to_smiles(bp['helm'])

    for bp in BENCHMARKS:
        name = bp['name']
        helm = bp['helm']
        ref_smi = bp['ref_smiles']
        note = bp['note']
        ref_src = bp['ref_source']

        print(f'\n── {name} ({note}) ──')

        # ── HELM → SMILES ──────────────────────────────────────────────────
        gen_smi = _helm_to_smiles(helm)
        if gen_smi:
            print(f'  HELM→SMILES: {gen_smi[:80]}{"..." if len(gen_smi)>80 else ""}')
        else:
            print('  HELM→SMILES: FAILED')

        tanimoto_h2s: float | str = 'n/a'
        if gen_smi and ref_smi:
            t = _tanimoto(gen_smi, ref_smi)
            tanimoto_h2s = round(t, 4) if t is not None else 'error'
            print(f'  Tanimoto vs {ref_src}: {tanimoto_h2s}')

        # ── SMILES → HELM ──────────────────────────────────────────────────
        s2h_match: str = 'n/a'
        s2h_unknown: int = 0
        if ref_smi:
            recovered_helm, unknowns = smiles_to_helm(
                ref_smi,
                cyclic=not bp['linear'],
                db=db,
                new_monomer_log=ROOT / 'data' / 'new_monomers.json',
            )
            s2h_unknown = len(unknowns)

            if recovered_helm:
                s2h_match = _compare_sequences(helm, recovered_helm, db)
                print(f'  SMILES→HELM recovered: {recovered_helm[:80]}')
                print(f'  Sequence match: {s2h_match}  unknowns: {s2h_unknown}')
            else:
                s2h_match = 'failed'
                print(f'  SMILES→HELM: FAILED  unknowns: {s2h_unknown}')

        rows.append({
            'name': name,
            'note': note,
            'helm_to_smiles': 'ok' if gen_smi else 'failed',
            'ref_source': ref_src,
            'tanimoto': tanimoto_h2s,
            'smiles_to_helm': s2h_match,
            'unknown_monomers': s2h_unknown,
        })

    if out_csv:
        fields = ['name', 'note', 'helm_to_smiles', 'ref_source', 'tanimoto',
                  'smiles_to_helm', 'unknown_monomers']
        with open(out_csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f'\nWritten → {out_csv}')

    return rows


if __name__ == '__main__':
    out = ROOT / 'data' / 'roundtrip_benchmark.csv'
    run_benchmark(out_csv=out)
