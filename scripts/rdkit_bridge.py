"""
RDKit bridge for peptide informatics.

Public API
----------
monomer_to_mol(entry)          → Chem.Mol   (R-groups capped, ready for descriptors)
mol_to_descriptors(mol)        → dict        (MW, LogP, TPSA, HBD, HBA, RotBonds, Rings, QED)
validate_roundtrip(smiles, mol) → dict       (tanimoto, substructure_match, pass)
"""

from __future__ import annotations

import re
import logging
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs, QED
from rdkit.Chem import RDKFingerprint

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# R-group capping
# ---------------------------------------------------------------------------

def _cap_from_capgroup(cap_smiles: str) -> str:
    """
    Extract the cap atom/fragment from a capGroupSmiles string.

    capGroupSmiles examples:
      '[*:1][H]'   → '[H]'
      'O[*:2]'     → 'O'
      '[*:3][H]'   → '[H]'
    """
    cap = re.sub(r'\[\*:\d+\]', '', cap_smiles).strip()
    return cap if cap else '[H]'


def _cap_rgroups(smiles: str, rgroups: list[dict]) -> str:
    """
    Replace every [*:N] attachment point in *smiles* with its cap atom.

    Uses capGroupSmiles from the rgroups list; falls back to [H] for any
    unmapped attachment points.
    """
    caps: dict[int, str] = {}
    for rg in rgroups:
        m = re.match(r'R(\d+)', rg.get('label', ''))
        if not m:
            continue
        caps[int(m.group(1))] = _cap_from_capgroup(rg.get('capGroupSmiles', '') or '[H]')

    result = smiles
    for n, cap in caps.items():
        result = result.replace(f'[*:{n}]', cap)

    # Any remaining unmapped attachment points → H
    result = re.sub(r'\[\*:\d+\]', '[H]', result)
    return result


# ---------------------------------------------------------------------------
# monomer_to_mol
# ---------------------------------------------------------------------------

def monomer_to_mol(entry: dict) -> Optional[Chem.Mol]:
    """
    Convert a MonomerDB entry to an RDKit Mol suitable for descriptor calculation.

    Handles two SMILES conventions:
      Pistoia   — attachment points as [H:N]/[OH:N], converted to [*:N] then capped
      CycPeptMPDB — replaced_SMILES (fully capped, no attachment markers)

    R-groups are capped using the entry's capGroupSmiles values; fallback is [H].
    """
    from monomer_db.monomer_db import _helm_to_star_smiles

    smiles = entry.get('smiles', '')
    if not smiles:
        return None

    # Convert Pistoia [H:N]/[OH:N] notation → [*:N]
    star = _helm_to_star_smiles(smiles)

    if '[*:' in star:
        rgroups = entry.get('rgroups', [])
        capped = _cap_rgroups(star, rgroups)
    else:
        capped = star  # already a complete molecule (CycPeptMPDB replaced_SMILES)

    mol = Chem.MolFromSmiles(capped)
    if mol is None:
        # fallback: strip all attachment points with H
        fallback = re.sub(r'\[\*:\d+\]', '[H]', star)
        mol = Chem.MolFromSmiles(fallback)
        if mol is not None:
            logger.debug("monomer_to_mol: used fallback capping for %r", entry.get('symbol'))

    return mol


# ---------------------------------------------------------------------------
# Descriptors
# ---------------------------------------------------------------------------

def mol_to_descriptors(mol: Chem.Mol) -> dict:
    """
    Calculate key physicochemical descriptors for a molecule.

    Returns
    -------
    dict with keys: MW, ExactMW, LogP, TPSA, HBD, HBA, RotBonds, Rings, AromaticRings, QED
    """
    if mol is None:
        return {}
    try:
        return {
            'MW':           round(Descriptors.MolWt(mol), 2),
            'ExactMW':      round(Descriptors.ExactMolWt(mol), 4),
            'LogP':         round(Descriptors.MolLogP(mol), 2),
            'TPSA':         round(rdMolDescriptors.CalcTPSA(mol), 2),
            'HBD':          rdMolDescriptors.CalcNumHBD(mol),
            'HBA':          rdMolDescriptors.CalcNumHBA(mol),
            'RotBonds':     rdMolDescriptors.CalcNumRotatableBonds(mol),
            'Rings':        rdMolDescriptors.CalcNumRings(mol),
            'AromaticRings':rdMolDescriptors.CalcNumAromaticRings(mol),
            'QED':          round(QED.qed(mol), 4),
        }
    except Exception as e:
        logger.warning("mol_to_descriptors failed: %s", e)
        return {}


# ---------------------------------------------------------------------------
# Roundtrip validation
# ---------------------------------------------------------------------------

def validate_roundtrip(input_smiles: str, computed_mol: Chem.Mol) -> dict:
    """
    Validate structural equivalence between a reference SMILES and a computed Mol.

    Returns
    -------
    dict: {tanimoto, substructure_match, pass, reason}
    """
    if computed_mol is None:
        return {'tanimoto': 0.0, 'substructure_match': False, 'pass': False,
                'reason': 'computed_mol is None'}

    ref_mol = Chem.MolFromSmiles(input_smiles)
    if ref_mol is None:
        return {'tanimoto': 0.0, 'substructure_match': False, 'pass': False,
                'reason': f'input_smiles failed parse: {input_smiles[:80]}'}

    ref_fp  = RDKFingerprint(ref_mol)
    comp_fp = RDKFingerprint(computed_mol)
    tanimoto = DataStructs.TanimotoSimilarity(ref_fp, comp_fp)

    fwd = computed_mol.HasSubstructMatch(ref_mol)
    rev = ref_mol.HasSubstructMatch(computed_mol)
    substructure_match = fwd and rev
    passed = tanimoto >= 0.99 or substructure_match

    return {
        'tanimoto':           round(tanimoto, 4),
        'substructure_match': substructure_match,
        'pass':               passed,
        'reason':             'ok' if passed else f'tanimoto={tanimoto:.3f}',
    }
