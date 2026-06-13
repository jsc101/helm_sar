"""Tests for SMILES → HELM fragmentation (smiles_to_helm.py)."""
import logging
logging.disable(logging.WARNING)

import pytest
from rdkit import Chem

from scripts.smiles_to_helm import (
    _find_backbone_amide_bonds,
    _assign_rgroups,
    _try_add_r3,
    fragment_smiles,
    smiles_to_helm,
)
from monomer_db.monomer_db import MonomerDB


# ── fixtures ──────────────────────────────────────────────────────────────────

GLY_ALA = 'NCC(=O)N[C@@H](C)C(=O)O'           # Gly-Ala linear
ALA3    = 'N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O'  # Ala-Ala-Ala

# Gramicidin S: head-to-tail cyclic decapeptide (PubChem CID 73357)
GRAM_S = (
    'CC(C)[C@H]1C(=O)N[C@H](CCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H]'
    '(Cc2ccccc2)C(=O)N2CCC[C@H]2C(=O)N[C@@H](C(C)C)C(=O)N[C@H]'
    '(CCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](Cc2ccccc2)C(=O)N2CCC'
    '[C@H]2C(=O)N1'
)

# Cyclosporin A: cyclic 11-mer with N-methyl AAs (PubChem CID 5284373)
CSA = (
    'C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)'
    'N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)'
    'C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)'
    '[C@@H](CC(C)C)C(=O)N(C)[C@@H](C(C)C)C(=O)N1C'
)


# ── backbone bond detection ────────────────────────────────────────────────────

def test_backbone_bonds_gly_ala():
    mol = Chem.MolFromSmiles(GLY_ALA)
    bonds = _find_backbone_amide_bonds(mol)
    assert len(bonds) == 1


def test_backbone_bonds_ala3():
    mol = Chem.MolFromSmiles(ALA3)
    bonds = _find_backbone_amide_bonds(mol)
    assert len(bonds) == 2


def test_backbone_bonds_gram_s():
    mol = Chem.MolFromSmiles(GRAM_S)
    bonds = _find_backbone_amide_bonds(mol)
    assert len(bonds) == 10, f'expected 10, got {len(bonds)}'


def test_backbone_bonds_csa_nmethyl():
    """N-methyl backbone amides (D2=3) are detected in CsA."""
    mol = Chem.MolFromSmiles(CSA)
    bonds = _find_backbone_amide_bonds(mol)
    assert len(bonds) == 11, f'expected 11, got {len(bonds)}'


# ── R3 fallback ───────────────────────────────────────────────────────────────

def test_r3_adds_to_primary_amine():
    """Orn-like fragment: free sidechain NH2 → [*:3] neighbour."""
    frag = 'NCCC[C@@H](N[*:1])C(=O)[*:2]'
    r3 = _try_add_r3(frag)
    assert r3 is not None
    assert '[*:3]' in r3
    # N atom must still be present (not replaced)
    assert 'N' in r3


def test_r3_does_not_touch_backbone_n():
    """Backbone [*:1] N must not be converted to R3."""
    frag = '[*:1]N[C@@H](CC)C(=O)[*:2]'   # Abu
    r3 = _try_add_r3(frag)
    assert r3 is None   # no free sidechain primary amine


# ── fragment_smiles ────────────────────────────────────────────────────────────

def test_fragment_gly_ala_inner():
    """Middle residue of AGA tripeptide is identified as Gly."""
    smi = 'N[C@@H](C)C(=O)NCC(=O)N[C@@H](C)C(=O)O'  # Ala-Gly-Ala
    frags = fragment_smiles(smi)
    symbols = [f.symbol for f in frags]
    assert 'G' in symbols


def test_fragment_gram_s_count():
    """Gramicidin S produces exactly 10 fragments."""
    frags = fragment_smiles(GRAM_S)
    assert len(frags) == 10


def test_fragment_gram_s_known_residues():
    """Val, Leu, D-Phe, Pro are all resolved in Gramicidin S."""
    frags = fragment_smiles(GRAM_S)
    syms = [f.symbol for f in frags]
    for expected in ['V', 'L', 'dF', 'P']:
        assert expected in syms, f'{expected} not found in {syms}'


def test_fragment_csa_all_resolved():
    """All 11 CsA residues are resolved (including MeBmt and N-methyl AAs)."""
    frags = fragment_smiles(CSA)
    assert len(frags) == 11
    unknowns = [f for f in frags if f.symbol is None]
    assert len(unknowns) == 0, f'Unresolved: {[u.smiles[:40] for u in unknowns]}'


# ── smiles_to_helm ─────────────────────────────────────────────────────────────

def test_smiles_to_helm_gram_s_length():
    """Gramicidin S → HELM has 10 residues in PEPTIDE1, no new monomers."""
    import re
    helm, new_mons = smiles_to_helm(GRAM_S, cyclic=True)
    assert helm is not None
    assert len(new_mons) == 0
    m = re.search(r'PEPTIDE1\{([^}]+)\}', helm)
    assert m is not None
    assert len(m.group(1).split('.')) == 10


def test_smiles_to_helm_csa_cyclic_closure():
    """CsA SMILES→HELM includes head-to-tail cyclization."""
    helm, new_mons = smiles_to_helm(CSA, cyclic=True)
    assert helm is not None
    assert '1:R1-11:R2' in helm


def test_smiles_to_helm_csa_exact_sequence():
    """CsA SMILES→HELM recovers the correct 11-residue sequence by canonical SMILES."""
    from scripts.roundtrip_benchmark import _compare_sequences
    db = MonomerDB()
    known_helm = (
        'PEPTIDE1{[MeBmt].[Abu].[Sar].[MeLeu].V.[MeLeu].A.[dA].[MeLeu].[MeLeu].[MeVal]}'
        '$PEPTIDE1,PEPTIDE1,1:R1-11:R2$$V2.0'
    )
    recovered_helm, new_mons = smiles_to_helm(CSA, cyclic=True, db=db)
    assert recovered_helm is not None
    assert len(new_mons) == 0
    match = _compare_sequences(known_helm, recovered_helm, db)
    assert match == 'exact', f'Sequence match was: {match}'


def test_smiles_to_helm_auto_registers_unknowns():
    """Unknown fragments get placeholder symbols and the HELM is complete."""
    import re
    # Use a peptide with a disulfide (produces unfragmentable bridge)
    from scripts.helm_parser import HELMParser
    helm_pram = (
        'PEPTIDE1{K.C.N.T.A}$PEPTIDE1,PEPTIDE1,2:R3-5:R3$$V2.0'
    )
    smi = HELMParser.parse(helm_pram).to_smiles()
    db = MonomerDB()
    helm, new_mons = smiles_to_helm(smi, cyclic=False, db=db)
    assert helm is not None
    assert 'UNKNOWN' not in helm
    # All new monomers have content-addressed symbols
    for nm in new_mons:
        assert nm['symbol'].startswith('UNK_')
        assert nm['needs_review'] is True
    # Second call with same DB: zero new monomers (reuses registered symbols)
    helm2, new2 = smiles_to_helm(smi, cyclic=False, db=db)
    assert len(new2) == 0
    assert helm == helm2


def test_smiles_to_helm_new_monomer_log(tmp_path):
    """Auto-registered monomers are persisted to the JSON log file."""
    import json
    from scripts.helm_parser import HELMParser
    helm_src = 'PEPTIDE1{K.C.N.T.A}$PEPTIDE1,PEPTIDE1,2:R3-5:R3$$V2.0'
    smi = HELMParser.parse(helm_src).to_smiles()
    log = tmp_path / 'new_monomers.json'
    db = MonomerDB()
    _, new_mons = smiles_to_helm(smi, cyclic=False, db=db, new_monomer_log=log)
    if new_mons:
        assert log.exists()
        registry = json.loads(log.read_text())
        assert isinstance(registry, dict)
        assert len(registry) == len(new_mons)
        assert all(v.get('needs_review') for v in registry.values())
        assert all(v.get('assigned_symbol') is None for v in registry.values())
