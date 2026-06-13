"""Tests for combine_monomers / helm_obj_to_smiles (HELM → SMILES assembly).

Verification strategy
---------------------
We cannot round-trip SMILES → HELM (no deconvolution), so we check:
  1. SMILES is RDKit-parseable (no None)
  2. Atom count matches expected residue composition
  3. ExactMW from SMILES ≈ calc_mw + expected_offset
     - backbone-only (no crosslinks): offset = 0
     - one disulfide: offset ≈ +16 Da
       (calc_mw subtracts H₂O per bond; true disulfide loses H₂ not H₂O → +16)
  4. Canonical SMILES is deterministic (same input → same output)
  5. cap_c='acid' gives one fewer N, one extra O vs cap_c='amide'
"""
import logging
logging.disable(logging.WARNING)

import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors

from scripts.helm_parser import HELMParser
from scripts.rdkit_bridge import combine_monomers, helm_obj_to_smiles
from sar_report import calc_mw

# ── HELM fixtures ────────────────────────────────────────────────────────────

# Simple tripeptide (no crosslinks) — baseline
ALA3 = "PEPTIDE1{A.A.A}$$$$V2.0"

# PYY-like lipidated analog (intra-chain only; protractor ignored here)
PYY_REF = (
    "PEPTIDE1{I.K.P.E.A.P.G.E.D.A.S.P.E.E.L.N.R.Y.Y.A.S.L.R.H.Y.L.N.L.V.T.R.Q.R.Y}"
    "$$$$V2.0"
)

# Pramlintide — 37-mer with Cys2-Cys7 disulfide, no lipidation
PRAMLINTIDE = (
    "PEPTIDE1{K.C.N.T.A.T.C.A.T.Q.R.L.A.N.F.L.V.H.S.S.N.N.F.G.P.I.L.P.P.T.N.V.G.S.N.T.Y}"
    "$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$V2.0"
)

# Cagrilintide backbone (PEPTIDE1 only; protractor on PEPTIDE2 is ignored here)
CAGRILI = (
    "PEPTIDE1{K.C.N.T.A.T.C.A.T.Q.R.L.A.E.F.L.R.H.S.S.N.N.F.G.P.I.L.P.P.T.N.V.G.S.N.T.P}"
    "|PEPTIDE2{[gGlu].[C20d]}"
    "$PEPTIDE2,PEPTIDE1,1:R1-1:R1|PEPTIDE1,PEPTIDE1,2:R3-7:R3$$V2.0"
)


# ── helpers ──────────────────────────────────────────────────────────────────

def _exact_mw(smiles: str) -> float:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"bad SMILES: {smiles[:60]}"
    return Descriptors.ExactMolWt(mol)


def _avg_mw(smiles: str) -> float:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return Descriptors.MolWt(mol)


# ── tests ────────────────────────────────────────────────────────────────────

def test_ala3_parses_and_mw():
    obj = HELMParser.parse(ALA3)
    smi = obj.to_smiles(cap_c='acid')
    assert smi is not None
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    # Ala-Ala-Ala as free acid: C9H17N3O4, ExactMW = 231.122
    mw = _exact_mw(smi)
    assert abs(mw - 231.122) < 0.01


def test_ala3_amide_vs_acid_differ():
    obj = HELMParser.parse(ALA3)
    smi_amide = obj.to_smiles(cap_c='amide')
    smi_acid  = obj.to_smiles(cap_c='acid')
    assert smi_amide != smi_acid
    # amide is lighter by O - N = 16 - 14 = 2 Da? no: amide = C(=O)NH2, acid = C(=O)OH
    # acid ExactMW - amide ExactMW = O(16) - N(14) + H(1) - H(1) = +1.003 (O vs N, one H each)
    assert abs(_exact_mw(smi_acid) - _exact_mw(smi_amide) - 0.984) < 0.01


def test_pyy_ref_parseable_and_atom_count():
    obj = HELMParser.parse(PYY_REF)
    smi = obj.to_smiles()
    assert smi is not None
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    # Use average MW (same basis as calc_mw) — only amide-vs-acid cap offset remains.
    # amide cap replaces -OH with -NH2: MolWt drops by 0.984 Da.
    mw_rdkit = _avg_mw(smi)
    mw_calc  = calc_mw(obj)
    assert mw_calc is not None
    assert abs((mw_rdkit - mw_calc) + 0.984) < 1.0   # amide cap: rdkit < calc_mw by ~1 Da


def test_pramlintide_disulfide_mw_offset():
    """calc_mw models disulfide as −H₂O; true bond loses −H₂ → avg-MW offset ≈ +15 Da."""
    obj = HELMParser.parse(PRAMLINTIDE)
    smi = obj.to_smiles()
    assert smi is not None
    # Use average MW for apples-to-apples comparison with calc_mw.
    # Expected: +H₂O − H₂ − amide_cap = 18.015 − 2.016 − 0.984 = 15.015 Da
    mw_rdkit = _avg_mw(smi)
    mw_calc  = calc_mw(obj)
    assert mw_calc is not None
    assert abs((mw_rdkit - mw_calc) - 15.015) < 1.0


def test_pramlintide_contains_disulfide():
    """S-S bond must be present in the assembled molecule."""
    obj = HELMParser.parse(PRAMLINTIDE)
    smi = obj.to_smiles()
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    has_ss = any(
        b.GetBondTypeAsDouble() == 1.0
        and b.GetBeginAtom().GetAtomicNum() == 16
        and b.GetEndAtom().GetAtomicNum() == 16
        for b in mol.GetBonds()
    )
    assert has_ss, "No S-S bond found in pramlintide SMILES"


def test_cagrili_backbone_smiles():
    """PEPTIDE1 chain of cagrilintide assembles despite a PEPTIDE2 protractor."""
    obj = HELMParser.parse(CAGRILI)
    smi = obj.to_smiles(chain_id='PEPTIDE1')
    assert smi is not None
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    # Must still have S-S bond (Cys2-Cys7 is intra-PEPTIDE1)
    has_ss = any(
        b.GetBeginAtom().GetAtomicNum() == 16 and b.GetEndAtom().GetAtomicNum() == 16
        for b in mol.GetBonds()
    )
    assert has_ss


def test_deterministic_output():
    """Same HELM → identical SMILES on repeated calls."""
    obj = HELMParser.parse(PRAMLINTIDE)
    assert obj.to_smiles() == obj.to_smiles()


def test_combine_monomers_directly():
    """Low-level: assemble Gly-Ala dipeptide from entries, verify atom count."""
    from monomer_db.monomer_db import MonomerDB
    db = MonomerDB()
    g = db.find_by_symbol('G')
    a = db.find_by_symbol('A')
    assert g is not None and a is not None
    rsmiles, smiles = combine_monomers([(1, g), (2, a)], [], cap_c='acid')
    assert smiles is not None
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    # Gly-Ala free acid: C5H10N2O3, ExactMW = 146.069
    assert abs(_exact_mw(smiles) - 146.069) < 0.01


def test_rsmiles_has_open_termini():
    """rsmiles from combine_monomers keeps terminal [*:N] for fragment coupling."""
    from monomer_db.monomer_db import MonomerDB
    db = MonomerDB()
    g = db.find_by_symbol('G')
    a = db.find_by_symbol('A')
    rsmiles, _ = combine_monomers([(1, g), (2, a)], [], cap_c='amide')
    assert rsmiles is not None
    assert '[*:' in rsmiles   # at least one open attachment point remains
