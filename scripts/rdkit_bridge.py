"""
RDKit bridge for peptide informatics.

Public API
----------
monomer_to_mol(entry)                        → Chem.Mol   (R-groups capped, ready for descriptors)
mol_to_descriptors(mol)                      → dict        (MW, LogP, TPSA, HBD, HBA, RotBonds, Rings, QED)
validate_roundtrip(smiles, mol)              → dict        (tanimoto, substructure_match, pass)
combine_monomers(pos_entry_pairs, conns, …)  → (rsmiles, smiles)  assemble chain SMILES
helm_obj_to_smiles(helm_obj, …)              → str | None  full HELM → SMILES
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


# ---------------------------------------------------------------------------
# HELM → SMILES  (combine_monomers + helm_obj_to_smiles)
# ---------------------------------------------------------------------------
# Global map-number convention: position p (effective) → R1=10p+1, R2=10p+2, R3=10p+3.
# All capping and remapping uses the RDKit atom API (SetAtomicNum / SetAtomMapNum)
# so no SMILES strings are manipulated after the initial mol parse.

# Pistoia single-atom cap strings → RDKit atomic numbers.
# Covers all PEPTIDE backbone caps in HELMCoreLibrary and custom_monomers.
_CAP_ATOMICNUM: dict[str, int] = {'[H]': 1, 'H': 1, 'O': 8, 'N': 7, 'S': 16}


def _remap_mol_mapnums(mol, eff_pos: int):
    """Return a new Mol with [*:N] map numbers remapped to 10*eff_pos + N."""
    from rdkit.Chem import RWMol
    rw = RWMol(mol)
    for atom in rw.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn > 0:
            atom.SetAtomMapNum(10 * eff_pos + mn)
    return rw.GetMol()


def _has_mapnum(mol, mapnum: int) -> bool:
    return any(a.GetAtomMapNum() == mapnum for a in mol.GetAtoms())


def _form_bond(rwmol, from_mapnum: int, to_mapnum: int) -> None:
    """
    Bond the real neighbors of the two [*:N] dummies, then remove both dummies.

    Re-finds atoms by map number on each call so indices stay valid across
    prior RemoveAtom operations.
    """
    from rdkit.Chem import rdchem

    fd = next((a for a in rwmol.GetAtoms() if a.GetAtomMapNum() == from_mapnum), None)
    td = next((a for a in rwmol.GetAtoms() if a.GetAtomMapNum() == to_mapnum),   None)
    if fd is None or td is None:
        logger.warning("_form_bond: map nums %d/%d not found", from_mapnum, to_mapnum)
        return
    if not fd.GetNeighbors() or not td.GetNeighbors():
        logger.warning("_form_bond: dummy %d or %d has no neighbor", from_mapnum, to_mapnum)
        return

    rwmol.AddBond(fd.GetNeighbors()[0].GetIdx(),
                  td.GetNeighbors()[0].GetIdx(),
                  rdchem.BondType.SINGLE)

    for idx in sorted([fd.GetIdx(), td.GetIdx()], reverse=True):
        rwmol.RemoveAtom(idx)


def _cap_dummy(rwmol, mapnum: int, atomic_num: int) -> None:
    """Replace [*:mapnum] dummy in-place with a real atom of atomic_num.

    SetNoImplicit(False) lets RDKit fill the remaining valence with implicit H
    during the next SanitizeMol (e.g. C-terminal N becomes NH₂, not bare N).
    """
    atom = next((a for a in rwmol.GetAtoms() if a.GetAtomMapNum() == mapnum), None)
    if atom is None:
        return
    atom.SetAtomicNum(atomic_num)
    atom.SetAtomMapNum(0)
    atom.SetNoImplicit(False)


def combine_monomers(
    pos_entry_pairs: list[tuple[int, dict]],
    intra_connections: list[dict],
    cap_c: str = 'amide',
    chain_offset: int = 0,
) -> tuple[str | None, str | None]:
    """
    Assemble a monomer sequence into SMILES via explicit bond formation.

    All operations use the RDKit atom graph (SetAtomicNum / SetAtomMapNum /
    AddBond / RemoveAtom).  No SMILES strings are manipulated after the initial
    Chem.MolFromSmiles call.

    Parameters
    ----------
    pos_entry_pairs : [(pos, entry), …] sorted ascending by HELM position.
    intra_connections : connectivity_graph entries filtered to this chain/fragment.
        Keys: from_pos, from_rgroup ('R1'/'R2'/'R3'), to_pos, to_rgroup.
    cap_c : 'amide' (−CONH₂) | 'acid' (−COOH)
    chain_offset : shift all positions by this amount to avoid map-number
        collisions when assembling multi-chain molecules.

    Returns
    -------
    (rsmiles, smiles)
        rsmiles : Pistoia-format SMILES ([*:1] N-term, [*:2] C-term open)
                  ready to be stored as a custom monomer entry.
        smiles  : fully capped, RDKit-canonical SMILES.
        Either field is None on failure (missing DB entry, parse error, …).
    """
    from rdkit import Chem
    from rdkit.Chem import RWMol
    from monomer_db.monomer_db import _helm_to_star_smiles

    if not pos_entry_pairs:
        return None, None

    # ── Phase 1: parse each monomer, remap to global map numbers ──────────────
    mapnum_to_atomicnum: dict[int, int] = {}   # global_mapnum → atomic number for capping
    mols: list[Chem.Mol] = []

    for pos, entry in pos_entry_pairs:
        if entry is None:
            logger.warning("combine_monomers: no DB entry at pos %d", pos)
            return None, None
        raw = entry.get('smiles', '')
        if not raw:
            logger.warning("combine_monomers: empty SMILES for %r", entry.get('symbol'))
            return None, None

        mol = Chem.MolFromSmiles(_helm_to_star_smiles(raw))
        if mol is None:
            logger.warning("combine_monomers: cannot parse SMILES for %r", entry.get('symbol'))
            return None, None

        eff_pos = pos + chain_offset
        mol = _remap_mol_mapnums(mol, eff_pos)

        for rg in entry.get('rgroups', []):
            m = re.match(r'R(\d+)', rg.get('label', ''))
            if not m:
                continue
            cap_str = _cap_from_capgroup(rg.get('capGroupSmiles', '') or '[*:1][H]')
            mapnum_to_atomicnum[10 * eff_pos + int(m.group(1))] = (
                _CAP_ATOMICNUM.get(cap_str, 1)   # default to H for unknown cap strings
            )
        mols.append(mol)

    # ── Phase 2: combine + form backbone + explicit bonds ─────────────────────
    combined = mols[0]
    for m in mols[1:]:
        combined = Chem.CombineMols(combined, m)
    rwmol = RWMol(combined)

    eff_positions = [pos + chain_offset for pos, _ in pos_entry_pairs]
    for i in range(len(eff_positions) - 1):
        _form_bond(rwmol, 10 * eff_positions[i] + 2, 10 * eff_positions[i + 1] + 1)

    for conn in intra_connections:
        fp = conn['from_pos'] + chain_offset
        tp = conn['to_pos']   + chain_offset
        fm_re = re.search(r'\d+', conn['from_rgroup'])
        tm_re = re.search(r'\d+', conn['to_rgroup'])
        if not fm_re or not tm_re:
            continue
        fm, tm = 10 * fp + int(fm_re.group()), 10 * tp + int(tm_re.group())
        if _has_mapnum(rwmol, fm) and _has_mapnum(rwmol, tm):
            _form_bond(rwmol, fm, tm)

    # ── Phase 3a: rsmiles — cap non-terminal dummies; normalize terminals ──────
    # Two copies: rsmiles keeps terminals open (Pistoia [*:1]/[*:2]);
    #             smiles caps everything.
    first_r1 = 10 * eff_positions[0]  + 1
    last_r2  = 10 * eff_positions[-1] + 2

    rw_r = RWMol(rwmol)
    for atom in rw_r.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn > 0 and mn not in (first_r1, last_r2):
            _cap_dummy(rw_r, mn, mapnum_to_atomicnum.get(mn, 1))
    # Normalize terminal map numbers → standard Pistoia [*:1] / [*:2]
    for atom in rw_r.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn == first_r1:
            atom.SetAtomMapNum(1)
        elif mn == last_r2:
            atom.SetAtomMapNum(2)
    try:
        Chem.SanitizeMol(rw_r)
        rsmiles = Chem.MolToSmiles(rw_r)
    except Exception as e:
        logger.warning("combine_monomers: rsmiles sanitize failed: %s", e)
        rsmiles = None

    # ── Phase 3b: smiles — cap all remaining dummies including termini ─────────
    c_term_atomicnum = 7 if cap_c == 'amide' else 8   # N=amide, O=acid
    rw_s = RWMol(rwmol)
    _cap_dummy(rw_s, first_r1, 1)              # N-term: free amine (H)
    _cap_dummy(rw_s, last_r2, c_term_atomicnum)
    for atom in rw_s.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn > 0:
            _cap_dummy(rw_s, mn, mapnum_to_atomicnum.get(mn, 1))
    try:
        Chem.SanitizeMol(rw_s)
        smiles = Chem.MolToSmiles(rw_s)
    except Exception as e:
        logger.warning("combine_monomers: smiles sanitize failed: %s", e)
        smiles = None

    return rsmiles, smiles


def helm_obj_to_smiles(
    helm_obj,
    chain_id: str | None = None,
    cap_c: str = 'amide',
) -> str | None:
    """
    Convert a HELMObject chain to a fully capped canonical SMILES string.

    Handles intra-chain connections (disulfide, lactam, etc.).
    Inter-chain connections (protractor/lipid) are not yet assembled — only
    the specified chain is returned.

    Returns None if any monomer in the chain is missing from the DB.
    """
    chain = helm_obj.get_chain(chain_id)
    if chain is None:
        return None

    cid = chain['chain_id']
    pos_entry_pairs = [(m['pos'], m['entry']) for m in chain['monomers']]

    intra = [
        c for c in helm_obj.data.get('connectivity_graph', [])
        if c['from_chain'] == cid and c['to_chain'] == cid
    ]

    _, smiles = combine_monomers(pos_entry_pairs, intra, cap_c=cap_c)
    return smiles
