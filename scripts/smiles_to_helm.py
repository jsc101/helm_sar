"""SMILES → HELM fragmentation and assembly.

Algorithm
---------
1. Find all backbone amide N-C(=O) bonds via SMARTS.
2. Fragment at those bonds with RDKit FragmentOnBonds.
3. For each fragment, identify which dummy atoms are R1 (N-terminal) vs
   R2 (C-terminal) by examining their neighbor's atom type.
4. Canonicalise each fragment and look it up in MonomerDB.
5. Unknown fragments are logged to a JSON file and returned as new_monomers.
6. Reassemble the ordered sequence into a HELM string.

Limitations
-----------
- Cyclic peptides: backbone ring is detected by checking if FragmentOnBonds
  produced ring-closure dummies; works for simple head-to-tail cyclic.
- N-methyl backbone amides are treated as backbone bonds.
- Sidechain amides (Asn, Gln, Arg) are excluded via neighbor valence check.
"""
from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path
from typing import NamedTuple

from rdkit import Chem

from monomer_db.monomer_db import MonomerDB

# Backbone amide N-C(=O):
#   D2 = NH secondary (two heavy-atom neighbours: amide-C + alpha-C)
#   D3 = N-methyl tertiary (amide-C + alpha-C + methyl)
# Excludes primary sidechain amides (Asn/Gln NH2, D1) and ureas/guanidines
# (which lack the C=O oxygen requirement).
_AMIDE_SMARTS = Chem.MolFromSmarts('[N;D2,D3:1]-[C;X3:2](=[O;X1])')


class CrossBond(NamedTuple):
    """Asymmetric inter-fragment connection (e.g. PEPTIDE↔CHEM for click chemistry)."""
    partner_idx: int      # index into the fragment list
    my_rgroup: str        # e.g. 'R3'
    partner_rgroup: str   # e.g. 'R1'


class Fragment(NamedTuple):
    smiles: str                        # canonical SMILES with [*:1]/[*:2]/[*:3] attachment points
    symbol: str | None                 # MonomerDB symbol, or None if unknown
    entry: dict | None                 # full MonomerDB entry, or None
    r3_partners: tuple[int, ...] = ()  # symmetric R3–R3 indices (disulfide)
    r3_type: str | None = None         # 'disulfide', 'click', etc.
    cross_bonds: tuple[CrossBond, ...] = ()  # asymmetric connections (click chem PEPTIDE↔CHEM)


def _find_backbone_amide_bonds(mol: Chem.Mol) -> list[int]:
    """Return bond indices for backbone amide N-C(=O) bonds.

    Uses N degree to exclude primary sidechain amides (Asn/Gln NH2, degree 1
    in terms of heavy-atom connections).  N-methyl backbone amides (degree 3)
    are included.
    """
    if _AMIDE_SMARTS is None:
        return []

    bond_ids: list[int] = []
    seen: set[int] = set()

    for match in mol.GetSubstructMatches(_AMIDE_SMARTS):
        n_idx, c_idx = match[0], match[1]
        bond = mol.GetBondBetweenAtoms(n_idx, c_idx)
        if bond is None or bond.GetIdx() in seen:
            continue
        seen.add(bond.GetIdx())
        bond_ids.append(bond.GetIdx())

    return bond_ids


def _dummy_role(mol: Chem.Mol, dummy_idx: int) -> str:
    """Return 'R1' if dummy connects to N-side, 'R2' if to C=O side."""
    dummy = mol.GetAtomWithIdx(dummy_idx)
    for nb in dummy.GetNeighbors():
        if nb.GetAtomicNum() == 7:
            return 'R1'   # dummy → N → backbone N-terminus side
        if nb.GetAtomicNum() == 6:
            for nb2 in nb.GetNeighbors():
                b = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                if b and nb2.GetAtomicNum() == 8 and b.GetBondTypeAsDouble() == 2.0:
                    return 'R2'   # dummy → C=O → backbone C-terminus side
    return 'R1'   # fallback


def _assign_rgroups(frag_mol: Chem.Mol) -> str | None:
    """Convert a fragment with [N*] dummies (FragmentOnBonds isotope labels)
    to clean [*:1]/[*:2] attachment-point SMILES."""
    rw = Chem.RWMol(frag_mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            role = _dummy_role(frag_mol, atom.GetIdx())
            atom.SetIsotope(0)              # clear FragmentOnBonds isotope tag
            atom.SetAtomMapNum(1 if role == 'R1' else 2)
    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw)
    except Exception:
        return None


# Free primary amine or thiol at end of sidechain → promote to [*:3] for DB lookup.
# Covers Lys (CCCCN), Orn (CCCN), Dab (CCN), Dap (CN), Cys (CS).
_R3_PATTERNS: list[tuple] = [
    (Chem.MolFromSmarts('[N;H2;D1]'), 7),    # primary amine NH2
    (Chem.MolFromSmarts('[S;H1;D1]'), 16),   # free thiol SH
]

# Sidechain COOH pattern — Pistoia DB stores the OH as [OH:3] → [*:3], meaning
# the OH is REPLACED by the dummy.  Covers Asp (D), Glu (E), gGlu.
_COOH_SMARTS = Chem.MolFromSmarts('[C;X3](=[O;X1])[O;H1;D1]')

# 1,4-disubstituted 1,2,3-triazole (CuAAC product):
#   N1 (X3 — bonded to Aha sidechain + 2 ring atoms) → R1 of Triazole14
#   C4 (X3 — bonded to Hpg sidechain + 2 ring atoms) → R2 of Triazole14
#   N2, N3 (X2 — ring bonds only), C5 (H1 — one H, two ring bonds)
_TRIAZOLE_SMARTS = Chem.MolFromSmarts('[n;r5;X3:1]1[n;r5;X2][n;r5;X2][c;r5;X3:4][c;r5;H1]1')


def _try_add_r3(frag_smiles: str) -> str | None:
    """Return a variant of frag_smiles with [*:3] added as a neighbour of any
    free sidechain primary amine or thiol (mirroring the DB format N[*:3]).

    Does NOT replace the heteroatom — it adds a dummy atom bonded to it, which
    matches the Pistoia convention where R-groups are departing-H placeholders.
    """
    mol = Chem.MolFromSmiles(frag_smiles)
    if mol is None:
        return None
    for pat, _ in _R3_PATTERNS:
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        if not matches:
            continue
        for (idx,) in matches:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomMapNum() != 0:
                continue
            rw = Chem.RWMol(Chem.RWMol(mol))
            dummy_idx = rw.AddAtom(Chem.Atom(0))
            rw.GetAtomWithIdx(dummy_idx).SetAtomMapNum(3)
            rw.GetAtomWithIdx(dummy_idx).SetNoImplicit(True)
            rw.AddBond(idx, dummy_idx, Chem.rdchem.BondType.SINGLE)
            try:
                Chem.SanitizeMol(rw)
                return Chem.MolToSmiles(rw)
            except Exception:
                pass
    return None


def _try_promote_cooh_r3(frag_smiles: str) -> str | None:
    """Replace a free sidechain carboxylic acid OH with [*:3].

    Asp/Glu store their sidechain COOH as [OH:3] in the Pistoia DB, which
    canonicalises to C(=O)[*:3].  This transform converts the free -COOH
    in a fragment to that form so the DB lookup can match.
    """
    if _COOH_SMARTS is None:
        return None
    mol = Chem.MolFromSmiles(frag_smiles)
    if mol is None:
        return None
    matches = mol.GetSubstructMatches(_COOH_SMARTS)
    if not matches:
        return None
    for match in matches:
        _c_idx, _o_dbl_idx, oh_idx = match
        oh_atom = mol.GetAtomWithIdx(oh_idx)
        if oh_atom.GetAtomMapNum() != 0:
            continue  # already mapped
        rw = Chem.RWMol(mol)
        oh = rw.GetAtomWithIdx(oh_idx)
        oh.SetAtomicNum(0)
        oh.SetAtomMapNum(3)
        oh.SetNoImplicit(True)
        oh.SetNumExplicitHs(0)
        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw)
        except Exception:
            pass
    return None


def _try_cap_terminals(frag_smiles: str) -> str | None:
    """For terminal residue fragments, add the missing backbone R-group.

    N-terminal (has [*:2], no [*:1]): add [*:1] bonded to free NH2.
    C-terminal (has [*:1], no [*:2]): replace free -COOH or -C(=O)NH2 with [*:2].

    This allows terminal residues to match DB entries that carry both R1 and R2.
    """
    mol = Chem.MolFromSmiles(frag_smiles)
    if mol is None:
        return None

    map_nums = {a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum() > 0}
    has_r1, has_r2 = 1 in map_nums, 2 in map_nums
    if has_r1 == has_r2:
        return None  # both present or both absent — nothing to do

    rw = Chem.RWMol(mol)

    if has_r2 and not has_r1:
        # N-terminal: find free primary amine (unmapped, ≥2 Hs).
        # Prefer the alpha-amine (closest to [*:2]) to avoid mis-capping Lys
        # epsilon-NH2 when two primary amines are present.
        r2_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == 2)
        candidates = [
            a for a in mol.GetAtoms()
            if a.GetAtomicNum() == 7 and a.GetAtomMapNum() == 0 and a.GetTotalNumHs() >= 2
        ]
        candidates.sort(key=lambda a: len(Chem.GetShortestPath(mol, a.GetIdx(), r2_idx)))
        for atom in candidates:
            rw = Chem.RWMol(mol)
            d_idx = rw.AddAtom(Chem.Atom(0))
            rw.GetAtomWithIdx(d_idx).SetAtomMapNum(1)
            rw.GetAtomWithIdx(d_idx).SetNoImplicit(True)
            rw.AddBond(atom.GetIdx(), d_idx, Chem.rdchem.BondType.SINGLE)
            try:
                Chem.SanitizeMol(rw)
                return Chem.MolToSmiles(rw)
            except Exception:
                pass

    if has_r1 and not has_r2:
        # C-terminal: replace free -COOH or primary C-amide with [*:2]
        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() != 6 or atom.GetAtomMapNum() != 0:
                continue
            has_dbl_o = any(
                nb.GetAtomicNum() == 8
                and rw.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                for nb in atom.GetNeighbors()
            )
            if not has_dbl_o:
                continue
            for nb in list(atom.GetNeighbors()):
                b = rw.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if b.GetBondTypeAsDouble() != 1.0:
                    continue
                if nb.GetAtomicNum() in (7, 8) and nb.GetAtomMapNum() == 0 and nb.GetTotalNumHs() >= 1:
                    nb.SetAtomicNum(0)
                    nb.SetAtomMapNum(2)
                    nb.SetNoImplicit(True)
                    nb.SetNumExplicitHs(0)
                    try:
                        Chem.SanitizeMol(rw)
                        return Chem.MolToSmiles(rw)
                    except Exception:
                        rw = Chem.RWMol(mol)  # reset on failure

    return None


def _lookup_fragment(frag_smiles: str, db: MonomerDB) -> tuple[str | None, dict | None]:
    """Try to find monomer using a cascade of SMILES transforms.

    Order: direct stereo → direct nostereo → R3 amine/thiol → R3 sidechain COOH
           → terminal cap → terminal cap + R3 amine → terminal cap + COOH R3
    """
    def _find(smi: str) -> dict | None:
        return db.find(smi, stereo=True) or db.find(smi, stereo=False)

    entry = _find(frag_smiles)
    if entry:
        return entry['symbol'], entry

    r3_smi = _try_add_r3(frag_smiles)
    if r3_smi:
        entry = _find(r3_smi)
        if entry:
            return entry['symbol'], entry

    cooh_smi = _try_promote_cooh_r3(frag_smiles)
    if cooh_smi:
        entry = _find(cooh_smi)
        if entry:
            return entry['symbol'], entry

    term_smi = _try_cap_terminals(frag_smiles)
    if term_smi:
        entry = _find(term_smi)
        if entry:
            return entry['symbol'], entry

        r3_term = _try_add_r3(term_smi)
        if r3_term:
            entry = _find(r3_term)
            if entry:
                return entry['symbol'], entry

        cooh_term = _try_promote_cooh_r3(term_smi)
        if cooh_term:
            entry = _find(cooh_term)
            if entry:
                return entry['symbol'], entry

    return None, None


def _helm_symbol(symbol: str) -> str:
    return symbol if len(symbol) == 1 and symbol.isalpha() else f'[{symbol}]'


def _unk_symbol(canonical_smiles: str) -> str:
    """Return a deterministic 6-char hex symbol for an unknown monomer fragment.

    The hash is derived from the canonical SMILES, so the same structure always
    gets the same symbol across sessions and runs — no counter drift.
    """
    digest = hashlib.sha256(canonical_smiles.encode()).hexdigest()[:6]
    return f'UNK_{digest}'


def _split_disulfide(
    frag_smiles: str, db: MonomerDB
) -> list[tuple[str, str | None, dict | None]] | None:
    """If a backbone fragment contains an S-S bond, cut it and return sub-fragments.

    Each resulting Cys-like fragment gets [*:3] placed on the now-free sulfur,
    matching the Pistoia DB format ``CS[*:3]``.  Returns None if no S-S bond found.
    """
    mol = Chem.MolFromSmiles(frag_smiles)
    if mol is None:
        return None
    ss_bonds = [
        b.GetIdx() for b in mol.GetBonds()
        if b.GetBeginAtom().GetAtomicNum() == 16 and b.GetEndAtom().GetAtomicNum() == 16
    ]
    if not ss_bonds:
        return None

    cut_mol = Chem.FragmentOnBonds(mol, ss_bonds, addDummies=True)
    sub_frags = Chem.GetMolFrags(cut_mol, asMols=True, sanitizeFrags=False)

    results: list[tuple[str, str | None, dict | None]] = []
    for sf in sub_frags:
        rw = Chem.RWMol(sf)
        # Newly created dummies (no map num yet) bonded to S → assign R3
        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == 0:
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 16:
                        atom.SetIsotope(0)
                        atom.SetAtomMapNum(3)
                        break
        try:
            Chem.SanitizeMol(rw)
            fsmi = Chem.MolToSmiles(rw)
            symbol, entry = _lookup_fragment(fsmi, db)
            results.append((fsmi, symbol, entry))
        except Exception:
            pass

    return results if len(results) >= 2 else None


def _split_triazole(
    frag_smiles: str, db: MonomerDB
) -> list[tuple[str, str | None, dict | None]] | None:
    """If a backbone fragment contains a 1,4-triazole, split into three pieces.

    Cuts the two exocyclic bonds of the triazole ring:
      - N1–C(Aha sidechain): assigns iso-1 dummy to R1 on the triazole side, R3 on Aha side
      - C4–C(Hpg sidechain): assigns iso-2 dummy to R2 on the triazole side, R3 on Hpg side

    Returns a list of exactly 3 (smiles, symbol, entry) tuples in order:
        [Aha-piece, Triazole-piece, Hpg-piece]
    or None if no 1,4-triazole is detected.

    Note: correct backbone ordering of the returned triplet assumes Aha appears
    before Hpg in the original SMILES (standard N→C notation).
    """
    if _TRIAZOLE_SMARTS is None:
        return None
    mol = Chem.MolFromSmiles(frag_smiles)
    if mol is None:
        return None
    matches = mol.GetSubstructMatches(_TRIAZOLE_SMARTS)
    if not matches:
        return None

    # match = (n1_idx, n2_idx, n3_idx, c4_idx, c5_idx)
    match = matches[0]
    n1_idx, _n2, _n3, c4_idx, _c5 = match
    ring_atoms = set(match)

    # Find the exocyclic C bonded to N1 (Aha's CH2) and C4 (Hpg's CH2)
    n1_ext = next(
        (nb.GetIdx() for nb in mol.GetAtomWithIdx(n1_idx).GetNeighbors()
         if nb.GetIdx() not in ring_atoms and nb.GetAtomicNum() == 6),
        None
    )
    c4_ext = next(
        (nb.GetIdx() for nb in mol.GetAtomWithIdx(c4_idx).GetNeighbors()
         if nb.GetIdx() not in ring_atoms and nb.GetAtomicNum() == 6),
        None
    )
    if n1_ext is None or c4_ext is None:
        return None

    bond_n1 = mol.GetBondBetweenAtoms(n1_idx, n1_ext).GetIdx()
    bond_c4 = mol.GetBondBetweenAtoms(c4_idx, c4_ext).GetIdx()

    # Cut both bonds: bond_n1 → isotope 1 dummies, bond_c4 → isotope 2 dummies
    cut_mol = Chem.FragmentOnBonds(mol, [bond_n1, bond_c4], addDummies=True)
    sub_frags = Chem.GetMolFrags(cut_mol, asMols=True, sanitizeFrags=False)
    if len(sub_frags) != 3:
        return None

    # FragmentOnBonds isotope convention: when atom A-B is cut, the dummy on
    # A's side has isotope = B's atom index; the dummy on B's side has isotope
    # = A's atom index.  Track the four expected isotope labels:
    #   triazole N1 side: loses n1_ext → dummy iso = n1_ext
    #   Aha side:         loses n1_idx → dummy iso = n1_idx
    #   triazole C4 side: loses c4_ext → dummy iso = c4_ext
    #   Hpg side:         loses c4_idx → dummy iso = c4_idx
    iso_tri_n1 = n1_ext   # on triazole side of N1-cut: bonded to N1, iso=n1_ext
    iso_tri_c4 = c4_ext   # on triazole side of C4-cut: bonded to C4, iso=c4_ext
    iso_aha    = n1_idx   # on Aha side of N1-cut:    bonded to C_ext, iso=n1_idx
    iso_hpg    = c4_idx   # on Hpg side of C4-cut:   bonded to C_ext2, iso=c4_idx

    aha_frag = triazole_frag = hpg_frag = None

    for sf in sub_frags:
        iso_set = {a.GetIsotope() for a in sf.GetAtoms()
                   if a.GetAtomicNum() == 0 and a.GetIsotope() != 0}
        rw = Chem.RWMol(sf)

        if iso_tri_n1 in iso_set and iso_tri_c4 in iso_set:
            # Triazole piece: N1-side dummy → R1, C4-side dummy → R2
            for atom in rw.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    iso = atom.GetIsotope()
                    atom.SetIsotope(0)
                    if iso == iso_tri_n1:
                        atom.SetAtomMapNum(1)  # R1 on N1
                    elif iso == iso_tri_c4:
                        atom.SetAtomMapNum(2)  # R2 on C4
            try:
                Chem.SanitizeMol(rw)
                fsmi = Chem.MolToSmiles(rw)
                sym, entry = _lookup_fragment(fsmi, db)
                triazole_frag = (fsmi, sym, entry)
            except Exception:
                return None

        elif iso_aha in iso_set:
            # Aha piece: single dummy → R3
            for atom in rw.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetIsotope() == iso_aha:
                    atom.SetIsotope(0)
                    atom.SetAtomMapNum(3)
            try:
                Chem.SanitizeMol(rw)
                fsmi = Chem.MolToSmiles(rw)
                sym, entry = _lookup_fragment(fsmi, db)
                aha_frag = (fsmi, sym, entry)
            except Exception:
                return None

        elif iso_hpg in iso_set:
            # Hpg piece: single dummy → R3
            for atom in rw.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetIsotope() == iso_hpg:
                    atom.SetIsotope(0)
                    atom.SetAtomMapNum(3)
            try:
                Chem.SanitizeMol(rw)
                fsmi = Chem.MolToSmiles(rw)
                sym, entry = _lookup_fragment(fsmi, db)
                hpg_frag = (fsmi, sym, entry)
            except Exception:
                return None

    if aha_frag is None or triazole_frag is None or hpg_frag is None:
        return None

    return [aha_frag, triazole_frag, hpg_frag]


def fragment_smiles(smiles: str, db: MonomerDB | None = None) -> list[Fragment]:
    """Decompose a peptide SMILES into ordered Fragment objects.

    For linear peptides the fragments are ordered N→C.
    For cyclic peptides (ring-backbone) the start position is arbitrary
    (first amide bond cut site), but the sequence is consistent.
    """
    if db is None:
        db = MonomerDB()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    bond_ids = _find_backbone_amide_bonds(mol)
    if not bond_ids:
        return []

    # FragmentOnBonds returns a mol with dummy atoms numbered per cut site
    frag_mol = Chem.FragmentOnBonds(mol, bond_ids, addDummies=True)
    frag_mols = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)

    fragments: list[Fragment] = []
    for fm in frag_mols:
        fsmi = _assign_rgroups(fm)
        if fsmi is None:
            continue

        # Layer 2a: S-S bond → split into two Cys-like monomers (symmetric R3–R3).
        ss_sub = _split_disulfide(fsmi, db)
        if ss_sub:
            base = len(fragments)
            n = len(ss_sub)
            for j, (sub_smi, sub_sym, sub_entry) in enumerate(ss_sub):
                partners = tuple(base + k for k in range(n) if k != j)
                fragments.append(Fragment(
                    smiles=sub_smi, symbol=sub_sym, entry=sub_entry,
                    r3_partners=partners, r3_type='disulfide',
                ))
            continue

        # Layer 2b: 1,4-triazole → split into (Aha-piece, Triazole14, Hpg-piece).
        # CrossBonds encode the asymmetric PEPTIDE↔CHEM connections:
        #   Aha:R3 → Triazole14:R1,  Triazole14:R2 → Hpg:R3
        triazole_sub = _split_triazole(fsmi, db)
        if triazole_sub:
            base = len(fragments)
            aha_t, tri_t, hpg_t = triazole_sub
            fragments.append(Fragment(
                smiles=aha_t[0], symbol=aha_t[1], entry=aha_t[2],
                cross_bonds=(CrossBond(base + 1, 'R3', 'R1'),),
            ))
            fragments.append(Fragment(
                smiles=tri_t[0], symbol=tri_t[1], entry=tri_t[2],
                cross_bonds=(CrossBond(base, 'R1', 'R3'), CrossBond(base + 2, 'R2', 'R3')),
            ))
            fragments.append(Fragment(
                smiles=hpg_t[0], symbol=hpg_t[1], entry=hpg_t[2],
                cross_bonds=(CrossBond(base + 1, 'R3', 'R2'),),
            ))
            continue

        symbol, entry = _lookup_fragment(fsmi, db)
        fragments.append(Fragment(smiles=fsmi, symbol=symbol, entry=entry))

    return fragments


def smiles_to_helm(
    smiles: str,
    chain_id: str = 'PEPTIDE1',
    cyclic: bool = False,
    db: MonomerDB | None = None,
    new_monomer_log: str | Path | None = None,
) -> tuple[str | None, list[dict]]:
    """Convert a peptide SMILES to a HELM string.

    Unknown fragments are auto-assigned content-addressed symbols like
    ``UNK_a3f9c2`` derived from a SHA256 hash of their canonical SMILES.
    The same structure always gets the same symbol across sessions — no
    counter drift between runs.  Symbols are registered live in *db* so
    repeated occurrences within a run resolve immediately.

    Parameters
    ----------
    smiles : str
        Canonical SMILES of the peptide (linear or head-to-tail cyclic).
    chain_id : str
        HELM chain identifier (default 'PEPTIDE1').
    cyclic : bool
        True for head-to-tail backbone-cyclic peptides.
    db : MonomerDB | None
        Monomer database.  A fresh instance is created if None; pass the
        same instance across calls to accumulate discoveries.
    new_monomer_log : path | None
        JSON file (dict keyed by symbol) to persist pending monomers for
        chemist review.

    Returns
    -------
    (helm_str, new_monomers)
        helm_str     — complete HELM string (None only if no backbone bonds).
        new_monomers — list of entry dicts for newly registered monomers;
                       each has 'symbol', 'smiles', 'name' (None),
                       'assigned_symbol' (None), 'needs_review' (True).
    """
    if db is None:
        db = MonomerDB()

    fragments = fragment_smiles(smiles, db)
    if not fragments:
        return None, []

    new_monomers: list[dict] = []
    resolved: list[Fragment] = []

    for frag in fragments:
        if frag.symbol is not None:
            resolved.append(frag)
            continue
        # Same SMILES already registered this session or in DB
        existing = db.find(frag.smiles, stereo=True) or db.find(frag.smiles, stereo=False)
        if existing:
            resolved.append(Fragment(
                smiles=frag.smiles, symbol=existing['symbol'], entry=existing,
                r3_partners=frag.r3_partners, r3_type=frag.r3_type,
                cross_bonds=frag.cross_bonds,
            ))
            continue
        sym = _unk_symbol(frag.smiles)
        entry = {
            'symbol':          sym,
            'name':            None,
            'polymerType':     'PEPTIDE',
            'monomerType':     'Backbone',
            'smiles':          frag.smiles,
            'assigned_symbol': None,
            'needs_review':    True,
        }
        db.register(entry)
        new_monomers.append(entry)
        resolved.append(Fragment(
            smiles=frag.smiles, symbol=sym, entry=entry,
            r3_partners=frag.r3_partners, r3_type=frag.r3_type,
            cross_bonds=frag.cross_bonds,
        ))

    if new_monomer_log and new_monomers:
        _log_new_monomers(new_monomers, Path(new_monomer_log))

    # ── Separate PEPTIDE and CHEM monomers ────────────────────────────────
    def _is_chem(frag: Fragment) -> bool:
        return bool(frag.entry and frag.entry.get('polymerType') == 'CHEM')

    peptide_idxs = [i for i, f in enumerate(resolved) if not _is_chem(f)]
    chem_idxs    = [i for i, f in enumerate(resolved) if _is_chem(f)]

    # Position maps: orig resolved index → chain position (1-based)
    peptide_pos_map: dict[int, int] = {orig: pos + 1 for pos, orig in enumerate(peptide_idxs)}
    chem_chain_map: dict[int, str]  = {orig: f'CHEM{i + 1}' for i, orig in enumerate(chem_idxs)}

    # ── Build polymer chain strings ────────────────────────────────────────
    peptide_tokens = [_helm_symbol(resolved[i].symbol) for i in peptide_idxs]
    monomer_str = '.'.join(peptide_tokens)

    chem_chain_strs: list[str] = []
    for orig in chem_idxs:
        cid = chem_chain_map[orig]
        sym = _helm_symbol(resolved[orig].symbol)
        chem_chain_strs.append(f'{cid}{{{sym}}}')

    # ── Build connection lines ─────────────────────────────────────────────
    connections: list[str] = []
    if cyclic:
        connections.append(f'{chain_id},{chain_id},1:R1-{len(peptide_idxs)}:R2')

    # Symmetric R3–R3 connections (disulfide) — both partners in PEPTIDE1
    seen_r3: set[frozenset] = set()
    for i, frag in enumerate(resolved):
        for partner_idx in frag.r3_partners:
            pair: frozenset = frozenset([i, partner_idx])
            if pair not in seen_r3:
                seen_r3.add(pair)
                pos1, pos2 = sorted([peptide_pos_map[i], peptide_pos_map[partner_idx]])
                connections.append(f'{chain_id},{chain_id},{pos1}:R3-{pos2}:R3')

    # Asymmetric cross-bonds (click chem: PEPTIDE↔CHEM) — emit each pair once
    seen_xb: set[frozenset] = set()
    for i, frag in enumerate(resolved):
        for bond in frag.cross_bonds:
            j = bond.partner_idx
            key: frozenset = frozenset([(i, bond.my_rgroup), (j, bond.partner_rgroup)])
            if key in seen_xb:
                continue
            seen_xb.add(key)

            if i in peptide_pos_map:
                chain_i, pos_i = chain_id, peptide_pos_map[i]
            else:
                chain_i, pos_i = chem_chain_map[i], 1

            if j in peptide_pos_map:
                chain_j, pos_j = chain_id, peptide_pos_map[j]
            else:
                chain_j, pos_j = chem_chain_map[j], 1

            connections.append(
                f'{chain_i},{chain_j},{pos_i}:{bond.my_rgroup}-{pos_j}:{bond.partner_rgroup}'
            )

    # ── Assemble HELM string ───────────────────────────────────────────────
    all_chains = [f'{chain_id}{{{monomer_str}}}'] + chem_chain_strs
    helm = f'{"|".join(all_chains)}${"|".join(connections)}$$V2.0'
    return helm, new_monomers


def _log_new_monomers(new_entries: list[dict], log_path: Path) -> None:
    """Persist pending monomers to a registry JSON (dict keyed by UNK_* symbol).

    Format::

        {
          "UNK_a3f9c2": {
            "smiles": "...",
            "name": null,
            "assigned_symbol": null,
            "needs_review": true,
            "first_seen": "2026-06-13"
          }
        }

    Existing entries are never overwritten; new symbols are merged in.
    """
    registry: dict[str, dict] = {}
    if log_path.exists():
        with open(log_path) as f:
            try:
                registry = json.load(f)
            except json.JSONDecodeError:
                registry = {}

    today = date.today().isoformat()
    for entry in new_entries:
        sym = entry['symbol']
        if sym not in registry:
            registry[sym] = {
                'smiles':          entry['smiles'],
                'name':            None,
                'assigned_symbol': None,
                'needs_review':    True,
                'first_seen':      today,
            }

    with open(log_path, 'w') as f:
        json.dump(registry, f, indent=2)
