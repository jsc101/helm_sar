"""
Chemistry-aware Zappo colour scheme for arbitrary monomer symbols.

Strategy
--------
1. Single-letter standard amino acids → canonical pyMSAviz Zappo colours
   (sourced from moshi4/pyMSAviz/config/color_schemes.tsv)
2. Multi-letter / non-standard → look up MonomerDB entry, build RDKit mol,
   classify by chemical structure (aromatic rings, charged groups, LogP, etc.)
3. D-amino acids and N-methylated residues → base colour of the sidechain class,
   slightly desaturated (D) or darkened (N-methyl) as a visual modifier
4. Name-heuristic fallback when no DB entry is available

Public API
----------
residue_color(symbol, db=None)  → hex string
category_of(symbol, db=None)    → category name string
LEGEND_CATEGORIES               → list[(hex, display_label, example_str)]
clear_cache()                   → wipe memoised results
"""

from __future__ import annotations

import re
from typing import Optional

# ──────────────────────────────────────────────────────────────────────────────
# Canonical pyMSAviz Zappo colours (single-letter standard AA)
# Source: moshi4/pyMSAviz/src/pymsaviz/config/color_schemes.tsv
# ──────────────────────────────────────────────────────────────────────────────

ZAPPO: dict[str, str] = {
    'A': '#FFAFAF',  # aliphatic
    'C': '#FFFF00',  # cysteine
    'D': '#FF0000',  # negative
    'E': '#FF0000',  # negative
    'F': '#FFC800',  # aromatic
    'G': '#FF00FF',  # conformational / small
    'H': '#6464FF',  # positive
    'I': '#FFAFAF',  # aliphatic
    'K': '#6464FF',  # positive
    'L': '#FFAFAF',  # aliphatic
    'M': '#FFAFAF',  # aliphatic
    'N': '#00FF00',  # polar
    'P': '#FF00FF',  # conformational / cyclic
    'Q': '#00FF00',  # polar
    'R': '#6464FF',  # positive
    'S': '#00FF00',  # polar
    'T': '#00FF00',  # polar
    'V': '#FFAFAF',  # aliphatic
    'W': '#FFC800',  # aromatic
    'Y': '#FFC800',  # aromatic
}

# ── Category → canonical colour ───────────────────────────────────────────────

_CAT_COLOR: dict[str, str] = {
    'aliphatic':      '#FFAFAF',
    'aromatic':       '#FFC800',
    'positive':       '#6464FF',
    'negative':       '#FF0000',
    'polar':          '#00FF00',
    'cysteine':       '#FFFF00',
    'conformational': '#FF00FF',
    'unknown':        '#CCCCCC',
}

# ── Legend entries: (hex, label, short_examples_str) ─────────────────────────

LEGEND_CATEGORIES: list[tuple[str, str, str]] = [
    ('#FFAFAF', 'Aliphatic',       'A, V, L, I, M, meL, Me_dL'),
    ('#FFC800', 'Aromatic',        'F, W, Y, Phe_4F, 1Nal'),
    ('#6464FF', 'Positive (+)',    'K, R, H, Orn, Dab'),
    ('#FF0000', 'Negative (−)',    'D, E, gGlu, Asp_OtBu'),
    ('#00FF00', 'Polar',           'S, T, N, Q, Hyp, Thr_tBu'),
    ('#FFFF00', 'Cysteine',        'C, Pen, dC'),
    ('#FF00FF', 'Conform./Small',  'G, P, Pip, Sar, Aib, meG'),
    ('#CCCCCC', 'Non-standard',    'PEG2, C18diacid, …'),
]

# ──────────────────────────────────────────────────────────────────────────────
# RDKit chemical classifier
# ──────────────────────────────────────────────────────────────────────────────

# Pre-compiled SMARTS (module-level to avoid repeated compilation)
try:
    from rdkit import Chem as _Chem
    _PAT = {
        'cooh':    _Chem.MolFromSmarts('C(=O)[OH]'),
        'nh2':     _Chem.MolFromSmarts('[NH2]'),
        'guan':    _Chem.MolFromSmarts('NC(=N)N'),          # guanidinium (Arg)
        'imid':    _Chem.MolFromSmarts('[nH]1ccnc1'),        # imidazole (His)
        'thiol':   _Chem.MolFromSmarts('[SH]'),
        'n_ring':  _Chem.MolFromSmarts('[N;R;!a]'),          # N in non-aromatic ring (Pro/Pip)
        'oh_sc':   _Chem.MolFromSmarts('[OH;!$(OC=O)]'),     # OH not in carboxyl (Ser/Thr/Hyp)
        'prim_am': _Chem.MolFromSmarts('[NH2]C=O'),          # primary amide sidechain (Asn/Gln)
    }
    _RDKIT_OK = True
except Exception:
    _PAT = {}
    _RDKIT_OK = False


def _classify_mol(mol) -> str:
    """
    Classify a capped RDKit Mol into a Zappo category.

    Rules (priority order):
      1. Aromatic ring present                 → aromatic
      2. Thiol present                         → cysteine
      3. >1 carboxylic acid group              → negative  (Asp/Glu-like)
      4. >1 primary amine, guanidinium, imidaz → positive  (Lys/Arg/His-like)
      5. N in non-aromatic ring OR MW < 96     → conformational (Pro/Pip/Gly-like)
      6. Sidechain OH or primary amide         → polar
      7. LogP ≥ 0.0                            → aliphatic
      8. otherwise                             → polar
    """
    if mol is None or not _RDKIT_OK:
        return 'unknown'
    try:
        from rdkit.Chem import rdMolDescriptors, Descriptors
        if rdMolDescriptors.CalcNumAromaticRings(mol) > 0:
            return 'aromatic'
        if mol.HasSubstructMatch(_PAT['thiol']):
            return 'cysteine'
        if len(mol.GetSubstructMatches(_PAT['cooh'])) > 1:
            return 'negative'
        n_nh2 = len(mol.GetSubstructMatches(_PAT['nh2']))
        if (n_nh2 > 1
                or mol.HasSubstructMatch(_PAT['guan'])
                or mol.HasSubstructMatch(_PAT['imid'])):
            return 'positive'
        if mol.HasSubstructMatch(_PAT['n_ring']) or Descriptors.MolWt(mol) < 96:
            return 'conformational'
        if mol.HasSubstructMatch(_PAT['oh_sc']) or mol.HasSubstructMatch(_PAT['prim_am']):
            return 'polar'
        return 'aliphatic' if Descriptors.MolLogP(mol) >= 0.0 else 'polar'
    except Exception:
        return 'unknown'


# ──────────────────────────────────────────────────────────────────────────────
# Colour modifier helpers
# ──────────────────────────────────────────────────────────────────────────────

def _blend(hex_color: str, target: tuple, fraction: float) -> str:
    """Blend hex_color toward an RGB target by `fraction` (0=original, 1=target)."""
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)
    tr, tg, tb = target
    r2 = int(r + (tr - r) * fraction)
    g2 = int(g + (tg - g) * fraction)
    b2 = int(b + (tb - b) * fraction)
    return f'#{max(0,min(255,r2)):02x}{max(0,min(255,g2)):02x}{max(0,min(255,b2)):02x}'


_GREY = (180, 180, 180)
_DARK = (60, 60, 60)


def _apply_modifier(base_color: str, is_d: bool, is_nmethyl: bool) -> str:
    color = base_color
    if is_d:
        color = _blend(color, _GREY, 0.22)    # subtle desaturation for D-AAs
    if is_nmethyl:
        color = _blend(color, _DARK, 0.15)    # slightly darker for N-methylated
    return color


# ──────────────────────────────────────────────────────────────────────────────
# Symbol parsing helpers
# ──────────────────────────────────────────────────────────────────────────────

_D_PREFIX  = re.compile(r'^(d([A-Z])|D-)', )   # dL, dA, D-Phe
_ME_PREFIX = re.compile(r'^[Mm]e([A-Z_])',  )   # meL, meA, Me_Bmt


def _parse_symbol(sym: str) -> tuple[str, bool, bool]:
    """
    Return (base_symbol, is_d_aminoacid, is_n_methylated).
    base_symbol strips the D-/me- prefix; the original symbol is kept for DB lookup.
    """
    s = sym.strip()
    is_d = bool(_D_PREFIX.match(s))
    is_me = bool(_ME_PREFIX.match(s))
    # Base for single-letter AA lookup: strip me/D prefixes
    m_d  = _D_PREFIX.match(s)
    m_me = _ME_PREFIX.match(s)
    if m_d:
        base = m_d.group(2) if m_d.group(2) else s[len(m_d.group(0)):]
    elif m_me:
        base = m_me.group(1)
    else:
        base = s
    return base, is_d, is_me


# ──────────────────────────────────────────────────────────────────────────────
# Name-heuristic fallback (no DB entry)
# ──────────────────────────────────────────────────────────────────────────────

def _name_heuristic(sym: str) -> str:
    s = sym.lower().strip('_')
    # Linkers / fatty acids
    if any(k in s for k in ('peg', 'c18', 'fatty', 'lipid', 'glu_c18', 'diacid')):
        return _CAT_COLOR['unknown']
    base, is_d, is_me = _parse_symbol(sym)
    std_color = ZAPPO.get(base[:1].upper() if base else 'X')
    if std_color:
        return _apply_modifier(std_color, is_d, is_me)
    return _CAT_COLOR['unknown']


# ──────────────────────────────────────────────────────────────────────────────
# Main public API
# ──────────────────────────────────────────────────────────────────────────────

_color_cache: dict[str, str] = {}
_cat_cache:   dict[str, str] = {}


def residue_color(symbol: str, db=None) -> str:
    """
    Return hex background colour for *symbol* using the chemistry-aware Zappo scheme.

    Parameters
    ----------
    symbol : monomer symbol (e.g. 'A', 'dL', 'meL', 'Phe_4F', 'gGlu')
    db     : MonomerDB instance; if None, falls back to name heuristics for
             non-standard monomers
    """
    if symbol in _color_cache:
        return _color_cache[symbol]

    # Fast path: unmodified single-letter standard AA
    if len(symbol) == 1 and symbol in ZAPPO:
        _color_cache[symbol] = ZAPPO[symbol]
        _cat_cache[symbol]   = _zappo_category(symbol)
        return ZAPPO[symbol]

    base, is_d, is_me = _parse_symbol(symbol)

    # Fast path: single-letter base after prefix stripping
    if len(base) == 1 and base in ZAPPO:
        base_color = ZAPPO[base]
        color = _apply_modifier(base_color, is_d, is_me)
        _cat_cache[symbol] = _zappo_category(base)
        _color_cache[symbol] = color
        return color

    # DB + RDKit path
    entry = db.find_by_symbol(symbol) if db is not None else None
    if entry is not None and _RDKIT_OK:
        try:
            from scripts.rdkit_bridge import monomer_to_mol
            mol = monomer_to_mol(entry)
            cat = _classify_mol(mol)
            base_color = _CAT_COLOR[cat]
            color = _apply_modifier(base_color, is_d, is_me)
            _cat_cache[symbol] = cat
            _color_cache[symbol] = color
            return color
        except Exception:
            pass

    # Name-heuristic fallback
    color = _name_heuristic(symbol)
    _color_cache[symbol] = color
    return color


def category_of(symbol: str, db=None) -> str:
    """Return the category name for *symbol* (triggers colour computation if needed)."""
    if symbol not in _cat_cache:
        residue_color(symbol, db)
    return _cat_cache.get(symbol, 'unknown')


def clear_cache() -> None:
    _color_cache.clear()
    _cat_cache.clear()


def _zappo_category(single_letter: str) -> str:
    color = ZAPPO.get(single_letter, '#CCCCCC')
    for cat, c in _CAT_COLOR.items():
        if c == color:
            return cat
    return 'unknown'


# ──────────────────────────────────────────────────────────────────────────────
# Text contrast helper (shared with alignment_viz)
# ──────────────────────────────────────────────────────────────────────────────

def text_color(bg_hex: str) -> str:
    """Return '#111111' or '#ffffff' for good contrast on bg_hex."""
    r = int(bg_hex[1:3], 16)
    g = int(bg_hex[3:5], 16)
    b = int(bg_hex[5:7], 16)
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    return '#111111' if lum > 140 else '#ffffff'
