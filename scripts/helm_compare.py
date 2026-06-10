"""
Structural equivalence comparison for HELM strings.

For linear peptides: compare monomer canonical SMILES sequences directly.
For cyclic peptides: try all N rotations (and reverse) to find a match.
Comparison uses MonomerDB canonical SMILES, not symbol strings — handles
synonyms and D/L naming differences between databases.
"""

from __future__ import annotations

from typing import Optional

from monomer_db.monomer_db import MonomerDB

# Module-level singleton to avoid reloading JSON on every call
_DB: Optional[MonomerDB] = None


def _get_db(monomer_db=None, extra_sources=None) -> MonomerDB:
    """Return the provided db, or fall back to the module-level singleton."""
    if monomer_db is not None:
        return monomer_db
    global _DB
    if _DB is None:
        _DB = MonomerDB(extra_sources=extra_sources)
    return _DB


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _canonicalize_smiles(smiles: str) -> Optional[str]:
    """Return RDKit canonical SMILES (isomericSmiles=True) or None."""
    try:
        from rdkit import Chem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _get_canonical_sequence(
    helm_obj, db: MonomerDB
) -> list[tuple[str, Optional[str]]]:
    """
    Extract monomers from the primary PEPTIDE chain of a HELMObject.

    For each monomer symbol:
      1. Try ``db.find_by_symbol(symbol)`` to get the entry.
      2. If the entry has a SMILES, canonicalize with RDKit (stereo=True)
         and use it as the comparison key.
      3. If no entry or no parseable SMILES, use the symbol string as a
         fallback key (the symbol will also appear in the unresolved list).

    Returns
    -------
    list of (symbol, canonical_smiles_or_None)
        None means the symbol could not be resolved to a canonical SMILES.
    """
    chains = helm_obj.data.get('_chains') or helm_obj.data.get('chains', [])
    peptide_chains = [c for c in chains if c.get('type') == 'PEPTIDE']
    if not peptide_chains:
        return []

    primary = max(peptide_chains, key=lambda c: len(c['monomers']))
    result: list[tuple[str, Optional[str]]] = []

    for m in primary['monomers']:
        symbol = m['symbol']
        # Prefer the entry already attached by the parser; fall back to DB lookup
        entry = m.get('entry') or db.find_by_symbol(symbol)
        canonical: Optional[str] = None
        if entry:
            raw_smiles = entry.get('smiles') or entry.get('SMILES', '')
            if raw_smiles:
                canonical = _canonicalize_smiles(raw_smiles)
        result.append((symbol, canonical))

    return result


def _compare_pair(
    pair1: tuple[str, Optional[str]],
    pair2: tuple[str, Optional[str]],
) -> tuple[bool, Optional[str]]:
    """
    Compare two (symbol, canonical_smiles) pairs.

    Returns (match: bool, unresolved_symbol: str | None).
    A symbol is "unresolved" when either side has no canonical SMILES.
    """
    sym1, can1 = pair1
    sym2, can2 = pair2

    if can1 is not None and can2 is not None:
        return (can1 == can2), None
    # Fallback: compare symbols case-insensitively
    unresolved = sym1 if can1 is None else sym2
    return (sym1.lower() == sym2.lower()), unresolved


def _compare_sequences(
    seq1: list[tuple[str, Optional[str]]],
    seq2: list[tuple[str, Optional[str]]],
) -> tuple[bool, list[str]]:
    """
    Compare two monomer sequences element-by-element.

    Returns
    -------
    (all_match, unresolved_symbols)
        ``all_match`` is True only if every position matches.
        ``unresolved_symbols`` lists symbols that fell back to string comparison.
    """
    if len(seq1) != len(seq2):
        return False, []

    all_match = True
    unresolved: list[str] = []

    for p1, p2 in zip(seq1, seq2):
        match, unresolv = _compare_pair(p1, p2)
        if not match:
            all_match = False
        if unresolv and unresolv not in unresolved:
            unresolved.append(unresolv)

    return all_match, unresolved


def _find_cyclic_match(
    seq1: list[tuple[str, Optional[str]]],
    seq2: list[tuple[str, Optional[str]]],
) -> dict:
    """
    Try all N rotations of seq2, and all N rotations of the reversed seq2,
    looking for a sequence that matches seq1.

    Returns a dict with keys:
        found     : bool
        offset    : int | None   (rotation index into seq2 that matched)
        reversed  : bool | None  (True if the match is in reverse orientation)
        unresolved: list[str]    (symbols that fell back to string comparison)
    """
    n = len(seq1)

    # Forward rotations
    for offset in range(n):
        rotated = seq2[offset:] + seq2[:offset]
        match, unresolved = _compare_sequences(seq1, rotated)
        if match:
            return {"found": True, "offset": offset, "reversed": False,
                    "unresolved": unresolved}

    # Reverse rotations
    seq2_rev = list(reversed(seq2))
    for offset in range(n):
        rotated = seq2_rev[offset:] + seq2_rev[:offset]
        match, unresolved = _compare_sequences(seq1, rotated)
        if match:
            return {"found": True, "offset": offset, "reversed": True,
                    "unresolved": unresolved}

    return {"found": False, "offset": None, "reversed": None, "unresolved": []}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#: Topology classes that require rotation-aware comparison
_CYCLIC_TOPOLOGIES = frozenset({
    'monocyclic_disulfide',
    'monocyclic_backbone',
    'bicyclic',
})


def compare(
    helm1: str,
    helm2: str,
    monomer_db=None,
    extra_sources=None,
) -> dict:
    """
    Structurally compare two HELM strings.

    Parameters
    ----------
    helm1, helm2 : str
        HELM V2 strings to compare.
    monomer_db : MonomerDB, optional
        Pre-built MonomerDB instance.  If None, the module-level singleton is
        used (created lazily).
    extra_sources : list[str], optional
        Extra JSON paths forwarded to MonomerDB if the singleton is being
        created for the first time.

    Returns
    -------
    dict with keys:
        match               : bool
        topology_match      : bool   — both strings parse to the same topology class
        n_residues_match    : bool
        sequence_match      : bool   — monomer-by-monomer canonical SMILES agree
        cyclic_offset       : int | None
        cyclic_reversed     : bool | None
        unresolved_monomers : list[str]
        details             : str
    """
    from scripts.helm_parser import HELMParser

    db = _get_db(monomer_db, extra_sources)

    obj1 = HELMParser.parse(helm1)
    obj2 = HELMParser.parse(helm2)

    topo1 = obj1.data.get('topology_class', 'linear')
    topo2 = obj2.data.get('topology_class', 'linear')
    topology_match = (topo1 == topo2)

    seq1 = _get_canonical_sequence(obj1, db)
    seq2 = _get_canonical_sequence(obj2, db)

    n_residues_match = (len(seq1) == len(seq2))

    # Template for a failed result
    def _fail(details: str, **extra) -> dict:
        base = {
            "match": False,
            "topology_match": topology_match,
            "n_residues_match": n_residues_match,
            "sequence_match": False,
            "cyclic_offset": None,
            "cyclic_reversed": None,
            "unresolved_monomers": [],
            "details": details,
        }
        base.update(extra)
        return base

    if not topology_match:
        return _fail(
            f"Topology mismatch: '{topo1}' vs '{topo2}'"
        )

    if not n_residues_match:
        return _fail(
            f"Length mismatch: {len(seq1)} vs {len(seq2)} residues"
        )

    is_cyclic = topo1 in _CYCLIC_TOPOLOGIES

    if not is_cyclic:
        seq_match, unresolved = _compare_sequences(seq1, seq2)
        if not seq_match:
            return _fail(
                "Linear sequences differ (canonical SMILES mismatch)",
                unresolved_monomers=unresolved,
            )
        details = "Linear sequences are identical"
        if unresolved:
            details += f" (symbol fallback used for: {unresolved})"
        return {
            "match": True,
            "topology_match": True,
            "n_residues_match": True,
            "sequence_match": True,
            "cyclic_offset": None,
            "cyclic_reversed": None,
            "unresolved_monomers": unresolved,
            "details": details,
        }

    # Cyclic comparison
    best = _find_cyclic_match(seq1, seq2)
    if not best["found"]:
        return _fail(
            f"Cyclic ({topo1}) sequences differ across all {len(seq1)} "
            "rotations and reverse orientations",
        )

    offset = best["offset"]
    rev = best["reversed"]
    unresolved = best["unresolved"]
    dir_str = "reversed + " if rev else ""
    details = (
        f"Cyclic ({topo1}) match found at {dir_str}rotation offset {offset}"
    )
    if unresolved:
        details += f" (symbol fallback used for: {unresolved})"

    return {
        "match": True,
        "topology_match": True,
        "n_residues_match": True,
        "sequence_match": True,
        "cyclic_offset": offset,
        "cyclic_reversed": rev,
        "unresolved_monomers": unresolved,
        "details": details,
    }


def canonical_cyclic_helm(
    helm_string: str,
    monomer_db=None,
    extra_sources=None,
) -> str:
    """
    Return a HELM string rotated to canonical (lexicographically smallest)
    position, useful for deduplication of cyclic compound libraries.

    The canonical position is determined by the lexicographically smallest
    tuple of per-monomer keys, where each key is:
      - the canonical SMILES if the monomer is resolved, or
      - the symbol string as a fallback.

    Only the primary PEPTIDE chain is reordered; connections and other
    sections are rebuilt to reflect the new numbering.

    Parameters
    ----------
    helm_string : str
        HELM V2 string (should be a monocyclic or bicyclic topology).
    monomer_db : MonomerDB, optional
    extra_sources : list[str], optional

    Returns
    -------
    str
        HELM string with the primary PEPTIDE chain rotated to canonical form.
        Non-cyclic HELM strings are returned unchanged.
    """
    from scripts.helm_parser import HELMParser

    db = _get_db(monomer_db, extra_sources)
    obj = HELMParser.parse(helm_string)

    topo = obj.data.get('topology_class', 'linear')
    if topo not in _CYCLIC_TOPOLOGIES:
        return helm_string

    seq = _get_canonical_sequence(obj, db)
    n = len(seq)
    if n == 0:
        return helm_string

    # Build sort keys: use canonical SMILES when available, else symbol
    keys = [can if can is not None else sym for sym, can in seq]

    # Find rotation with lexicographically smallest key tuple
    best_offset = 0
    best_tuple = tuple(keys)
    for offset in range(1, n):
        rotated = tuple(keys[offset:] + keys[:offset])
        if rotated < best_tuple:
            best_tuple = rotated
            best_offset = offset

    if best_offset == 0:
        return helm_string

    # Reconstruct the HELM string with rotated primary chain
    chains = obj.data.get('_chains') or obj.data.get('chains', [])
    peptide_chains = [c for c in chains if c.get('type') == 'PEPTIDE']
    if not peptide_chains:
        return helm_string

    primary = max(peptide_chains, key=lambda c: len(c['monomers']))
    old_monomers = primary['monomers']
    rotated_monomers = old_monomers[best_offset:] + old_monomers[:best_offset]
    rotated_symbols = [m['symbol'] for m in rotated_monomers]

    # Build new monomer notation (bracket non-standard symbols)
    def _fmt(sym: str) -> str:
        # Standard single-letter amino acids need no brackets in HELM
        # Multi-char or non-alpha symbols need square brackets
        if len(sym) == 1 and sym.isalpha():
            return sym
        return f'[{sym}]'

    chain_id = primary['chain_id']
    new_block = f"{chain_id}" + "{" + ".".join(_fmt(s) for s in rotated_symbols) + "}"

    # Rebuild all polymer blocks; replace primary chain with rotated version
    blocks = []
    for chain in chains:
        if chain['chain_id'] == primary['chain_id']:
            blocks.append(new_block)
        else:
            syms = [_fmt(m['symbol']) for m in chain['monomers']]
            blocks.append(
                f"{chain['chain_id']}" + "{" + ".".join(syms) + "}"
            )
    polymer_section = "|".join(blocks)

    # Rebuild connections, adjusting positions within the rotated primary chain
    conn_graph = obj.data.get('connectivity_graph', [])
    conn_parts = []
    for conn in conn_graph:
        fc = conn['from_chain']
        tc = conn['to_chain']
        fp = conn['from_pos']
        tp = conn['to_pos']

        # Adjust 1-based positions for rotated primary chain
        if fc == chain_id:
            fp = (fp - 1 - best_offset) % n + 1
        if tc == chain_id:
            tp = (tp - 1 - best_offset) % n + 1

        conn_parts.append(
            f"{fc},{tc},{fp}:{conn['from_rgroup']}-{tp}:{conn['to_rgroup']}"
        )
    conn_section = "|".join(conn_parts)

    # Reassemble full HELM string (groups and annotations left empty, V2.0)
    return f"{polymer_section}${conn_section}$$$V2.0"
