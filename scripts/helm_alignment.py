"""
HELMAlignment — cyclic peptide alignment via per-residue descriptor similarity.

Rotation model
--------------
A cyclic peptide of length n has n valid rotations (and optionally n reversed).
Rotation k shifts the sequence so that original position k+1 becomes position 1,
with all connection positions remapped modulo n.

After alignment, the rotated HELM has the same underlying molecule as the original
(verifiable via Tanimoto = 1.0 or helm_compare.compare returning match=True).

Scoring
-------
Each residue → descriptor vector [MW, LogP, TPSA, HBD, HBA, RotBonds, AromaticRings, Chiral]
Each value is divided by its scale factor before the cosine dot-product, so each
dimension contributes comparably. Cosine similarity is direction-only (magnitude-
invariant per vector), so the scale factors control the relative weight each
descriptor has on the angular distance between residues.
Score(rotation k) = mean cosine similarity between corresponding position vectors.

Supports monocyclic_backbone, monocyclic_disulfide, and bicyclic topologies.
For bicyclics the full rotation space (0..n-1) is searched; connectivity is
preserved automatically because connection positions are remapped with the chain.

Public API
----------
rotate_chain(obj, k)            → HELMObject   (rotated by k)
helm_from_obj(obj)              → str           (reconstruct HELM V2 string)
HELMAlignment(reference)        → aligner
  .score_rotation(query, k)     → float
  .align(query)                 → AlignmentResult
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from scripts.helm_engine import HELMObject

# Descriptor keys and per-descriptor scale factors (approximate max for amino acids)
# 'Chiral' is a synthetic feature: 1.0 for D-amino acids (symbol starts with 'd' or 'D-'),
# 0.0 for L/achiral.  Helps distinguish enantiomers whose RDKit descriptors are identical.
_DESC_KEYS  = ['MW',    'LogP', 'TPSA',  'HBD', 'HBA', 'RotBonds', 'AromaticRings', 'Chiral']
_DESC_SCALE = [500.0,    8.0,   200.0,    5.0,   10.0,   10.0,       3.0,             1.0    ]


def _chirality_flag(symbol: str) -> float:
    """Return 1.0 for D-amino acids, 0.0 for L/achiral."""
    s = symbol.strip()
    # Explicit D-prefix conventions used in Pistoia / CycPeptMPDB
    if s.startswith('d') and len(s) > 1 and s[1].isupper():
        return 1.0   # e.g. dL, dA, dK
    if s.startswith('D-'):
        return 1.0   # e.g. D-Phe, D-Pro
    return 0.0


# ---------------------------------------------------------------------------
# Chain rotation
# ---------------------------------------------------------------------------

def rotate_chain(obj: HELMObject, k: int, chain_id: Optional[str] = None) -> HELMObject:
    """
    Return a new HELMObject with the primary (or specified) chain rotated by k positions.

    Rotation k means: original monomer at position k+1 becomes position 1.
    All connection positions referencing the rotated chain are remapped modulo n.
    Secondary chains are not rotated.

    Parameters
    ----------
    obj      : HELMObject to rotate
    k        : rotation offset (0 = identity; wraps modulo chain length)
    chain_id : target chain id; None → primary (longest PEPTIDE) chain
    """
    chains  = obj.data.get('_chains', [])
    primary = obj.get_chain(chain_id)
    if primary is None or len(primary['monomers']) == 0:
        return obj

    pid = primary['chain_id']
    n   = len(primary['monomers'])
    k   = k % n
    if k == 0:
        return obj

    # Rotate monomer list, re-index positions 1..n
    old_monomers = primary['monomers']
    rotated_monomers = [
        {**m, 'pos': i + 1}
        for i, m in enumerate(old_monomers[k:] + old_monomers[:k])
    ]

    new_chains = [
        ({**c, 'monomers': rotated_monomers} if c['chain_id'] == pid else c)
        for c in chains
    ]

    # Remap connection positions in the rotated chain
    conn_graph = obj.data.get('connectivity_graph', [])
    new_conn_graph = []
    for conn in conn_graph:
        fp = (conn['from_pos'] - 1 - k) % n + 1 if conn['from_chain'] == pid else conn['from_pos']
        tp = (conn['to_pos']   - 1 - k) % n + 1 if conn['to_chain']   == pid else conn['to_pos']
        new_conn_graph.append({**conn, 'from_pos': fp, 'to_pos': tp})

    # Rebuild flat connectivity_map
    offset_map: dict[str, int] = {}
    cursor = 0
    for c in new_chains:
        offset_map[c['chain_id']] = cursor
        cursor += len(c['monomers'])

    new_conn_map = [
        (offset_map.get(c['from_chain'], 0) + c['from_pos'] - 1,
         offset_map.get(c['to_chain'],   0) + c['to_pos']   - 1,
         c['bond_type'])
        for c in new_conn_graph
    ]

    new_data = {
        **obj.data,
        '_chains':            new_chains,
        'connectivity_graph': new_conn_graph,
        'connectivity_map':   new_conn_map,
    }
    return HELMObject(new_data)


# ---------------------------------------------------------------------------
# HELM string reconstruction
# ---------------------------------------------------------------------------

def helm_from_obj(obj: HELMObject) -> str:
    """
    Reconstruct a HELM V2 string from a HELMObject.

    Single-letter alpha symbols are written bare; everything else is bracketed.
    Connections use standard HELM format: CHAIN,CHAIN,pos:Rg-pos:Rg
    """
    def _fmt(sym: str) -> str:
        return sym if (len(sym) == 1 and sym.isalpha()) else f'[{sym}]'

    chains = obj.data.get('_chains', [])
    blocks = [
        '{}{{{}}}' .format(c['chain_id'], '.'.join(_fmt(m['symbol']) for m in c['monomers']))
        for c in chains
    ]
    conn_graph = obj.data.get('connectivity_graph', [])
    conn_parts = [
        f"{c['from_chain']},{c['to_chain']},{c['from_pos']}:{c['from_rgroup']}-{c['to_pos']}:{c['to_rgroup']}"
        for c in conn_graph
    ]
    return f"{'|'.join(blocks)}${'|'.join(conn_parts)}$$$V2.0"


# ---------------------------------------------------------------------------
# Descriptor matrix
# ---------------------------------------------------------------------------

def _descriptor_matrix(obj: HELMObject, chain_id: Optional[str] = None) -> np.ndarray:
    """
    Return an (n, 8) float matrix: one row per residue in chain position order.
    Columns: [MW, LogP, TPSA, HBD, HBA, RotBonds, AromaticRings, Chiral], each scale-normalized.
    Residues with no DB entry (unknown monomers) get the zero vector.
    """
    all_desc = obj.all_monomer_descriptors(chain_id)
    chain    = obj.get_chain(chain_id)
    positions = sorted(all_desc.keys())
    monomer_by_pos = {m['pos']: m for m in (chain['monomers'] if chain else [])}
    rows = []
    for pos in positions:
        d   = all_desc[pos]['descriptors']
        sym = monomer_by_pos.get(pos, {}).get('symbol', '')
        row = [d.get(k, 0.0) / s for k, s in zip(_DESC_KEYS[:-1], _DESC_SCALE[:-1])]
        row.append(_chirality_flag(sym) / _DESC_SCALE[-1])
        rows.append(row)
    return np.array(rows, dtype=float) if rows else np.empty((0, len(_DESC_KEYS)))


def _mean_cosine(A: np.ndarray, B: np.ndarray) -> float:
    """
    Mean cosine similarity between corresponding rows of A and B (both n×d).
    Zero vectors contribute 0.0 to the mean.
    """
    if A.shape != B.shape or A.shape[0] == 0:
        return 0.0
    na = np.linalg.norm(A, axis=1, keepdims=True)
    nb = np.linalg.norm(B, axis=1, keepdims=True)
    na = np.where(na == 0, 1.0, na)
    nb = np.where(nb == 0, 1.0, nb)
    dots = np.sum((A / na) * (B / nb), axis=1)
    return float(np.mean(dots))


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class AlignmentResult:
    best_k:       int
    best_score:   float
    all_scores:   list[float]
    rotated:      HELMObject
    rotated_helm: str
    rotated_jpv:  str
    verification: dict        # from helm_compare.compare
    n:            int         # query chain length


# ---------------------------------------------------------------------------
# HELMAlignment
# ---------------------------------------------------------------------------

class HELMAlignment:
    """
    Align a query cyclic peptide to a fixed reference by finding the rotation
    that maximises per-position descriptor cosine similarity.

    Parameters
    ----------
    reference : HELMObject
        The fixed reference peptide. Must have a primary PEPTIDE chain.

    Example
    -------
    >>> aligner = HELMAlignment(ref_obj)
    >>> result  = aligner.align(query_obj)
    >>> print(result.best_k, result.best_score)
    >>> print(result.rotated_jpv)
    """

    def __init__(self, reference: HELMObject, chain_id: Optional[str] = None) -> None:
        self.reference  = reference
        self.chain_id   = chain_id
        self._ref_mat   = _descriptor_matrix(reference, chain_id)

    # ------------------------------------------------------------------ #

    def score_rotation(self, query: HELMObject, k: int) -> float:
        """Cosine similarity score for query rotated by k against reference."""
        rotated = rotate_chain(query, k, self.chain_id)
        q_mat   = _descriptor_matrix(rotated, self.chain_id)
        if q_mat.shape[0] != self._ref_mat.shape[0]:
            return -1.0
        return _mean_cosine(q_mat, self._ref_mat)

    # ------------------------------------------------------------------ #

    def align(
        self,
        query: HELMObject,
        try_reverse: bool = False,
    ) -> AlignmentResult:
        """
        Find the best rotation of query to match reference.

        Parameters
        ----------
        query       : HELMObject to align
        try_reverse : if True, also score the reversed sequence
                      (useful for non-directional ring topologies)

        Returns
        -------
        AlignmentResult with the best rotation, score, JPV, and verification.
        """
        from scripts.helm_compare import compare

        n = len(query.positions(self.chain_id))
        if n == 0:
            raise ValueError("Query has no residues in the target chain.")

        ref_n = self._ref_mat.shape[0]
        if n != ref_n:
            raise ValueError(
                f"Query length {n} ≠ reference length {ref_n}. "
                "Alignment requires equal-length chains."
            )

        # Score all rotations
        scores: list[float] = [self.score_rotation(query, k) for k in range(n)]

        if try_reverse:
            rev_scores = self._score_reversed(query, n)
            if max(rev_scores) > max(scores):
                scores = rev_scores

        best_k     = int(np.argmax(scores))
        best_score = scores[best_k]

        rotated      = rotate_chain(query, best_k, self.chain_id)
        rotated_helm = helm_from_obj(rotated)
        rotated_jpv  = rotated.get_jpv()

        # Verify structural identity: original and rotated should match
        original_helm = helm_from_obj(query)
        verification  = compare(original_helm, rotated_helm)

        return AlignmentResult(
            best_k       = best_k,
            best_score   = best_score,
            all_scores   = scores,
            rotated      = rotated,
            rotated_helm = rotated_helm,
            rotated_jpv  = rotated_jpv,
            verification = verification,
            n            = n,
        )

    def _score_reversed(self, query: HELMObject, n: int) -> list[float]:
        """Score all rotations of the reversed primary chain sequence."""
        chain  = query.get_chain(self.chain_id)
        if chain is None:
            return [-1.0] * n
        rev_monomers = list(reversed(chain['monomers']))
        rev_monomers = [{**m, 'pos': i + 1} for i, m in enumerate(rev_monomers)]
        new_chains   = [
            ({**c, 'monomers': rev_monomers}
             if c['chain_id'] == chain['chain_id'] else c)
            for c in query.data.get('_chains', [])
        ]
        rev_obj = HELMObject({**query.data, '_chains': new_chains})
        return [self.score_rotation(rev_obj, k) for k in range(n)]

    # ------------------------------------------------------------------ #

    @staticmethod
    def score_matrix(
        queries: list[HELMObject],
        reference: HELMObject,
        chain_id: Optional[str] = None,
    ) -> np.ndarray:
        """
        Compute pairwise best-rotation scores for a list of queries vs one reference.
        Returns a 1D array of shape (len(queries),).
        """
        aligner = HELMAlignment(reference, chain_id)
        return np.array([aligner.align(q).best_score for q in queries])
