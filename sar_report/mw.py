"""Molecular-weight calculation for a HELM structure."""
from __future__ import annotations


def calc_mw(obj) -> float | None:
    """MW = Σ(residue MW) − n_bonds × 18.015, where
    n_bonds = (n_residues − 1) + n_explicit_connections.

    Returns None if any residue is missing from the MonomerDB.
    """
    desc = obj.all_monomer_descriptors()
    if not desc:
        return None
    if any(not v['in_db'] for v in desc.values()):
        return None
    total = sum(v['descriptors']['MW'] for v in desc.values())
    n = len(desc)
    n_explicit = len(obj.data.get('connectivity_graph', []))
    return total - ((n - 1) + n_explicit) * 18.015
