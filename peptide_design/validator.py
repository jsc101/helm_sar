from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.helm_parser import HELMParser
from monomer_db.monomer_db import MonomerDB

_DB: MonomerDB | None = None

def _db() -> MonomerDB:
    global _DB
    if _DB is None:
        _DB = MonomerDB()
    return _DB


def validate_helm(helm: str) -> tuple[bool, str]:
    """
    Parse the HELM string and check all monomers are in the DB.
    Returns (True, "OK") or (False, error_message).
    """
    try:
        obj = HELMParser.parse(helm)
    except Exception as e:
        return False, f"Parse error: {e}"

    chains = obj.data.get("_chains", [])
    if not chains:
        return False, "Parse error: no polymer chains found in HELM string"

    db = _db()
    unknown = []
    for chain in chains:
        for m in chain.get("monomers", []):
            sym = m.get("symbol", "")
            if sym and db.find_by_symbol(sym) is None:
                unknown.append(sym)

    if unknown:
        return False, f"Unknown monomers: {', '.join(sorted(set(unknown)))}"
    return True, "OK"
