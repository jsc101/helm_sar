from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from peptide_design.core import Compound, Sidechain

_DB = None

def _db():
    global _DB
    if _DB is None:
        from monomer_db.monomer_db import MonomerDB
        _DB = MonomerDB()
    return _DB


def fmt_sym(s: str) -> str:
    return f"[{s}]" if len(s) > 1 else s


def _polymer_type(symbol: str) -> str:
    entry = _db().find_by_symbol(symbol)
    if entry is None:
        return "UNKNOWN"
    return entry.get("polymerType", "PEPTIDE")


def _can_collapse(sc: Sidechain) -> bool:
    if not sc.bonds or not sc.monomers:
        return False
    if any(_polymer_type(m) != "PEPTIDE" for m in sc.monomers):
        return False
    if sc.bonds[0].split("-")[1] != "R1":
        return False
    if any(b != "R2-R1" for b in sc.bonds[1:]):
        return False
    return True


def _build_collapsed(sc: Sidechain, chain_n: int) -> tuple[str, str]:
    body = ".".join(fmt_sym(m) for m in sc.monomers)
    chain_id = f"PEPTIDE{chain_n}"
    chain_def = f"{chain_id}{{{body}}}"
    prox_rg = sc.bonds[0].split("-")[0]
    conn = f"{chain_id},PEPTIDE1,1:R1-{sc.site}:{prox_rg}"
    return chain_def, conn


def _build_expanded(sc: Sidechain, peptide_n: int, chem_n: int
                    ) -> tuple[list[str], list[str], int, int]:
    chain_defs: list[str] = []
    connections: list[str] = []
    prev_chain_id: str | None = None

    for k, (bond, monomer) in enumerate(zip(sc.bonds, sc.monomers)):
        ptype = _polymer_type(monomer)
        if ptype == "CHEM":
            chain_id = f"CHEM{chem_n}"
            chem_n += 1
        else:
            chain_id = f"PEPTIDE{peptide_n}"
            peptide_n += 1

        chain_defs.append(f"{chain_id}{{{fmt_sym(monomer)}}}")

        prox_rg, dist_rg = bond.split("-")
        if k == 0:
            conn = f"{chain_id},PEPTIDE1,1:{dist_rg}-{sc.site}:{prox_rg}"
        else:
            conn = f"{chain_id},{prev_chain_id},1:{dist_rg}-1:{prox_rg}"
        connections.append(conn)
        prev_chain_id = chain_id

    return chain_defs, connections, peptide_n, chem_n


def generate_helm(compound: Compound) -> str:
    body = ".".join(fmt_sym(s) for s in compound.main_chain)
    chains = [f"PEPTIDE1{{{body}}}"]
    connections: list[str] = []
    peptide_n = 2
    chem_n = 1

    for sc in compound.sidechains:
        if _can_collapse(sc):
            chain_def, conn = _build_collapsed(sc, peptide_n)
            chains.append(chain_def)
            connections.append(conn)
            peptide_n += 1
        else:
            cdefs, conns, peptide_n, chem_n = _build_expanded(sc, peptide_n, chem_n)
            chains.extend(cdefs)
            connections.extend(conns)

    chain_part = "|".join(chains)
    if connections:
        return f"{chain_part}${'|'.join(connections)}$$$V2.0"
    return f"{chain_part}$$$$V2.0"
