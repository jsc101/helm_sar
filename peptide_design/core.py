from __future__ import annotations
from dataclasses import dataclass, field
from typing import List


@dataclass
class Sidechain:
    site: int           # HELM position (1-based)
    bonds: List[str]    # b1..bN, each "proxRg-distRg"; len == len(monomers)
    monomers: List[str] # SC1..SCN symbols, proximal→distal


@dataclass
class Compound:
    name: str
    main_chain: List[str]
    sidechains: List[Sidechain] = field(default_factory=list)


def parse_row(row: dict) -> Compound:
    name = str(row.get("Name", ""))
    mc_cols = sorted([k for k in row if str(k).isdigit()], key=int)
    main_chain = [str(row[c]).strip() for c in mc_cols if str(row.get(c, "")).strip()]
    sidechains = []
    for suffix in ("", "_2", "_3"):
        site_key = f"Site{suffix}"
        raw = str(row.get(site_key, "")).strip()
        if not raw:
            continue
        site = int(raw)
        bonds, monomers = [], []
        k = 1
        while True:
            bval = str(row.get(f"b{k}{suffix}", "")).strip()
            sval = str(row.get(f"SC{k}{suffix}", "")).strip()
            if not bval or not sval:
                break
            bonds.append(bval)
            monomers.append(sval)
            k += 1
        if bonds:
            sidechains.append(Sidechain(site=site, bonds=bonds, monomers=monomers))
    return Compound(name=name, main_chain=main_chain, sidechains=sidechains)
