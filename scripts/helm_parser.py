"""
HELMParser — parse a HELM V2 string into a HELMObject.

HELM string format:
    <polymer_blocks>$<connections>$<groups>$<annotations>$<version>

Sections are split on '$'.  Sections 2, 3, 4 are ignored.
"""

from __future__ import annotations

import re
import sys
import os
from pathlib import Path
from typing import Optional

# Make monomer_db importable when running from repo root
_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from monomer_db.monomer_db import MonomerDB
from scripts.helm_engine import HELMObject

# Module-level singleton to avoid reloading JSON on every parse call
_DB: Optional[MonomerDB] = None


def _get_db() -> MonomerDB:
    global _DB
    if _DB is None:
        _DB = MonomerDB()
    return _DB


# ---------------------------------------------------------------------------
# Regex constants
# ---------------------------------------------------------------------------

_CHAIN_RE = re.compile(r'(\w+)\{([^}]+)\}')


def _normalize_symbol(sym: str) -> str:
    """
    Sanitize a monomer symbol so it is safe inside JPV notation.
    - Dashes: strip leading/trailing cap-indicator dashes; replace internal dashes with '_'
    - Parens: replace '(' with '_', drop ')'; parens are reserved for JPV branch notation
    - Arrows (->): replace with 'to' before dash removal
    - Commas inside substituent lists: replace with '_'
    """
    sym = sym.strip()
    sym = sym.lstrip('-').rstrip('-')
    sym = sym.replace('->', 'to')
    sym = sym.replace('(', '_').replace(')', '')
    sym = sym.replace('-', '_')
    sym = sym.replace(',', '_')
    while '__' in sym:
        sym = sym.replace('__', '_')
    return sym.strip('_')
# Standard HELM V2: CHAIN1,CHAIN2,pos1:Rg1-pos2:Rg2
_CONN_RE = re.compile(r'(\w+),(\w+),(\d+):(\w+)-(\d+):(\w+)')
# Alternate format (e.g. Semaglutide): CHAIN1:Rg1(pos1)-CHAIN2:Rg2(pos2)
# Position is optional; when omitted R1→1, R2→last monomer of chain.
_CONN_RE_ALT = re.compile(r'(\w+):(\w+)(?:\((\d+)\))?-(\w+):(\w+)(?:\((\d+)\))?')


# ---------------------------------------------------------------------------
# Bond type inference
# ---------------------------------------------------------------------------

def _infer_bond_type(
    from_chain: str,
    from_rgroup: str,
    to_chain: str,
    to_rgroup: str,
) -> str:
    """
    Infer bond type from rgroup labels and chain types.

    Rules (in priority order):
    1. Both R3 → "disulfide" (or "thioether" if one side is a CHEM chain)
    2. R1-R2 on the same chain → "head_to_tail"
    3. One side is a CHEM chain → "custom"
    4. Default → "peptide"
    """
    from_is_chem = from_chain.startswith("CHEM")
    to_is_chem = to_chain.startswith("CHEM")

    rgs = {from_rgroup, to_rgroup}

    if rgs == {"R3", "R3"}:
        if from_is_chem or to_is_chem:
            return "thioether"
        return "disulfide"

    if from_chain == to_chain and rgs == {"R1", "R2"}:
        return "head_to_tail"

    if from_is_chem or to_is_chem:
        return "custom"

    return "peptide"


# ---------------------------------------------------------------------------
# Topology classification
# ---------------------------------------------------------------------------

def _classify_topology(
    chains: list[dict],
    connections: list[dict],
) -> str:
    """
    Infer topology_class from the parsed chains and connections.

    Priority order:
    - Multiple PEPTIDE chains (even without connections) → "multichain"
    - No connections → "linear"
    - Multiple cyclic features → "bicyclic"
    - R3-R3 intrachain → "monocyclic_disulfide"
    - R1-R2 same chain → "monocyclic_backbone"
    - R1-R3 or R3-R2 same chain → "lariat" (side-chain-to-terminus ring)
    - Connection to CHEM chain → "branched_conjugate"
    - Multiple PEPTIDE chains with interchain connections → "multichain"
    - Default → "linear"
    """
    peptide_chains = [c for c in chains if c.get("type") == "PEPTIDE"]
    if len(peptide_chains) > 1 and not connections:
        return "multichain"

    if not connections:
        return "linear"

    has_disulfide_loop = False
    has_backbone_loop = False
    has_lariat = False
    has_chem_connection = False
    has_interchain = False

    for conn in connections:
        fc = conn["from_chain"]
        tc = conn["to_chain"]
        frg = conn["from_rgroup"]
        trg = conn["to_rgroup"]

        if fc.startswith("CHEM") or tc.startswith("CHEM"):
            has_chem_connection = True
            continue

        if fc != tc:
            has_interchain = True
            continue

        # Same chain
        rgs = {frg, trg}
        if rgs == {"R3", "R3"}:
            if not has_disulfide_loop:
                has_disulfide_loop = True
            else:
                return "bicyclic"
        elif rgs == {"R1", "R2"}:
            if not has_backbone_loop:
                has_backbone_loop = True
            else:
                return "bicyclic"
        elif rgs in ({"R1", "R3"}, {"R2", "R3"}):
            # Side-chain terminus connection → lariat
            has_lariat = True

    # Two different cyclic types → bicyclic
    if has_disulfide_loop and has_backbone_loop:
        return "bicyclic"

    if has_disulfide_loop:
        return "monocyclic_disulfide"

    if has_backbone_loop:
        return "monocyclic_backbone"

    if has_lariat:
        return "lariat"

    if has_chem_connection:
        return "branched_conjugate"

    if has_interchain:
        return "multichain"

    return "linear"


# ---------------------------------------------------------------------------
# HELMParser
# ---------------------------------------------------------------------------

class HELMParser:
    """
    Parse a HELM V2 notation string into a HELMObject.
    """

    @classmethod
    def parse(cls, helm_string: str) -> HELMObject:
        """
        Parse a HELM string and return a HELMObject.

        Parameters
        ----------
        helm_string : str
            Full HELM V2 string, e.g.::

                PEPTIDE1{A.K.C}$CHEM1{[Triazole14]}$PEPTIDE1,CHEM1,3:R3-1:R1$$V2.0

        Returns
        -------
        HELMObject
        """
        db = _get_db()
        sections = helm_string.split("$")

        # ---- Section 0: polymer blocks ----
        chains: list[dict] = cls._parse_polymer_blocks(sections[0], db)

        # ---- Section 1 or 2: connections ----
        # HELM V2 can have polymer blocks in section 0 (|‑separated) or
        # split across sections 0 and 1 ($-separated) for multi-polymer strings.
        # Connections land in the first section that contains a connection
        # pattern (from_chain,to_chain,pos:rg-pos:rg).  We scan sections 1
        # and 2 to handle both conventions.
        connections: list[dict] = []
        for sec_idx in (1, 2):
            if len(sections) > sec_idx and sections[sec_idx].strip():
                candidate = sections[sec_idx].strip()
                # If the candidate looks like a polymer block (contains '{'),
                # it is another polymer chain — fold it into chains and skip.
                if "{" in candidate and _CHAIN_RE.search(candidate):
                    extra_chains = cls._parse_polymer_blocks(candidate, db)
                    chains.extend(extra_chains)
                    continue
                # Otherwise treat as connections section
                connections = cls._parse_connections(candidate, chains)
                break

        # ---- Topology classification ----
        topology_class = _classify_topology(chains, connections)

        # ---- Build primary_string ----
        # Use the longest PEPTIDE chain (or first chain) as primary
        peptide_chains = [c for c in chains if c["type"] == "PEPTIDE"]
        if peptide_chains:
            primary_chain = max(peptide_chains, key=lambda c: len(c["monomers"]))
        else:
            primary_chain = chains[0] if chains else {"monomers": []}

        primary_string = ".".join(m["symbol"] for m in primary_chain["monomers"])

        # ---- Build global monomer index offset map ----
        # Maps chain_id → start index in a flat global sequence
        offset_map: dict[str, int] = {}
        cursor = 0
        for chain in chains:
            offset_map[chain["chain_id"]] = cursor
            cursor += len(chain["monomers"])

        # ---- Build connectivity_map (global indices) ----
        connectivity_map: list[tuple[int, int, str]] = []
        for conn in connections:
            fc = conn["from_chain"]
            tc = conn["to_chain"]
            fp = conn["from_pos"] - 1  # 1-based → 0-based
            tp = conn["to_pos"] - 1

            # Validate positions are within chain bounds
            from_chain_len = next((len(c["monomers"]) for c in chains if c["chain_id"] == fc), 0)
            to_chain_len = next((len(c["monomers"]) for c in chains if c["chain_id"] == tc), 0)
            if from_chain_len > 0 and fp >= from_chain_len:
                raise ValueError(f"Connection from_pos {conn['from_pos']} exceeds chain {fc} length {from_chain_len}")
            if to_chain_len > 0 and tp >= to_chain_len:
                raise ValueError(f"Connection to_pos {conn['to_pos']} exceeds chain {tc} length {to_chain_len}")

            u = offset_map.get(fc, 0) + fp
            v = offset_map.get(tc, 0) + tp
            connectivity_map.append((u, v, conn["bond_type"]))

        # ---- Build connectivity_graph (schema-compliant list of dicts) ----
        connectivity_graph = [
            {
                "from_chain": conn["from_chain"],
                "from_pos": conn["from_pos"],
                "from_rgroup": conn["from_rgroup"],
                "to_chain": conn["to_chain"],
                "to_pos": conn["to_pos"],
                "to_rgroup": conn["to_rgroup"],
                "bond_type": conn["bond_type"],
            }
            for conn in connections
        ]

        # ---- Assemble HELMObject payload ----
        helm_json = {
            "metadata": {
                "name": helm_string[:50],
            },
            "sequence": {
                "primary_string": primary_string,
            },
            "rdkit_data": {},
            "helm_string": helm_string,
            "connectivity_map": connectivity_map,
            "connectivity_graph": connectivity_graph,
            "topology_class": topology_class,
            # Full chain data for multi-chain peptides
            "_chains": chains,
        }

        return HELMObject(helm_json)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_polymer_blocks(section: str, db: MonomerDB) -> list[dict]:
        """
        Parse section 0 of a HELM string into a list of chain dicts.

        Each chain dict::

            {
                "chain_id": str,
                "type": "PEPTIDE" | "CHEM" | "RNA" | "DNA" | "UNKNOWN",
                "monomers": [
                    {"pos": int, "symbol": str, "entry": dict | None},
                    ...
                ]
            }
        """
        chains: list[dict] = []
        for block in section.split("|"):
            block = block.strip()
            if not block:
                continue
            m = _CHAIN_RE.match(block)
            if not m:
                continue

            chain_id, monomer_string = m.group(1).upper(), m.group(2)

            # Determine polymer type from chain_id prefix
            if chain_id.startswith("PEPTIDE"):
                chain_type = "PEPTIDE"
            elif chain_id.startswith("CHEM"):
                chain_type = "CHEM"
            elif chain_id.startswith("RNA"):
                chain_type = "RNA"
            elif chain_id.startswith("DNA"):
                chain_type = "DNA"
            else:
                chain_type = "UNKNOWN"

            monomers: list[dict] = []
            for pos, raw_sym in enumerate(monomer_string.split("."), start=1):
                # Strip surrounding square brackets if present (non-standard monomers)
                symbol = raw_sym.strip()
                if symbol.startswith("[") and symbol.endswith("]"):
                    symbol = symbol[1:-1].strip()
                symbol = _normalize_symbol(symbol)

                entry = db.find_by_symbol(symbol)
                monomers.append({
                    "pos": pos,
                    "symbol": symbol,
                    "entry": entry,
                })

            chains.append({
                "chain_id": chain_id,
                "type": chain_type,
                "monomers": monomers,
            })

        return chains

    @staticmethod
    def _resolve_pos(chain_id: str, rgroup: str, explicit: str | None, chains: list[dict]) -> int:
        """Return 1-based position when explicit position string may be absent."""
        if explicit is not None:
            return int(explicit)
        chain = next((c for c in chains if c['chain_id'] == chain_id), None)
        n = len(chain['monomers']) if chain else 1
        return n if rgroup == 'R2' else 1

    @staticmethod
    def _parse_connections(section: str, chains: list[dict]) -> list[dict]:
        """
        Parse section 1 of a HELM string into a list of connection dicts.
        Handles two formats:
          Standard:  CHAIN1,CHAIN2,pos1:Rg1-pos2:Rg2
          Alternate: CHAIN1:Rg1(pos1)-CHAIN2:Rg2(pos2)  (pos optional)

        Each connection dict::

            {
                "from_chain": str,
                "from_pos": int,   # 1-based
                "from_rgroup": str,
                "to_chain": str,
                "to_pos": int,     # 1-based
                "to_rgroup": str,
                "bond_type": str,
            }
        """
        connections: list[dict] = []
        for conn_str in section.split("|"):
            conn_str = conn_str.strip()
            if not conn_str:
                continue

            m = _CONN_RE.match(conn_str)
            if m:
                from_chain  = m.group(1).upper()
                to_chain    = m.group(2).upper()
                from_pos    = int(m.group(3))
                from_rgroup = m.group(4)
                to_pos      = int(m.group(5))
                to_rgroup   = m.group(6)
            else:
                m = _CONN_RE_ALT.match(conn_str)
                if not m:
                    continue
                from_chain  = m.group(1).upper()
                from_rgroup = m.group(2)
                to_chain    = m.group(4).upper()
                to_rgroup   = m.group(5)
                from_pos    = HELMParser._resolve_pos(from_chain, from_rgroup, m.group(3), chains)
                to_pos      = HELMParser._resolve_pos(to_chain,   to_rgroup,   m.group(6), chains)

            bond_type = _infer_bond_type(from_chain, from_rgroup, to_chain, to_rgroup)

            connections.append({
                "from_chain": from_chain,
                "from_pos":   from_pos,
                "from_rgroup": from_rgroup,
                "to_chain":   to_chain,
                "to_pos":     to_pos,
                "to_rgroup":  to_rgroup,
                "bond_type":  bond_type,
            })

        return connections
