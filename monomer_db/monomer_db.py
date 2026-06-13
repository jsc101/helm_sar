"""
MonomerDB — unified lookup for HELM monomers.

Loads HELMCoreLibrary.json (from PistoiaHELM/HELMMonomerSets) and
custom_monomers.json from the same directory, then builds three indexes:

    _by_symbol      : {symbol: entry}
    _index          : {(canonical_smiles, n_rgroups): entry}  stereo-aware
    _index_nostereo : {(stereo_free_smiles, n_rgroups): entry}  fallback
"""

from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_HERE = Path(__file__).parent


# ---------------------------------------------------------------------------
# SMILES normalisation helpers
# ---------------------------------------------------------------------------

def _helm_to_star_smiles(smiles: str) -> str:
    """
    Convert HELM attachment-point notation to RDKit dummy-atom notation.

    HELM uses atom-mapped caps inline, e.g.::

        C[C@H](N[H:1])C([OH:2])=O   ->  C[C@H](N[*:1])C([*:2])=O

    Patterns handled:
        [H:N]   -> [*:N]   (H-cap attachment)
        [OH:N]  -> [*:N]   (OH-cap attachment, discard O for indexing)
        [*:N]   unchanged  (already in star form, used by custom_monomers)
    """
    s = re.sub(r'\[H:(\d+)\]', r'[*:\1]', smiles)
    s = re.sub(r'\[OH:(\d+)\]', r'[*:\1]', s)
    return s


def _canonicalise(smiles: str, stereo: bool = True) -> Optional[tuple[str, int]]:
    """
    Return (canonical_smiles, n_rgroups) or None if RDKit cannot parse.

    n_rgroups is the count of atoms carrying a non-zero atom-map number,
    which corresponds to HELM attachment points.
    """
    try:
        from rdkit import Chem
    except ImportError:
        logger.warning("RDKit not available — SMILES index will be empty")
        return None

    star_smiles = _helm_to_star_smiles(smiles)
    mol = Chem.MolFromSmiles(star_smiles)
    if mol is None:
        return None

    n_rgroups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0)

    if not stereo:
        Chem.RemoveStereochemistry(mol)
        # Normalise tautomers in the stereo-free index so that, e.g., the two
        # imidazole tautomers of His (c1cnc[nH]1 vs c1c[nH]cn1) share one key.
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        except Exception:
            pass  # fall back to non-normalised form if unavailable

    canon = Chem.MolToSmiles(mol)
    return canon, n_rgroups


# ---------------------------------------------------------------------------
# Schema normalisation helpers
# ---------------------------------------------------------------------------

def _normalise_rgroups(raw_rgroups: list[dict]) -> list[dict]:
    """
    Return a normalised rgroups list where every entry has at least:
        label, capGroupSmiles
    Both HELMCoreLibrary and custom_monomers are handled.
    """
    out = []
    for r in raw_rgroups:
        label = r.get("label", "")
        # HELMCoreLibrary uses capGroupSmiles (lowercase s)
        # custom_monomers uses capGroupSmiles as well
        # monomerLib2.0 uses capGroupSMILES (uppercase)
        cap = (
            r.get("capGroupSmiles")
            or r.get("capGroupSMILES")
            or r.get("capGroupSmiles", "")
        )
        norm = {"label": label, "capGroupSmiles": cap}
        # pass through any extra keys unchanged
        for k, v in r.items():
            if k not in norm:
                norm[k] = v
        out.append(norm)
    return out


def _normalise_entry(raw: dict) -> dict:
    """
    Return a normalised monomer entry with consistent field names:
        symbol, name, polymerType, monomerType, smiles, rgroups
    """
    symbol = raw.get("symbol") or raw.get("id", "")
    name = raw.get("name", "")
    polymer_type = raw.get("polymerType", "")
    monomer_type = raw.get("monomerType", "")
    smiles = raw.get("smiles", "")
    rgroups = _normalise_rgroups(raw.get("rgroups", []))

    entry = dict(raw)  # keep all original fields
    entry["symbol"] = symbol
    entry["name"] = name
    entry["polymerType"] = polymer_type
    entry["monomerType"] = monomer_type
    entry["smiles"] = smiles
    entry["rgroups"] = rgroups
    return entry


# ---------------------------------------------------------------------------
# MonomerDB
# ---------------------------------------------------------------------------

class MonomerDB:
    """
    Unified monomer lookup for HELM-notation peptide informatics.

    Attributes
    ----------
    _by_symbol : dict[str, dict]
        Direct lookup by monomer symbol (e.g. "A", "Aha", "Triazole14").
    _index : dict[tuple[str, int], dict]
        Stereo-aware index keyed by (canonical_smiles, n_rgroups).
    _index_nostereo : dict[tuple[str, int], dict]
        Stereo-free fallback index keyed by (canonical_smiles, n_rgroups).
    """

    def __init__(
        self,
        core_path: Optional[str | Path] = None,
        custom_path: Optional[str | Path] = None,
        extra_sources: Optional[list[str]] = None,
    ) -> None:
        core_path = Path(core_path) if core_path else _HERE / "HELMCoreLibrary.json"
        custom_path = Path(custom_path) if custom_path else _HERE / "custom_monomers.json"

        self._by_symbol: dict[str, dict] = {}
        self._index: dict[tuple[str, int], dict] = {}
        self._index_nostereo: dict[tuple[str, int], dict] = {}

        raw_entries: list[dict] = []
        raw_entries.extend(self._load_json(core_path))
        raw_entries.extend(self._load_json(custom_path))

        # Load from extra_sources before processing
        if extra_sources:
            for path in extra_sources:
                path = os.path.expanduser(path)
                if not os.path.exists(path):
                    logger.warning("extra_sources path not found: %s", path)
                    continue
                try:
                    raw_entries.extend(self._load_json(Path(path)))
                except Exception as e:
                    logger.warning("Failed to load extra_sources %s: %s", path, e)

        self._load_entries(raw_entries)
        logger.info(
            "MonomerDB loaded: %d symbols, %d stereo-smiles, %d nostereo-smiles",
            len(self._by_symbol),
            len(self._index),
            len(self._index_nostereo),
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _load_json(path: Path) -> list[dict]:
        if not path.exists():
            logger.warning("Monomer file not found: %s", path)
            return []
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        if isinstance(data, list):
            return data
        # Some formats wrap in a dict
        if isinstance(data, dict):
            for key in ("monomers", "monomerList", "items"):
                if key in data:
                    return data[key]
            logger.warning("Unknown JSON structure in %s — skipping", path)
            return []
        return []

    def _load_entries(self, raw_entries: list[dict]) -> None:
        """
        Process a list of raw monomer entries, adding them to the internal indexes.
        Handles duplicate symbol resolution via polymer-type priority.
        """
        # Polymer-type priority for duplicate symbol resolution.
        # Higher value = higher priority.  PEPTIDE wins over RNA for shared
        # single-letter codes (A, C, G, T).
        _PRIORITY = {"PEPTIDE": 3, "CHEM": 2, "RNA": 1, "DNA": 1}

        for raw in raw_entries:
            entry = _normalise_entry(raw)
            symbol = entry["symbol"]

            if symbol in self._by_symbol:
                existing = self._by_symbol[symbol]
                existing_pri = _PRIORITY.get(existing.get("polymerType", ""), 0)
                new_pri = _PRIORITY.get(entry.get("polymerType", ""), 0)
                if new_pri <= existing_pri:
                    # Keep existing higher-priority entry; still index the SMILES
                    logger.debug(
                        "Duplicate symbol %r: keeping %s over %s",
                        symbol,
                        existing.get("polymerType"),
                        entry.get("polymerType"),
                    )
                else:
                    logger.debug(
                        "Duplicate symbol %r: replacing %s with higher-priority %s",
                        symbol,
                        existing.get("polymerType"),
                        entry.get("polymerType"),
                    )
                    self._by_symbol[symbol] = entry
            else:
                self._by_symbol[symbol] = entry

            smiles = entry.get("smiles", "")
            if not smiles:
                continue

            # Stereo-aware index
            result = _canonicalise(smiles, stereo=True)
            if result is None:
                logger.warning(
                    "Could not parse SMILES for %r (%r) — skipping SMILES index entry",
                    symbol,
                    smiles,
                )
            else:
                key = result
                if key not in self._index:
                    self._index[key] = entry
                else:
                    logger.warning(
                        "SMILES collision in stereo index: key %r already maps to %r, ignoring %r",
                        key, self._index[key].get("symbol"), entry.get("symbol")
                    )

            # Stereo-free index
            result_ns = _canonicalise(smiles, stereo=False)
            if result_ns is not None:
                key_ns = result_ns
                if key_ns not in self._index_nostereo:
                    self._index_nostereo[key_ns] = entry
                else:
                    logger.warning(
                        "SMILES collision in nostereo index: key %r already maps to %r, ignoring %r",
                        key_ns, self._index_nostereo[key_ns].get("symbol"), entry.get("symbol")
                    )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def find(self, smiles: str, stereo: bool = True) -> Optional[dict]:
        """
        Look up a monomer by SMILES.

        Tries the stereo-aware index first; if ``stereo=True`` and not found,
        falls back to the stereo-free index automatically.  Set ``stereo=False``
        to search only the stereo-free index.

        Parameters
        ----------
        smiles:
            Query SMILES (any format accepted by RDKit, including HELM notation
            with ``[H:N]`` / ``[OH:N]`` attachment points).
        stereo:
            When True (default) try stereo match first, then nostereo fallback.
            When False, skip stereo match and go straight to nostereo.

        Returns
        -------
        dict or None
        """
        if stereo:
            result = _canonicalise(smiles, stereo=True)
            if result is not None and result in self._index:
                return self._index[result]
        # fallback
        result_ns = _canonicalise(smiles, stereo=False)
        if result_ns is not None:
            return self._index_nostereo.get(result_ns)
        return None

    def find_by_symbol(self, symbol: str) -> Optional[dict]:
        """
        Look up a monomer by its HELM symbol.

        Parameters
        ----------
        symbol : str
            HELM symbol, e.g. ``"A"`` for Alanine, ``"Aha"``, ``"Triazole14"``.

        Returns
        -------
        dict or None
        """
        hit = self._by_symbol.get(symbol)
        if hit is not None:
            return hit
        # HELMParser converts [D-Dab] → D_Dab; try the dash form as fallback
        if '_' in symbol:
            dash_form = symbol.replace('_', '-')
            hit = self._by_symbol.get(dash_form)
            if hit is not None:
                return hit
            # Also try partial: only first underscore (D_Dab → D-Dab but Ala_tBu stays)
            first = symbol.replace('_', '-', 1)
            if first != dash_form:
                hit = self._by_symbol.get(first)
                if hit is not None:
                    return hit
        return None

    def get_rgroups(self, symbol: str) -> list[dict]:
        """
        Return the rgroups list for a monomer symbol.

        Parameters
        ----------
        symbol : str

        Returns
        -------
        list[dict]
            List of R-group dicts (normalised to always contain ``label``
            and ``capGroupSmiles`` keys).  Empty list if symbol not found.
        """
        entry = self._by_symbol.get(symbol)
        if entry is None:
            return []
        return entry.get("rgroups", [])

    def register(self, entry: dict) -> None:
        """Add a monomer entry to the live in-memory indices.

        Intended for auto-discovered monomers during SMILES→HELM fragmentation.
        The entry must have at minimum a ``symbol`` and ``smiles`` key; ``rgroups``
        defaults to R1+R2 if absent.  Changes are not persisted to disk.
        """
        entry = dict(entry)
        entry.setdefault('polymerType', 'PEPTIDE')
        entry.setdefault('monomerType', 'Backbone')
        entry.setdefault('rgroups', [
            {'label': 'R1', 'capGroupSmiles': '[H]',  'description': 'N-terminus'},
            {'label': 'R2', 'capGroupSmiles': '[OH]', 'description': 'C-terminus'},
        ])
        self._by_symbol[entry['symbol']] = entry
        smiles = entry.get('smiles', '')
        if smiles:
            key = _canonicalise(smiles, stereo=True)
            if key and key not in self._index:
                self._index[key] = entry
            key_ns = _canonicalise(smiles, stereo=False)
            if key_ns and key_ns not in self._index_nostereo:
                self._index_nostereo[key_ns] = entry

    def __len__(self) -> int:
        return len(self._by_symbol)

    def __repr__(self) -> str:
        return (
            f"MonomerDB(symbols={len(self._by_symbol)}, "
            f"smiles_entries={len(self._index)})"
        )
