# Peptide Design Excel→HELM Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a tool that lets a chemist author peptide designs in a flat one-row-per-compound Excel format (main-chain columns + alternating bond/monomer sidechain columns) and generate valid HELM V2 strings from it.

**Architecture:** A Python module (`peptide_design/`) parses the flat row format into structured dataclasses, then generates HELM using a collapse-or-expand rule (collapse to a single chain when all bonds are default head-to-tail within one polymer type; otherwise emit one chain per monomer with explicit bond lines). An Excel `.xlsm` template (populated by `make_template.py` + a `helm_builder.bas` VBA macro) provides the standalone ready-to-go deliverable using the same algorithm.

**Tech Stack:** Python 3.10+, openpyxl (Excel template), existing `monomer_db/monomer_db.py` (MonomerDB class), existing `scripts/helm_parser.py` (validator). VBA for the Excel macro.

---

## Column naming convention

Main chain: columns with headers that are digit-only strings (`"1"`, `"2"`, … `"N"`), sorted by integer value — these are HELM positions (1-based).

Sidechain (repeating up to 3 times): `Site`, `b1`, `SC1`, `b2`, `SC2`, … for sidechain 1; `Site_2`, `b1_2`, `SC1_2`, … for sidechain 2; `Site_3`, `b1_3`, `SC1_3`, … for sidechain 3.

- `Site` — integer HELM position where the sidechain attaches (must match the column header of the attachment residue)
- `bk` — bond entering `SCk`, written `proxRg-distRg` (e.g. `R3-R1`). proxRg = R-group from the previous piece (backbone for k=1, prior monomer for k>1); distRg = R-group on `SCk` that faces back.
- `SCk` — monomer symbol (proximal→distal, reading outward from backbone). Last monomer has no trailing bond — its distal terminus is free.

### Collapse rule

A sidechain collapses to a single-chain HELM notation when **all** of:
1. Every monomer has `polymerType == "PEPTIDE"`
2. `b1`'s distal Rg (right side of `-`) == `"R1"`
3. Every `bk` (k ≥ 2) == `"R2-R1"`

Collapsed output: `PEPTIDE{N}{[SC1].[SC2]…[SCK]}` + one connection line `PEPTIDE{N},PEPTIDE1,1:R1-{site}:{b1_prox}`.

### Expand rule

If collapse fails: emit one chain per monomer (PEPTIDE or CHEM per DB lookup), plus explicit connection lines. Bond bk = `proxRg-distRg` → connection line: `chain_k,chain_{k-1},1:{distRg}-1:{proxRg}`. For b1: `chain_1,PEPTIDE1,1:{distRg_b1}-{site}:{proxRg_b1}`.

---

## File structure

```
peptide_design/
├── __init__.py
├── core.py          # Compound + Sidechain dataclasses, parse_row()
├── generator.py     # generate_helm(), collapse/expand logic, fmt_sym()
├── validator.py     # validate_helm() → (ok: bool, message: str)
├── cli.py           # reads xlsx/csv → writes HELM + status column
└── make_template.py # creates .xlsx template (Design/HELM/MonomerDB/Example sheets)
helm_builder.bas     # VBA macro (lives at helm_sar/ root for easy import)
tests/
└── peptide_design/
    ├── test_core.py
    ├── test_generator.py
    └── test_validator.py
```

MonomerDB is loaded from `monomer_db/monomer_db.py` (already exists) using its `MonomerDB` class.

---

## Task 1: Scaffold + core data model

**Files:**
- Create: `peptide_design/__init__.py`
- Create: `peptide_design/core.py`
- Create: `tests/peptide_design/__init__.py`
- Create: `tests/peptide_design/test_core.py`

- [ ] **Step 1: Write failing tests for parse_row**

```python
# tests/peptide_design/test_core.py
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import Compound, Sidechain, parse_row

def _pyy_row():
    """Analog_28-pos30 from PYY library: C18 at HELM pos 30."""
    row = {"Name": "Analog_28-pos30"}
    # 34-residue PYY3-36 backbone (positions 1..34)
    seq = list("IKPEAPGEDASPEE") + list("LNRYYASLRHYLNL") + list("VTRQRY")
    assert len(seq) == 34
    for i, aa in enumerate(seq, 1):
        row[str(i)] = aa
    # Put K at position 30
    row["30"] = "K"
    row["Site"] = "30"
    row["b1"] = "R3-R1"
    row["SC1"] = "gGlu"
    row["b2"] = "R2-R1"
    row["SC2"] = "Ado"
    row["b3"] = "R2-R1"
    row["SC3"] = "Ado"
    row["b4"] = "R2-R1"
    row["SC4"] = "C18d"
    return row

def test_parse_row_name():
    c = parse_row(_pyy_row())
    assert c.name == "Analog_28-pos30"

def test_parse_row_main_chain_length():
    c = parse_row(_pyy_row())
    assert len(c.main_chain) == 34

def test_parse_row_main_chain_k_at_30():
    c = parse_row(_pyy_row())
    assert c.main_chain[29] == "K"  # 0-indexed

def test_parse_row_one_sidechain():
    c = parse_row(_pyy_row())
    assert len(c.sidechains) == 1

def test_parse_row_sidechain_site():
    c = parse_row(_pyy_row())
    assert c.sidechains[0].site == 30

def test_parse_row_sidechain_bonds():
    c = parse_row(_pyy_row())
    sc = c.sidechains[0]
    assert sc.bonds == ["R3-R1", "R2-R1", "R2-R1", "R2-R1"]

def test_parse_row_sidechain_monomers():
    c = parse_row(_pyy_row())
    sc = c.sidechains[0]
    assert sc.monomers == ["gGlu", "Ado", "Ado", "C18d"]

def test_parse_row_no_sidechain():
    row = {"Name": "ref"}
    for i, aa in enumerate(list("IKPEAPGE"), 1):
        row[str(i)] = aa
    c = parse_row(row)
    assert c.sidechains == []

def test_parse_row_two_sidechains():
    row = {"Name": "bis"}
    for i, aa in enumerate(list("IKPEAPGE"), 1):
        row[str(i)] = aa
    row.update({"Site": "3", "b1": "R3-R1", "SC1": "Ado", "b2": "R2-R1", "SC2": "C18d"})
    row.update({"Site_2": "6", "b1_2": "R3-R1", "SC1_2": "Ado", "b2_2": "R2-R1", "SC2_2": "C12d"})
    c = parse_row(row)
    assert len(c.sidechains) == 2
    assert c.sidechains[1].site == 6
    assert c.sidechains[1].monomers == ["Ado", "C12d"]
```

- [ ] **Step 2: Run tests — expect failure (module missing)**

```bash
cd /Users/claudedev/Projects/helm_sar
python -m pytest tests/peptide_design/test_core.py -v 2>&1 | head -20
```

Expected: `ModuleNotFoundError: No module named 'peptide_design'`

- [ ] **Step 3: Create scaffold and core.py**

```python
# peptide_design/__init__.py
# (empty)
```

```python
# peptide_design/core.py
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
    # main chain: all columns with digit-only headers, sorted by integer value
    mc_cols = sorted([k for k in row if str(k).isdigit()], key=int)
    main_chain = [str(row[c]).strip() for c in mc_cols if str(row.get(c, "")).strip()]
    sidechains = []
    for suffix in ("", "_2", "_3"):
        site_key = f"Site{suffix}"
        raw = str(row.get(site_key, "")).strip()
        if not raw:
            break
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
```

- [ ] **Step 4: Create test __init__ and run tests**

```bash
touch tests/__init__.py tests/peptide_design/__init__.py
python -m pytest tests/peptide_design/test_core.py -v
```

Expected: all 9 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add peptide_design/__init__.py peptide_design/core.py \
        tests/__init__.py tests/peptide_design/__init__.py \
        tests/peptide_design/test_core.py
git commit -m "feat: peptide_design core data model and row parser"
```

---

## Task 2: HELM generator — collapse path

**Files:**
- Create: `peptide_design/generator.py`
- Create: `tests/peptide_design/test_generator.py` (collapse cases)

- [ ] **Step 1: Write failing tests for collapse path**

```python
# tests/peptide_design/test_generator.py
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import Compound, Sidechain, parse_row
from peptide_design.generator import generate_helm

# Helper
def _simple_compound(backbone, site, bonds, monomers):
    sc = Sidechain(site=site, bonds=bonds, monomers=monomers)
    return Compound(name="test", main_chain=list(backbone), sidechains=[sc])

# ── Collapse tests ──────────────────────────────────────────────────────────

def test_no_sidechain():
    c = Compound(name="ref", main_chain=["I","K","P","E"])
    helm = generate_helm(c)
    assert helm == "PEPTIDE1{I.K.P.E}$$$$V2.0"

def test_collapse_c18_at_pos4():
    """Plain gGlu-Ado-Ado-C18d all PEPTIDE, all R2-R1 after b1."""
    c = _simple_compound(
        ["I","K","P","K"],
        site=4,
        bonds=["R3-R1","R2-R1","R2-R1","R2-R1"],
        monomers=["gGlu","Ado","Ado","C18d"],
    )
    helm = generate_helm(c)
    assert "PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}" in helm
    assert "PEPTIDE2,PEPTIDE1,1:R1-4:R3" in helm
    assert helm.count("PEPTIDE") == 2  # only PEPTIDE1 + PEPTIDE2

def test_collapse_n_terminal_r1():
    """Attachment via R1 (N-terminal, not epsilon-amine)."""
    c = _simple_compound(
        ["I","K","P","E"],
        site=1,
        bonds=["R1-R1","R2-R1","R2-R1","R2-R1"],
        monomers=["gGlu","Ado","Ado","C18d"],
    )
    helm = generate_helm(c)
    # backbone_prox = R1 (left of b1 = "R1-R1")
    assert "PEPTIDE2,PEPTIDE1,1:R1-1:R1" in helm

def test_collapse_single_monomer():
    c = _simple_compound(["A","K"], site=2,
                         bonds=["R3-R1"], monomers=["C18d"])
    helm = generate_helm(c)
    assert "PEPTIDE2{[C18d]}" in helm
    assert "PEPTIDE2,PEPTIDE1,1:R1-2:R3" in helm

def test_fmt_multi_char_monomer_gets_brackets():
    c = _simple_compound(["A","K"], site=2,
                         bonds=["R3-R1","R2-R1"], monomers=["Ado","C18d"])
    helm = generate_helm(c)
    assert "[Ado]" in helm
    assert "[C18d]" in helm

def test_fmt_single_char_monomer_no_brackets():
    c = _simple_compound(["A","K","G"], site=3,
                         bonds=["R3-R1"], monomers=["G"])
    helm = generate_helm(c)
    assert "PEPTIDE2{G}" in helm
```

- [ ] **Step 2: Run — expect failure**

```bash
python -m pytest tests/peptide_design/test_generator.py -v 2>&1 | head -15
```

Expected: `ModuleNotFoundError: No module named 'peptide_design.generator'`

- [ ] **Step 3: Implement generator.py (collapse path)**

```python
# peptide_design/generator.py
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from peptide_design.core import Compound, Sidechain

# Lazy-load MonomerDB to avoid import cost when not needed
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
    """True when sidechain can be expressed as a single PEPTIDE chain."""
    if not sc.bonds or not sc.monomers:
        return False
    # All monomers must be PEPTIDE
    if any(_polymer_type(m) != "PEPTIDE" for m in sc.monomers):
        return False
    # b1 distal Rg (right of '-') must be R1
    if sc.bonds[0].split("-")[1] != "R1":
        return False
    # All subsequent bonds must be R2-R1
    if any(b != "R2-R1" for b in sc.bonds[1:]):
        return False
    return True


def _build_collapsed(sc: Sidechain, chain_n: int) -> tuple[str, str]:
    """Return (chain_def, connection_line) for a collapsible sidechain."""
    body = ".".join(fmt_sym(m) for m in sc.monomers)
    chain_id = f"PEPTIDE{chain_n}"
    chain_def = f"{chain_id}{{{body}}}"
    prox_rg = sc.bonds[0].split("-")[0]  # left side of b1
    conn = f"{chain_id},PEPTIDE1,1:R1-{sc.site}:{prox_rg}"
    return chain_def, conn


def _build_expanded(sc: Sidechain, peptide_n: int, chem_n: int
                    ) -> tuple[list[str], list[str], int, int]:
    """
    Return (chain_defs, connections, next_peptide_n, next_chem_n).
    Emits one chain per monomer.
    """
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
            # Connect to backbone PEPTIDE1
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
```

- [ ] **Step 4: Run collapse tests**

```bash
python -m pytest tests/peptide_design/test_generator.py -v
```

Expected: all 6 collapse tests PASS.

- [ ] **Step 5: Commit**

```bash
git add peptide_design/generator.py tests/peptide_design/test_generator.py
git commit -m "feat: HELM generator collapse path"
```

---

## Task 3: HELM generator — expand path (click chemistry)

**Files:**
- Modify: `tests/peptide_design/test_generator.py` (add expand tests)

- [ ] **Step 1: Add expand tests**

Append to `tests/peptide_design/test_generator.py`:

```python
# ── Expand tests ─────────────────────────────────────────────────────────────

def test_expand_triggered_by_chem_monomer():
    """Triazole14 is CHEM — forces expand path."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    assert "CHEM1{[Triazole14]}" in helm

def test_expand_backbone_connection_uses_b1():
    """PEPTIDE2 (Hpg) connects to PEPTIDE1 at site 37 via b1=R3-R1."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    # Hpg is PEPTIDE2, attaches at backbone pos 37 R3, via Hpg R1
    assert "PEPTIDE2,PEPTIDE1,1:R1-37:R3" in helm

def test_expand_triazole_connects_to_hpg():
    """CHEM1 (Triazole14) connects to PEPTIDE2 (Hpg): b2=R3-R1."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    # b2=R3-R1: CHEM1 R1 ← Hpg R3
    assert "CHEM1,PEPTIDE2,1:R1-1:R3" in helm

def test_expand_aha_connects_to_triazole():
    """PEPTIDE3 (Aha) connects to CHEM1 (Triazole14): b3=R4-R3."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    # b3=R4-R3: PEPTIDE3 R3 ← CHEM1 R4
    assert "PEPTIDE3,CHEM1,1:R3-1:R4" in helm

def test_expand_non_default_bond_without_chem():
    """Non-default bond (R3-R1 on b2) forces expand even for all-PEPTIDE."""
    c = _simple_compound(
        ["A","K"],
        site=2,
        bonds=["R3-R1", "R3-R1"],   # b2 is not R2-R1 → expand
        monomers=["Hpg", "Aha"],
    )
    helm = generate_helm(c)
    # Must emit two separate PEPTIDE chains
    assert "PEPTIDE2{[Hpg]}" in helm
    assert "PEPTIDE3{[Aha]}" in helm
    assert "PEPTIDE3,PEPTIDE2,1:R1-1:R3" in helm
```

- [ ] **Step 2: Run — expect expand tests to fail**

```bash
python -m pytest tests/peptide_design/test_generator.py::test_expand_triggered_by_chem_monomer -v
```

Expected: FAIL — `AssertionError` (expand path not exercised yet, or `CHEM1` not emitted correctly).

- [ ] **Step 3: Verify the generator already handles expand (it was written in Task 2)**

The `_build_expanded` function and `_can_collapse` check are already in `generator.py`. Run all generator tests:

```bash
python -m pytest tests/peptide_design/test_generator.py -v
```

Expected: all tests PASS. If any expand test fails, check that `_polymer_type("Triazole14")` returns `"CHEM"` (requires MonomerDB to load correctly). Debug with:

```bash
python3 -c "
import sys; sys.path.insert(0,'.')
from monomer_db.monomer_db import MonomerDB
db = MonomerDB()
print(db.find_by_symbol('Triazole14'))
print(db.find_by_symbol('Hpg'))
"
```

- [ ] **Step 4: Commit**

```bash
git add tests/peptide_design/test_generator.py
git commit -m "test: expand path generator tests (click chemistry / non-default bonds)"
```

---

## Task 4: Validator

**Files:**
- Create: `peptide_design/validator.py`
- Create: `tests/peptide_design/test_validator.py`

- [ ] **Step 1: Write failing validator tests**

```python
# tests/peptide_design/test_validator.py
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import Compound, Sidechain
from peptide_design.generator import generate_helm
from peptide_design.validator import validate_helm

def _make_helm(backbone, site=None, bonds=None, monomers=None):
    scs = []
    if site is not None:
        scs = [Sidechain(site=site, bonds=bonds, monomers=monomers)]
    c = Compound(name="t", main_chain=list(backbone), sidechains=scs)
    return generate_helm(c)

def test_valid_no_sidechain():
    helm = _make_helm(["A","K","P","E"])
    ok, msg = validate_helm(helm)
    assert ok, msg

def test_valid_with_collapse_sidechain():
    helm = _make_helm(["A","K"], site=2,
                      bonds=["R3-R1","R2-R1"], monomers=["Ado","C18d"])
    ok, msg = validate_helm(helm)
    assert ok, msg

def test_invalid_helm_string():
    ok, msg = validate_helm("NOT_VALID_HELM")
    assert not ok
    assert msg  # has an error message

def test_unknown_monomer_flagged():
    # Generate HELM with a monomer not in any DB
    from peptide_design.core import Compound, Sidechain
    from peptide_design.generator import generate_helm
    c = Compound(name="t", main_chain=["A","K"],
                 sidechains=[Sidechain(site=2, bonds=["R3-R1"], monomers=["C16m"])])
    helm = generate_helm(c)
    ok, msg = validate_helm(helm)
    assert not ok
    assert "C16m" in msg or "unknown" in msg.lower()
```

- [ ] **Step 2: Run — expect failure**

```bash
python -m pytest tests/peptide_design/test_validator.py -v 2>&1 | head -10
```

Expected: `ModuleNotFoundError: No module named 'peptide_design.validator'`

- [ ] **Step 3: Implement validator.py**

```python
# peptide_design/validator.py
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

    db = _db()
    unknown = []
    for chain in obj.data.get("_chains", {}).values():
        for m in chain.get("monomers", []):
            sym = m.get("symbol", "")
            if sym and db.find_by_symbol(sym) is None:
                unknown.append(sym)

    if unknown:
        return False, f"Unknown monomers: {', '.join(sorted(set(unknown)))}"
    return True, "OK"
```

- [ ] **Step 4: Run validator tests**

```bash
python -m pytest tests/peptide_design/test_validator.py -v
```

Expected: all 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add peptide_design/validator.py tests/peptide_design/test_validator.py
git commit -m "feat: HELM validator (parse + monomer DB check)"
```

---

## Task 5: CLI

**Files:**
- Create: `peptide_design/cli.py`

- [ ] **Step 1: Implement cli.py**

```python
#!/usr/bin/env python3
# peptide_design/cli.py
"""
Read a CSV or XLSX in the peptide-design format and emit HELM + validation status.

Usage:
    python -m peptide_design.cli input.csv            # writes input_helm.csv
    python -m peptide_design.cli input.xlsx           # writes input_helm.xlsx
    python -m peptide_design.cli input.csv --stdout   # print to terminal
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

# ensure helm_sar root is on sys.path
sys.path.insert(0, str(Path(__file__).parent.parent))

from peptide_design.core import parse_row
from peptide_design.generator import generate_helm
from peptide_design.validator import validate_helm


def _read_csv(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    return reader.fieldnames or [], rows


def _read_xlsx(path: Path) -> tuple[list[str], list[dict]]:
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb.active
    rows_iter = ws.iter_rows(values_only=True)
    headers = [str(h) if h is not None else "" for h in next(rows_iter)]
    rows = []
    for row in rows_iter:
        if all(v is None for v in row):
            continue
        rows.append({h: (str(v).strip() if v is not None else "") for h, v in zip(headers, row)})
    return headers, rows


def process(rows: list[dict]) -> list[dict]:
    out = []
    for row in rows:
        compound = parse_row(row)
        helm = generate_helm(compound)
        ok, msg = validate_helm(helm)
        result = dict(row)
        result["HELM"] = helm
        result["HELM_status"] = "OK" if ok else f"ERROR: {msg}"
        out.append(result)
    return out


def _write_csv(path: Path, orig_headers: list[str], rows: list[dict]) -> None:
    extra = ["HELM", "HELM_status"]
    fieldnames = list(orig_headers) + [h for h in extra if h not in orig_headers]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description="Peptide design → HELM")
    p.add_argument("input", help="CSV or XLSX input file")
    p.add_argument("--stdout", action="store_true", help="Print HELM strings to stdout")
    args = p.parse_args(argv)

    path = Path(args.input)
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        headers, rows = _read_xlsx(path)
    else:
        headers, rows = _read_csv(path)

    results = process(rows)

    if args.stdout:
        for r in results:
            print(f"{r.get('Name','?')}\t{r['HELM']}\t{r['HELM_status']}")
        return

    out_path = path.with_stem(path.stem + "_helm").with_suffix(".csv")
    _write_csv(out_path, headers, results)
    print(f"Written {len(results)} rows → {out_path}")
    errors = [r for r in results if r["HELM_status"] != "OK"]
    if errors:
        print(f"  {len(errors)} validation errors:")
        for r in errors:
            print(f"    {r.get('Name','?')}: {r['HELM_status']}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test CLI on a hand-crafted CSV**

Create `/tmp/design_test.csv` manually:

```csv
Name,1,2,3,4,Site,b1,SC1,b2,SC2,b3,SC3,b4,SC4
ref,I,K,P,E,,,,,,,,,,
c18_at_4,I,K,P,K,4,R3-R1,gGlu,R2-R1,Ado,R2-R1,Ado,R2-R1,C18d
```

```bash
python -m peptide_design.cli /tmp/design_test.csv --stdout
```

Expected output (two lines):
```
ref     PEPTIDE1{I.K.P.E}$$$$V2.0      OK
c18_at_4    PEPTIDE1{I.K.P.K}|PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}$PEPTIDE2,PEPTIDE1,1:R1-4:R3$$$V2.0    OK
```

- [ ] **Step 3: Commit**

```bash
git add peptide_design/cli.py
git commit -m "feat: peptide_design CLI (csv/xlsx → HELM + validation)"
```

---

## Task 6: Integration test on PYY library

Prove the full pipeline on a real dataset. We re-encode the PYY CSV into the design format and verify the generated HELMs round-trip through HELMParser.

**Files:**
- Create: `tests/peptide_design/test_integration_pyy.py`

- [ ] **Step 1: Write integration test**

```python
# tests/peptide_design/test_integration_pyy.py
"""
Re-encode four representative PYY analogs into design-format rows and
verify that generate_helm + validate_helm produce valid HELM strings.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import parse_row
from peptide_design.generator import generate_helm
from peptide_design.validator import validate_helm

# PYY3-36 reference backbone (HELM positions 1..34)
_PYY_SEQ = list("IKPEAPGEDASPEE") + list("LNRYYASLRHYLNL") + list("VTRQRY")
assert len(_PYY_SEQ) == 34

def _base_row(name, backbone=None):
    row = {"Name": name}
    seq = backbone or _PYY_SEQ
    for i, aa in enumerate(seq, 1):
        row[str(i)] = aa
    return row

def _add_sc(row, site, bonds, monomers, suffix=""):
    row[f"Site{suffix}"] = str(site)
    for k, (b, m) in enumerate(zip(bonds, monomers), 1):
        row[f"b{k}{suffix}"] = b
        row[f"SC{k}{suffix}"] = m

def _run(row):
    c = parse_row(row)
    helm = generate_helm(c)
    ok, msg = validate_helm(helm)
    return helm, ok, msg

# ── Tests ───────────────────────────────────────────────────────────────────

def test_pyy_ref_no_sidechain():
    row = _base_row("ref")
    helm, ok, msg = _run(row)
    assert ok, msg
    assert helm.startswith("PEPTIDE1{I.K.P.E.A.P.G.E.D.A.S.P.E.E.L.N.R.Y.Y.A.S.L.R.H.Y.L.N.L.V.T.R.Q.R.Y}")

def test_pyy_analog_pos30_c18():
    """Analog_28: C18 gGlu-2xAdo at HELM pos 30 (K30)."""
    backbone = list(_PYY_SEQ); backbone[29] = "K"
    row = _base_row("Analog_28", backbone)
    _add_sc(row, 30, ["R3-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}" in helm
    assert "PEPTIDE2,PEPTIDE1,1:R1-30:R3" in helm

def test_pyy_analog_n_terminal_r1():
    """Analog_01-posNa: sidechain attaches to N-terminus via R1."""
    row = _base_row("Analog_01-posNa")
    _add_sc(row, 1, ["R1-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "PEPTIDE2,PEPTIDE1,1:R1-1:R1" in helm

def test_pyy_analog_c14d_no_gglu():
    """Analog_37: 2xAdo-C14d (no gGlu) at pos 30."""
    backbone = list(_PYY_SEQ); backbone[29] = "K"
    row = _base_row("Analog_37", backbone)
    _add_sc(row, 30, ["R3-R1","R2-R1","R2-R1"], ["Ado","Ado","C14d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "PEPTIDE2{[Ado].[Ado].[C14d]}" in helm

def test_pyy_combinatorial_backbone_plus_lipid():
    """Analog_52: Arg2 + Gln16 + K30 + C18 lipid (group C compound)."""
    backbone = list(_PYY_SEQ)
    backbone[1]  = "R"   # pos 2: Ile→Arg (paper: Arg4)
    backbone[15] = "Q"   # pos 16: Asn→Gln (paper: Gln18)
    backbone[29] = "K"   # pos 30: lipidation site
    row = _base_row("Analog_52", backbone)
    _add_sc(row, 30, ["R3-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "R." in helm   # Arg at pos 2
    assert ".Q." in helm  # Gln at pos 16
```

- [ ] **Step 2: Run integration tests**

```bash
python -m pytest tests/peptide_design/test_integration_pyy.py -v
```

Expected: all 5 tests PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/peptide_design/test_integration_pyy.py
git commit -m "test: PYY integration tests for peptide_design pipeline"
```

---

## Task 7: Excel template (make_template.py)

Creates a `.xlsx` with four sheets: **Design** (the authoring sheet), **HELM** (generated output), **MonomerDB** (for VBA lookup), and **Example** (pre-filled PYY analogs).

**Files:**
- Create: `peptide_design/make_template.py`

- [ ] **Step 1: Implement make_template.py**

```python
#!/usr/bin/env python3
# peptide_design/make_template.py
"""
Generate the Excel design template as a .xlsx file.
The MonomerDB sheet is pre-populated from the JSON files so the VBA
macro can look up polymer types without calling Python.

Usage:
    python -m peptide_design.make_template          # writes peptide_design_template.xlsx
    python -m peptide_design.make_template --out /path/to/file.xlsx
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

# ── Constants ────────────────────────────────────────────────────────────────

MAIN_CHAIN_POSITIONS = 34   # PYY3-36 length; user can extend
MAX_SC_MONOMERS = 8         # max sidechain depth
MAX_SIDECHAINS = 2          # sidechain blocks per compound

BLUE_FILL   = PatternFill("solid", fgColor="D9E1F2")
GREEN_FILL  = PatternFill("solid", fgColor="E2EFDA")
GREY_FILL   = PatternFill("solid", fgColor="F2F2F2")
ORANGE_FILL = PatternFill("solid", fgColor="FCE4D6")
BOLD        = Font(bold=True)


def _header(ws, col, row, value, fill=None, bold=True):
    cell = ws.cell(row=row, column=col, value=value)
    cell.font = Font(bold=bold)
    if fill:
        cell.fill = fill
    cell.alignment = Alignment(horizontal="center")
    return cell


def _build_design_headers(ws) -> list[str]:
    """Write row-1 headers and return the ordered list of header names."""
    col = 1
    headers = []

    # Name
    _header(ws, col, 1, "Name", GREY_FILL)
    ws.column_dimensions[get_column_letter(col)].width = 22
    headers.append("Name")
    col += 1

    # Main chain positions
    for pos in range(1, MAIN_CHAIN_POSITIONS + 1):
        _header(ws, col, 1, str(pos), BLUE_FILL)
        ws.column_dimensions[get_column_letter(col)].width = 6
        headers.append(str(pos))
        col += 1

    # Sidechain blocks
    for sc_idx in range(1, MAX_SIDECHAINS + 1):
        suffix = "" if sc_idx == 1 else f"_{sc_idx}"
        site_key = f"Site{suffix}"
        _header(ws, col, 1, site_key, ORANGE_FILL)
        ws.column_dimensions[get_column_letter(col)].width = 8
        headers.append(site_key)
        col += 1

        for k in range(1, MAX_SC_MONOMERS + 1):
            bk = f"b{k}{suffix}"
            sk = f"SC{k}{suffix}"
            _header(ws, col, 1, bk, GREEN_FILL)
            ws.column_dimensions[get_column_letter(col)].width = 8
            headers.append(bk)
            col += 1
            _header(ws, col, 1, sk, GREEN_FILL)
            ws.column_dimensions[get_column_letter(col)].width = 8
            headers.append(sk)
            col += 1

    ws.freeze_panes = "B2"
    return headers


def _write_example_rows(ws, headers):
    """Pre-fill a few PYY analogs so the format is self-documenting."""
    PYY_SEQ = list("IKPEAPGEDASPEE") + list("LNRYYASLRHYLNL") + list("VTRQRY")

    def row_dict(name, backbone=None, sc1=None, sc2=None):
        d = {"Name": name}
        seq = backbone or PYY_SEQ
        for i, aa in enumerate(seq, 1):
            d[str(i)] = aa
        if sc1:
            site, bonds, monomers = sc1
            d["Site"] = site
            for k, (b, m) in enumerate(zip(bonds, monomers), 1):
                d[f"b{k}"] = b
                d[f"SC{k}"] = m
        return d

    backbone_k30 = list(PYY_SEQ); backbone_k30[29] = "K"

    examples = [
        row_dict("PYY3-36-ref"),
        row_dict("Analog_28-pos30-C18", backbone_k30,
                 sc1=(30, ["R3-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])),
        row_dict("Analog_37-pos30-C14d-2xAdo", backbone_k30,
                 sc1=(30, ["R3-R1","R2-R1","R2-R1"], ["Ado","Ado","C14d"])),
        row_dict("Analog_01-posNa-N-term",
                 sc1=(1, ["R1-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])),
    ]

    for row_num, ex in enumerate(examples, start=2):
        for col_num, h in enumerate(headers, start=1):
            val = ex.get(h, "")
            if val:
                ws.cell(row=row_num, column=col_num, value=str(val))


def _write_monomer_db_sheet(wb):
    ws = wb.create_sheet("MonomerDB")
    _header(ws, 1, 1, "Symbol", GREY_FILL)
    _header(ws, 2, 1, "PolymerType", GREY_FILL)
    ws.column_dimensions["A"].width = 20
    ws.column_dimensions["B"].width = 14

    root = Path(__file__).parent.parent / "monomer_db"
    entries = []
    for fname in ("HELMCoreLibrary.json", "custom_monomers.json"):
        path = root / fname
        if path.exists():
            data = json.loads(path.read_text())
            if isinstance(data, list):
                entries.extend(data)

    seen = set()
    row = 2
    for e in entries:
        sym = e.get("symbol") or e.get("id", "")
        pt  = e.get("polymerType", "PEPTIDE")
        if sym and sym not in seen:
            ws.cell(row=row, column=1, value=sym)
            ws.cell(row=row, column=2, value=pt)
            seen.add(sym)
            row += 1


def _write_helm_sheet(wb):
    ws = wb.create_sheet("HELM")
    for col, (h, w) in enumerate([("Name", 22), ("HELM", 80), ("Status", 16)], 1):
        _header(ws, col, 1, h, GREY_FILL)
        ws.column_dimensions[get_column_letter(col)].width = w
    ws.freeze_panes = "A2"


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="peptide_design_template.xlsx")
    args = p.parse_args(argv)

    wb = openpyxl.Workbook()
    ws_design = wb.active
    ws_design.title = "Design"

    headers = _build_design_headers(ws_design)
    _write_example_rows(ws_design, headers)
    _write_monomer_db_sheet(wb)
    _write_helm_sheet(wb)

    # Example sheet (copy of Design pre-fill)
    ws_ex = wb.create_sheet("Example")
    for row in ws_design.iter_rows():
        for cell in row:
            ws_ex.cell(row=cell.row, column=cell.column, value=cell.value)
    ws_ex.freeze_panes = "B2"

    out = Path(args.out)
    wb.save(out)
    print(f"Template written → {out}")
    print("Next step: open in Excel and import helm_builder.bas via the VBE (Alt+F11 → File → Import)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run make_template.py and verify**

```bash
python -m peptide_design.make_template --out /tmp/peptide_design_template.xlsx
python3 -c "
import openpyxl
wb = openpyxl.load_workbook('/tmp/peptide_design_template.xlsx')
print('Sheets:', wb.sheetnames)
ws = wb['Design']
print('Design row 1 headers (first 8):', [ws.cell(1,c).value for c in range(1,9)])
print('Example row 2 name:', ws.cell(2,1).value)
wb2 = openpyxl.load_workbook('/tmp/peptide_design_template.xlsx')
mdb = wb2['MonomerDB']
print('MonomerDB row count:', mdb.max_row)
"
```

Expected:
```
Sheets: ['Design', 'MonomerDB', 'HELM', 'Example']
Design row 1 headers (first 8): ['Name', '1', '2', '3', '4', '5', '6', '7']
Example row 2 name: PYY3-36-ref
MonomerDB row count: 700+
```

- [ ] **Step 3: Commit**

```bash
git add peptide_design/make_template.py
git commit -m "feat: Excel template generator (make_template.py)"
```

---

## Task 8: VBA macro (helm_builder.bas)

Implements the same collapse/expand algorithm in VBA, reading the Design sheet and writing to the HELM sheet. The MonomerDB sheet (written by `make_template.py`) provides polymer-type lookups.

**Files:**
- Create: `helm_builder.bas`

- [ ] **Step 1: Write helm_builder.bas**

```vba
' helm_builder.bas
' Import into the Excel workbook via VBE (Alt+F11 → File → Import File)
' Then assign BuildHELM to a button on the Design sheet.
Option Explicit

' ── Public entry point ───────────────────────────────────────────────────────

Public Sub BuildHELM()
    Dim wsDesign As Worksheet, wsHelm As Worksheet
    Dim polyTypes As Object   ' Scripting.Dictionary: symbol → polymerType
    Set polyTypes = LoadMonomerDB()

    Set wsDesign = ThisWorkbook.Sheets("Design")
    Set wsHelm   = ThisWorkbook.Sheets("HELM")

    ' Clear HELM sheet below header
    If wsHelm.Cells(wsHelm.Rows.Count, 1).End(xlUp).Row > 1 Then
        wsHelm.Rows("2:" & wsHelm.Rows.Count).ClearContents
    End If

    ' Read Design sheet headers from row 1
    Dim headers() As String
    Dim lastCol As Long
    lastCol = wsDesign.Cells(1, wsDesign.Columns.Count).End(xlToLeft).Column
    ReDim headers(1 To lastCol)
    Dim c As Long
    For c = 1 To lastCol
        headers(c) = Trim(CStr(wsDesign.Cells(1, c).Value))
    Next c

    ' Process each data row
    Dim r As Long, outRow As Long
    outRow = 2
    r = 2
    Do While Trim(CStr(wsDesign.Cells(r, 1).Value)) <> ""
        Dim name As String
        name = CStr(wsDesign.Cells(r, 1).Value)

        ' Build column → value dict for this row
        Dim rowData As Object
        Set rowData = CreateObject("Scripting.Dictionary")
        For c = 1 To lastCol
            rowData(headers(c)) = Trim(CStr(wsDesign.Cells(r, c).Value))
        Next c

        Dim helm As String, status As String
        helm  = GenerateHELM(rowData, polyTypes)
        status = ValidateHELM(helm)

        wsHelm.Cells(outRow, 1).Value = name
        wsHelm.Cells(outRow, 2).Value = helm
        wsHelm.Cells(outRow, 3).Value = status
        outRow = outRow + 1
        r = r + 1
    Loop

    MsgBox "Done. " & (outRow - 2) & " compounds written to HELM sheet.", vbInformation
End Sub

' ── MonomerDB ────────────────────────────────────────────────────────────────

Private Function LoadMonomerDB() As Object
    Dim ws As Worksheet
    Dim db As Object
    Set db = CreateObject("Scripting.Dictionary")

    On Error Resume Next
    Set ws = ThisWorkbook.Sheets("MonomerDB")
    On Error GoTo 0
    If ws Is Nothing Then
        MsgBox "MonomerDB sheet not found. Run make_template.py to regenerate.", vbExclamation
        Set LoadMonomerDB = db
        Exit Function
    End If

    Dim lastR As Long
    lastR = ws.Cells(ws.Rows.Count, 1).End(xlUp).Row
    Dim i As Long
    For i = 2 To lastR
        Dim sym As String, pt As String
        sym = Trim(CStr(ws.Cells(i, 1).Value))
        pt  = Trim(CStr(ws.Cells(i, 2).Value))
        If sym <> "" Then db(sym) = pt
    Next i
    Set LoadMonomerDB = db
End Function

' ── Row → HELM ───────────────────────────────────────────────────────────────

Private Function GenerateHELM(rowData As Object, polyTypes As Object) As String
    ' 1. Build main chain
    Dim mainChain As String
    Dim pos As Long
    For pos = 1 To 200
        Dim posKey As String: posKey = CStr(pos)
        If Not rowData.Exists(posKey) Then Exit For
        Dim aa As String: aa = rowData(posKey)
        If aa = "" Then Exit For
        If Len(mainChain) > 0 Then mainChain = mainChain & "."
        mainChain = mainChain & FmtSym(aa)
    Next pos

    Dim chains As String, conns As String
    chains = "PEPTIDE1{" & mainChain & "}"
    conns  = ""

    Dim peptideN As Long, chemN As Long
    peptideN = 2: chemN = 1

    ' 2. Process up to 3 sidechain blocks
    Dim scIdx As Long
    For scIdx = 1 To 3
        Dim sfx As String
        sfx = IIf(scIdx = 1, "", "_" & scIdx)

        Dim siteKey As String: siteKey = "Site" & sfx
        If Not rowData.Exists(siteKey) Then Exit For
        Dim siteVal As String: siteVal = rowData(siteKey)
        If siteVal = "" Then Exit For
        Dim site As Long: site = CLng(siteVal)

        ' Collect bonds and monomers
        Dim bonds(1 To 20) As String, monomers(1 To 20) As String
        Dim nSC As Long: nSC = 0
        Dim k As Long
        For k = 1 To 20
            Dim bKey As String: bKey = "b" & k & sfx
            Dim mKey As String: mKey = "SC" & k & sfx
            If Not rowData.Exists(bKey) Then Exit For
            Dim bVal As String: bVal = rowData(bKey)
            Dim mVal As String: mVal = rowData(mKey)
            If bVal = "" Or mVal = "" Then Exit For
            nSC = nSC + 1
            bonds(nSC) = bVal
            monomers(nSC) = mVal
        Next k

        If nSC = 0 Then GoTo NextSC

        If CanCollapse(bonds, monomers, nSC, polyTypes) Then
            Dim scBody As String: scBody = ""
            For k = 1 To nSC
                If Len(scBody) > 0 Then scBody = scBody & "."
                scBody = scBody & FmtSym(monomers(k))
            Next k
            Dim chainId As String
            chainId = "PEPTIDE" & peptideN
            peptideN = peptideN + 1
            chains = chains & "|" & chainId & "{" & scBody & "}"
            Dim b1Parts() As String
            b1Parts = Split(bonds(1), "-")
            Dim proxRg As String: proxRg = b1Parts(0)
            If Len(conns) > 0 Then conns = conns & "|"
            conns = conns & chainId & ",PEPTIDE1,1:R1-" & site & ":" & proxRg
        Else
            ' Expand: one chain per monomer
            Dim prevId As String: prevId = "PEPTIDE1"
            Dim prevPos As String: prevPos = CStr(site)
            Dim isFirst As Boolean: isFirst = True

            For k = 1 To nSC
                Dim ptype As String
                ptype = "PEPTIDE"
                If polyTypes.Exists(monomers(k)) Then ptype = polyTypes(monomers(k))

                Dim curId As String
                If ptype = "CHEM" Then
                    curId = "CHEM" & chemN: chemN = chemN + 1
                Else
                    curId = "PEPTIDE" & peptideN: peptideN = peptideN + 1
                End If

                chains = chains & "|" & curId & "{" & FmtSym(monomers(k)) & "}"

                Dim bParts() As String
                bParts = Split(bonds(k), "-")
                Dim pRg As String: pRg = bParts(0)
                Dim dRg As String: dRg = bParts(1)

                Dim conn As String
                If isFirst Then
                    conn = curId & ",PEPTIDE1,1:" & dRg & "-" & site & ":" & pRg
                    isFirst = False
                Else
                    conn = curId & "," & prevId & ",1:" & dRg & "-1:" & pRg
                End If
                If Len(conns) > 0 Then conns = conns & "|"
                conns = conns & conn
                prevId = curId
            Next k
        End If

NextSC:
    Next scIdx

    If Len(conns) > 0 Then
        GenerateHELM = chains & "$" & conns & "$$$V2.0"
    Else
        GenerateHELM = chains & "$$$$V2.0"
    End If
End Function

' ── Helpers ──────────────────────────────────────────────────────────────────

Private Function FmtSym(sym As String) As String
    If Len(sym) = 1 Then
        FmtSym = sym
    Else
        FmtSym = "[" & sym & "]"
    End If
End Function

Private Function CanCollapse(bonds() As String, monomers() As String, _
                              nSC As Long, polyTypes As Object) As Boolean
    ' All PEPTIDE, b1 distal=R1, all bk>=2 == R2-R1
    Dim k As Long
    For k = 1 To nSC
        Dim pt As String: pt = "PEPTIDE"
        If polyTypes.Exists(monomers(k)) Then pt = polyTypes(monomers(k))
        If pt <> "PEPTIDE" Then
            CanCollapse = False: Exit Function
        End If
    Next k
    ' b1 distal Rg must be R1
    Dim b1p() As String: b1p = Split(bonds(1), "-")
    If b1p(1) <> "R1" Then
        CanCollapse = False: Exit Function
    End If
    ' All subsequent bonds R2-R1
    For k = 2 To nSC
        If bonds(k) <> "R2-R1" Then
            CanCollapse = False: Exit Function
        End If
    Next k
    CanCollapse = True
End Function

Private Function ValidateHELM(helm As String) As String
    ' Lightweight check: looks for required HELM V2.0 suffix and balanced braces
    If Right(helm, 6) <> "V2.0" Then
        ValidateHELM = "ERROR: missing V2.0 suffix"
        Exit Function
    End If
    Dim depth As Long: depth = 0
    Dim i As Long
    For i = 1 To Len(helm)
        Dim ch As String: ch = Mid(helm, i, 1)
        If ch = "{" Then depth = depth + 1
        If ch = "}" Then depth = depth - 1
        If depth < 0 Then
            ValidateHELM = "ERROR: unbalanced braces"
            Exit Function
        End If
    Next i
    If depth <> 0 Then
        ValidateHELM = "ERROR: unbalanced braces"
    Else
        ValidateHELM = "OK"
    End If
End Function
```

- [ ] **Step 2: Commit the VBA file**

```bash
git add helm_builder.bas
git commit -m "feat: VBA macro helm_builder.bas (Design sheet → HELM sheet)"
```

- [ ] **Step 3: Manual test in Excel**

```bash
python -m peptide_design.make_template --out /tmp/test_template.xlsx
```

Open `/tmp/test_template.xlsx` in Excel. In VBE (Alt+F11): File → Import File → select `helm_builder.bas`. Add a button to the Design sheet and assign it to `BuildHELM`. Click the button. Verify:
- HELM sheet row 1 (PYY3-36-ref): `PEPTIDE1{I.K.P...R.Y}$$$$V2.0`, Status: OK
- HELM sheet row 2 (Analog_28-pos30-C18): contains `PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}`, Status: OK

---

## Task 9: Final integration + push

- [ ] **Step 1: Run full test suite**

```bash
python -m pytest tests/peptide_design/ -v
```

Expected: all tests PASS (test_core, test_generator, test_validator, test_integration_pyy).

- [ ] **Step 2: Run CLI on PYY design CSV (round-trip)**

Convert `data/pyy_lipidation_scan.csv` to design format and verify HELM strings parse:

```bash
python3 -c "
import sys, csv
sys.path.insert(0,'.')
from scripts.helm_parser import HELMParser

errors = []
with open('data/pyy_lipidation_scan.csv') as f:
    for row in csv.DictReader(f):
        try:
            HELMParser.parse(row['HELM'])
        except Exception as e:
            errors.append((row['Name'], str(e)))

print(f'Checked {len(list(open(\"data/pyy_lipidation_scan.csv\")))-1} rows')
if errors:
    for name, err in errors: print(f'  ERR {name}: {err}')
else:
    print('All existing PYY HELMs parse OK (baseline confirmed)')
"
```

- [ ] **Step 3: Generate the final template**

```bash
python -m peptide_design.make_template --out peptide_design_template.xlsx
```

- [ ] **Step 4: Update .gitignore to exclude generated output but track template**

Add to `.gitignore`:
```
*_helm.csv
```

Do NOT gitignore `peptide_design_template.xlsx` — it should be tracked.

- [ ] **Step 5: Final commit and push**

```bash
git add peptide_design_template.xlsx .gitignore
git commit -m "feat: peptide design tool complete (Python CLI + Excel template + VBA)"
git push origin main
```

---

## Open items resolved by defaults

| Item | Default chosen |
|---|---|
| Second sidechain column headers | `Site_2`, `b1_2`, `SC1_2`, … |
| Python + VBA monomer DB source | Both read `monomer_db/` JSON files; VBA caches into MonomerDB sheet |
| `peptide_design/` location | Subdirectory of `helm_sar/` |
| Max sidechain depth in template | 8 monomers (extendable by user) |
