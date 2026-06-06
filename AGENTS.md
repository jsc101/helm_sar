# HELM SAR Agent — Mental Model & Tool Guide

This file tells an AI agent (Claude Code or similar) everything needed to use this repo
to generate sequence-alignment SAR reports for HELM-encoded peptide libraries.

---

## Full Pipeline

```
Paper (PDF / HTML)
      │
      ▼
[Step 1 — LLM agent, interactive]
  Read paper tables, figures, and footnotes.
  For each compound construct a HELM string and extract all assay values.
  Output: library.csv  (Name, HELM, assay_col_1, assay_col_2, …)
      │
      ▼
[Step 2 — run_sar.py]
  Parse HELM → HELMObject → alignment → colours → plots → HTML
      │
      ▼
  report.html  (interactive, self-contained)
```

**Step 1 is intentionally agent-driven** — academic papers encode SAR data in too many
formats (R-group substituent tables, multi-panel figures, supplementary spreadsheets,
footnote-defined modifications) for a rigid script to handle reliably.  An LLM agent reads
all of it, applies chemical reasoning to construct valid HELM strings, and writes a clean
CSV that Step 2 can consume directly.

**Step 2 is deterministic** — one command:
```
python run_sar.py --input library.csv --activity IC50_nM
```

---

## What This Repo Does (Step 2 detail)

A library.csv with HELM-encoded compounds and optional activity columns.  This repo:

1. Parses every HELM string → `HELMObject` (full structural graph)
2. Identifies the reference compound automatically
3. Rotates cyclic peptides to their best-aligned orientation (cosine descriptor similarity)
4. NW-aligns main-chain sequences to the reference; sidechain columns inserted structurally
5. Calculates MW and per-residue physicochemical descriptors via RDKit
6. Renders an interactive HTML report: two colour schemes, sortable columns, conservation bars,
   scatter plots (lipidation position vs activity; sidechain string vs activity)

---

## HELM V2 Notation — What the Agent Must Know

HELM is a compact text format for peptide structures.

```
PEPTIDE1{A.K.[meA].G.[D-Phe]}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$V2.0
          └───chain───────┘    └──────────connections──────────────┘
```

**Monomer encoding:**
| Notation        | Meaning                             | Example         |
|-----------------|-------------------------------------|-----------------|
| `A`             | Standard amino acid (single letter) | Alanine         |
| `[meA]`         | Modified residue (bracketed)        | N-methyl-Ala    |
| `[D-Phe]`       | D-amino acid                        | D-Phenylalanine |
| `[Pip]`         | Non-proteinogenic (symbol in DB)    | Pipecolinic acid|

**Connection syntax:** `CHAIN,CHAIN,from_pos:from_Rgroup-to_pos:to_Rgroup`
- `R1` = N-terminal attachment (alpha-amino group)
- `R2` = C-terminal attachment (alpha-carboxyl group)
- `R3` = sidechain attachment (e.g. Lys ε-amine, Asp sidechain COOH)

**Topology examples:**
```
# Monocyclic (backbone ring)
PEPTIDE1{A.K.G.F}$PEPTIDE1,PEPTIDE1,1:R1-4:R2$$$V2.0

# Bicycle (backbone ring + isopeptide sidechain)
PEPTIDE1{A.F.G.K.[meI].G.[meF].T.[meA].L.D.[Pip]}$PEPTIDE1,PEPTIDE1,1:R1-12:R2|PEPTIDE1,PEPTIDE1,4:R3-11:R3$$$V2.0

# Multi-chain with sidechain peptide (e.g. GLP-1 analog with fatty acid linker)
PEPTIDE1{H.[Aib].E.G.T.F.T.S.D.V.S.S.Y.L.E.G.Q.A.A.K.E.F.I.A.W.L.V.R.G.R.G}|PEPTIDE2{[gGlu].G.G.[Ph]}$PEPTIDE2,PEPTIDE1,1:R1-20:R3$$$V2.0
```

---

## JPV — Human-Readable Structure Notation

JPV (JSON Peptide Visualization) is the human-readable linearisation of a HELM structure.
Call `obj.get_jpv()` to get it.

```python
obj = HELMParser.parse(helm)
print(obj.get_jpv())
# → H-Aib-E-G-T-F-T-S-D-V-S-S-Y-L-E-G-Q-A-A-K(1,3)(gGlu(1,1)-G-G-Ph)-E-F-I-A-W-L-V-R-G-R-G
```

**Reading JPV:**
- Residues separated by `-`
- `K(1,3)` = Lys connected at R1 (backbone) and R3 (sidechain)
- `(gGlu(1,1)-G-G-Ph)` = sidechain branch: gGlu connected at R1 only, reading to Ph
- Cyclic rings shown as `(1,2)` at both ends — position 1 connects to last

JPV is the best way to visually verify a structure before running analysis.

---

## HELMObject — Core Data Structure

```python
from scripts.helm_parser import HELMParser
obj = HELMParser.parse(helm_string)   # Never pass db= kwarg; not supported
```

**Key methods:**

| Method | Returns | Use |
|--------|---------|-----|
| `obj.get_jpv()` | `str` | Human-readable sequence with topology annotations |
| `obj.get_jpv_flat()` | `list[dict]` | **Flat token list for alignment** — main chain + sidechains inserted in outward order; each token has `symbol`, `label`, `is_main`, `main_pos`, `sc_pos` |
| `obj.get_lipidation_pos()` | `int \| None` | **1-based main chain position bearing a secondary chain** (e.g. 28 for a Lys30-lipidated PYY analogue) |
| `obj.get_sidechain_string()` | `str` | **Dash-separated sidechain monomers reading outward** from backbone (e.g. `"gGlu-Ado-Ado-C18d"`); `""` for single-chain |
| `obj.get_chain(chain_id=None)` | `dict` | Primary (longest) PEPTIDE chain; monomers list |
| `obj.positions(chain_id=None)` | `list[int]` | Position indices in the chain |
| `obj.all_monomer_descriptors(chain_id=None)` | `dict[int, dict]` | Per-position RDKit descriptors |
| `obj.data['_chains']` | `list[dict]` | All chains; each has `chain_id`, `monomers` |
| `obj.data['connectivity_graph']` | `list[dict]` | All bonds (from_chain, from_pos, from_rgroup, ...) |

**Monomer dict structure:**
```python
{'pos': 3, 'symbol': 'meA', 'chain_id': 'PEPTIDE1'}
```

**Descriptor dict structure** (from `all_monomer_descriptors`):
```python
{
  'in_db': True,
  'descriptors': {'MW': 89.1, 'LogP': -0.49, 'TPSA': 63.3,
                  'HBD': 2, 'HBA': 3, 'RotBonds': 0,
                  'AromaticRings': 0, 'QED': 0.56}
}
```

---

## Sequence Flattening (Multi-Chain)

`obj.get_jpv_flat()` returns a **flat token list** with sidechain monomers inserted after
their backbone attachment position, reading **outward from the backbone**.  
`_per_pos_syms(obj)` in `scripts/report.py` is a thin wrapper that returns just the symbols.

For the GLP-1 analog above (PEPTIDE2 attaches at K pos 20):
```
Primary:   H Aib E G T F T S D V S S Y L E G Q A A K                   E F I A W L V R G R G
Sidechain at pos 20 (outward):                          gGlu G G Ph
Flattened: H Aib E G T F T S D V S S Y L E G Q A A K gGlu G G Ph E F I A W L V R G R G
Labels:    1   2 3 4 5 6 7 8 9 ...                  20 20.1 20.2 20.3 20.4 21 ...
```

**Token dict keys:**

| Key | Type | Example |
|-----|------|---------|
| `symbol` | `str` | `"gGlu"`, `"K"` |
| `label` | `str` | `"20"` (main chain), `"20.1"` (1st sidechain monomer at pos 20) |
| `is_main` | `bool` | `True` for backbone, `False` for sidechain |
| `main_pos` | `int` | 1-based main chain position (same for backbone and its sidechain tokens) |
| `sc_pos` | `int` | 0 for backbone tokens; 1, 2, … for sidechain tokens in outward order |

**Sidechain ordering:** position 1 in PEPTIDE2 must be the monomer closest to the backbone.
If the attachment is at PEPTIDE2's last position (e.g. reversed chain), `get_jpv_flat()` and
`get_sidechain_string()` automatically reverse the order so labels always read outward.

A compound without the sidechain (31 residues) aligns against this 35-position layout,
and NW places gaps at the 4 sidechain positions.

---

## Alignment Pipeline

### Step 1 — Rotation (cyclic peptides only)

For two cyclic peptides of the same length, `HELMAlignment` finds the rotation k
that maximises per-position cosine similarity of 8-descriptor vectors:

```python
from scripts.helm_alignment import HELMAlignment
aligner = HELMAlignment(ref_obj)
result  = aligner.align(query_obj)
# result.best_k      → int, rotation offset applied
# result.best_score  → float 0–1, cosine similarity
# result.rotated     → HELMObject with residues in aligned order
# result.rotated_jpv → str, JPV of rotated compound
```

Descriptor set and scale factors (controls angular weight in cosine space):
```
MW(500) LogP(8) TPSA(200) HBD(5) HBA(10) RotBonds(10) AromaticRings(3) Chiral(1)
```

Chirality flag: 1.0 for D-amino acids (symbol starts with `d` or `D-`), 0.0 for L/achiral.

### Step 2 — NW Sequence Alignment

After rotation, `_nw_align(ref_syms, query_syms)` performs global Needleman-Wunsch
alignment (match=2.0, mismatch=0.0, gap=-1.5).

`_align_group_rows(grp_ref_syms, rows, grp_ref_labels=None)` runs NW for all rows and
merges gap patterns into a shared master column layout.  Pass the reference position labels
(from `_per_pos_labels(ref_obj)`) to get labelled column headers in the report:

```python
n_cols, master_ref, master_labels, al_syms, al_desc = _align_group_rows(
    ref_syms, rows, ref_labels)
# master_ref[c]    = None → gap column (shown as grey · in header)
# master_labels[c] = None → gap column; "28" → main pos; "28.1" → sidechain sub-pos
# al_syms[i][c]    = None → this compound has a gap at column c
```

Column headers are styled by label type: plain integers for backbone positions, green
background for sidechain sub-positions (`N.1`, `N.2`, …).

### Display cutoff

All rows are aligned to the reference in one table. Rows are excluded if:
- Identity to reference < `min_identity` (default 0.3 = 30%) — too dissimilar to be informative
- Display count already reached `max_rows` (default 50) — table would be too wide to read

Remaining rows are sorted by identity descending. Excluded compounds are listed below the table.

Pass `--max-rows` and `--min-identity` to `run_sar.py` to override the defaults.

---

## MW Calculation

```python
from scripts.report import calc_mw
mw = calc_mw(obj)  # float | None
```

Formula: `MW = Σ(residue_MW) − n_bonds × 18.015`

Where `n_bonds = (n_residues − 1) + n_explicit_connections`.
- A 14-mer bicycle with 2 explicit connections: (13 + 2) = 15 bonds.
- Returns `None` if any residue is missing from the MonomerDB.

---

## MonomerDB

The toolkit ships with `monomer_db/HELMCoreLibrary.json` (Pistoia Alliance standard set,
~650 common proteinogenic and non-proteinogenic residues).

**If your library contains residues not in the standard set, you must supply an additional
JSON file.** Ask the collaborator: *"Do your HELM strings contain any non-standard residues
in brackets, e.g. `[meA]`, `[D-Phe]`, `[gGlu]`? If so, please provide a monomer JSON file
for those residues in Pistoia Alliance format."*

### Pistoia Alliance monomer JSON format

Each entry is a JSON object in an array:

```json
[
  {
    "symbol":      "meA",
    "name":        "N-Methyl Alanine",
    "polymerType": "PEPTIDE",
    "monomerType": "Backbone",
    "smiles":      "C[C@@H](N[H:1])C([OH:2])=O",
    "rgroups": [
      {"label": "R1", "capGroupSmiles": "[*:1][H]"},
      {"label": "R2", "capGroupSmiles": "O[*:2]"}
    ]
  }
]
```

| Field | Required | Notes |
|-------|----------|-------|
| `symbol` | yes | Must match HELM string exactly (e.g. `meA` for `[meA]`) |
| `polymerType` | yes | `"PEPTIDE"` for amino acids |
| `smiles` | yes | Attachment points as `[H:1]`/`[OH:1]` or `[*:1]` dummy atoms |
| `rgroups` | yes | R1 = N-term cap, R2 = C-term cap, R3 = sidechain |
| `name` | no | Human-readable name |
| `monomerType` | no | `"Backbone"` or `"Sidechain"` |

Loading extra monomers:
```python
db = MonomerDB(extra_sources=['path/to/my_monomers.json'])
entry = db.find_by_symbol('meA')   # → dict or None
entry = db.find('C[C@H](N[*:1])C(=O)[*:2]')  # SMILES lookup
```

**Symbol lookup fallback**: `D_Phe` → tries `D-Phe` (underscore→dash), handles HELMParser
normalisation transparently.

## Reference Compound

**Always ask the collaborator which compound is the reference / parent / wildtype.**

The agent auto-detects by scanning the Name column for: `parent`, `ref`, `wt`, `wildtype`,
`reference` (case-insensitive). If none match, the first row is used.

If the collaborator's names don't follow this convention, ask:
*"Which compound should I use as the alignment reference? Please give me the exact Name
as it appears in your table."* Then pass it with `--ref "ExactName"`.

---

## Colour Schemes

### Zappo (chemistry class)
| Class | Colour | Examples |
|-------|--------|---------|
| Aliphatic | Pink/salmon | A, V, L, I, M |
| Aromatic | Orange | F, W, Y |
| Positive | Blue | K, R, H |
| Negative | Red | D, E |
| Hydrophilic | Green | S, T, N, Q |
| Conformational | Magenta | G, P |
| Cysteine | Yellow | C |

### Charge / LogP
| Class | Colour | Notes |
|-------|--------|-------|
| Positive at pH 7 | Red | K, R, H, Dab, Orn (structural heuristic) |
| Negative at pH 7 | Blue | D, E, gGlu |
| Polar / Neutral | Near-white | S, T, N, Q, G, P, C |
| Aliphatic/Aromatic | Yellow→Orange | Scaled by RDKit LogP (t = LogP/4) |

*Note: charge classification is structural (SMARTS-based), not electrochemical. For true
pH-7 charge state you would need Epik or Henderson-Hasselbalch protonation.*

---

## Fatty Acid Lipidation — Multi-Chain HELM Convention

When encoding peptides with a fatty acid protractor (e.g. GLP-1, PYY, semaglutide analogs),
use a two-chain HELM where **PEPTIDE2 position 1 is always the monomer closest to the backbone**.

```
PEPTIDE1{backbone}|PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}$PEPTIDE2,PEPTIDE1,1:R1-28:R3$$V2.0
                                                           ^^^^                   ^^^^
                                          PEPTIDE2 pos 1 attaches            to Lys R3
```

Standard sidechain reading order (backbone → distal end): `gGlu → Ado → Ado → C18d`

| Connection type | Example connection string | Notes |
|-----------------|--------------------------|-------|
| Lys ε-amine (R3) | `PEPTIDE2,PEPTIDE1,1:R1-28:R3` | Standard; HELM pos 28 = PYY Lys30 |
| N-terminal (R1) | `PEPTIDE2,PEPTIDE1,1:R1-1:R1` | Nα acylation; R1→R1 convention |
| Reversed chain | `PEPTIDE2,PEPTIDE1,4:R2-1:R1` | gGlu is PEPTIDE2 pos 4; both `get_sidechain_string()` and `get_jpv_flat()` auto-reverse to report `gGlu-…` |

**Custom monomers in `monomer_db/custom_monomers.json`** (shipped with this repo):

| Symbol | Description |
|--------|-------------|
| `gGlu` | γ-glutamic acid linker (alpha-N attaches to Lys R3 or chain) |
| `Ado` | 8-amino-3,6-dioxaoctanoic acid (PEG spacer) |
| `C14d` – `C20d` | C14/C16/C18/C20 diacid terminal monomers |
| `meR` | N-methyl-arginine |
| `AcI` | N-acetyl-isoleucine (N-terminal acetylation cap) |

---

## Scatter Plots in the Report

`build_html()` calls `_make_scatter_plots(rows, activity_map, act_label)` which generates
two embedded PNG plots (base64, rendered inline):

1. **Activity vs. lipidation position** — x = `obj.get_lipidation_pos()` + 2 (PYY numbering),
   y = activity value. One dot per analogue with a sidechain.

2. **Activity vs. sidechain composition** — x = `obj.get_sidechain_string()` (categorical,
   rotated 90° labels), y = activity value.

Both plots appear above the alignment table. Compounds missing activity data are omitted.

---

## Step 1 — Extracting a Paper into library.csv (Agent Instructions)

This is an interactive task.  The user provides a paper (PDF, URL, or pasted text).
Your job is to produce a valid `library.csv` that `run_sar.py` can consume.

### Process

1. **Read the paper** — scan all tables, figures, and footnotes.  Note:
   - Substituent/R-group tables (e.g. "Table 1: R1 = H, Phe, Trp…")
   - Activity columns and units (IC50 nM, pIC50, %inh, EC50, half-life, etc.)
   - Which compound is the parent/reference

2. **Identify the backbone** — write out the parent sequence in one-letter or three-letter
   HELM symbols.  Map unusual residues to custom monomer symbols (see `monomer_db/custom_monomers.json`).

3. **Build HELM strings** — for each compound:
   - Start from the parent backbone
   - Apply substitutions from the SAR table
   - For lipidated analogues: add a PEPTIDE2 chain (see "Fatty Acid Lipidation" section below)
   - Verify the HELM string parses: `HELMParser.parse(helm_string)` should not raise

4. **Assemble the CSV** — columns: `Name`, `HELM`, then one column per assay.
   Use `''` (empty string) for missing values, never `None` or `NaN`.

5. **Show the user a preview** — print the first 5 rows and ask them to confirm before
   writing the file.

6. **Run Step 2** — once confirmed: `python run_sar.py --input library.csv --activity <col>`

### Common paper formats and how to handle them

| Format | Approach |
|--------|----------|
| R-group table (R1, R2 at defined positions) | Enumerate all combinations; substitute into backbone template |
| Scaffold + list of analogues | Use scaffold as parent; each analogue applies named changes |
| Lipidation scan (position × fatty acid) | Multi-chain HELM; PEPTIDE2 = protractor sequence |
| IC50 reported as range or inequality | Store as empty; note in a `Notes` column |
| Activity in figure only | Estimate from bar chart / curve; flag with `~` prefix or `Notes` |
| Stereochemistry specified | Use D-amino acid symbols (e.g. `dA`, `dF`) from monomer DB |

### Validation checklist before writing CSV

- [ ] Every HELM string parses without error
- [ ] Reference compound row present (`Name` contains "parent", "ref", or "wt")
- [ ] All activity values are numeric or empty string
- [ ] Monomer symbols that are not standard amino acids exist in `monomer_db/custom_monomers.json`
  (add them if needed — see "Custom Monomers" section in this file)
- [ ] Column headers match what you will pass to `--activity`

---

## Running the Agent — Common Tasks

### Generate a report from CSV
```
python run_sar.py --input library.csv --activity IC50_nM --out report.html
```

### Generate a report from Excel
```
python run_sar.py --input library.xlsx --activity pIC50
```

### Override the reference compound
```
python run_sar.py --input library.csv --ref "Compound_1"
```

### Use as a library in Python
```python
import sys; sys.path.insert(0, '.')
from monomer_db.monomer_db import MonomerDB
from scripts.helm_parser import HELMParser
from scripts.report import build_data, build_html

db = MonomerDB(extra_sources=['monomer_db/cycpeptmpdb_monomers.json'])
pairs = [('Parent', 'PEPTIDE1{A.K.G.F}$PEPTIDE1,PEPTIDE1,1:R1-4:R2$$$V2.0'),
         ('Analog_1', 'PEPTIDE1{A.R.G.F}$PEPTIDE1,PEPTIDE1,1:R1-4:R2$$$V2.0')]
ref_obj, rows = build_data(pairs, db)
html = build_html(ref_obj, rows, db=db,
                  activity_map={'Analog_1': 12.3}, act_label='IC50 (nM)')
open('report.html', 'w').write(html)
```

### Inspect a HELM string
```python
from scripts.helm_parser import HELMParser
obj = HELMParser.parse("PEPTIDE1{A.K.[meA]}$PEPTIDE1,PEPTIDE1,1:R1-3:R2$$$V2.0")
print(obj.get_jpv())
print(obj.all_monomer_descriptors())

# Multi-chain: inspect lipidation
helm = ("PEPTIDE1{I.K.P.E.A.P.G.E.D.A.S.P.E.E.L.N.R.Y.Y.A.S.L.R.H.Y.L.N.K.V.T.R.Q.R.Y}"
        "|PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}$PEPTIDE2,PEPTIDE1,1:R1-28:R3$$V2.0")
obj = HELMParser.parse(helm)
print(obj.get_lipidation_pos())    # → 28
print(obj.get_sidechain_string())  # → "gGlu-Ado-Ado-C18d"
flat = obj.get_jpv_flat()
# flat[27] = {'symbol':'K',    'label':'28',   'is_main':True,  'main_pos':28, 'sc_pos':0}
# flat[28] = {'symbol':'gGlu', 'label':'28.1', 'is_main':False, 'main_pos':28, 'sc_pos':1}
```

---

## SAR Interpretation Guide

When reading the HTML report, look for:

1. **% Identity column** — compounds with ≥93% identity (green bar) are close to parent.
   Sudden drops usually mean multiple simultaneous changes.

2. **Align score** — cosine similarity of descriptor vectors after rotation alignment.
   Score <0.85 typically means a chemically significant substitution.

3. **Changes vs Ref column** — explicit list of substitutions (P6:F→G means position 6
   changed from Phe to Gly).

4. **Conservation bars** — height = frequency of the most common residue at that position.
   Low conservation = position tolerates variation; high = position is critical.

5. **Zappo colour pattern** — scan rows for colour changes vs the reference (top row).
   A charge class change (blue↔red) at a position often signals a big SAR shift.

6. **Activity column** — colour-coded green (potent) to red (inactive).
   Sort by activity then compare sequence patterns in potent vs inactive rows.

7. **Excluded compounds** — listed below the table with their identity score.
   Very low identity (<30%) usually means a scaffold change; consider running those
   separately with their own reference.

---

## Input Data Checklist (for the collaborator)

- [ ] One row per compound
- [ ] Column named `Name` (or `ID`, `Compound`) — unique compound identifier
- [ ] Column named `HELM` — valid HELM V2 string
- [ ] (Optional) one or more activity columns (IC50, pIC50, etc.)
- [ ] Reference/parent compound identifiable by name ("Parent", "Ref", "WT") or listed first
- [ ] File format: CSV, Excel (.xlsx), or Apple Numbers (.numbers)

---

## File Structure

```
helm_sar/
├── run_sar.py              ← Entry point; run this
├── AGENTS.md               ← This file
├── README.md               ← Human quickstart
├── requirements.txt        ← Python dependencies
├── scripts/
│   ├── helm_engine.py      ← HELMObject data structure
│   ├── helm_parser.py      ← HELM V2 → HELMObject
│   ├── helm_alignment.py   ← Cosine rotation alignment + NW
│   ├── helm_compare.py     ← Structural identity check
│   ├── residue_colors.py   ← Zappo + charge/LogP colour schemes
│   ├── rdkit_bridge.py     ← RDKit descriptor calculation
│   └── report.py           ← HTML report generator (build_data, build_html)
└── monomer_db/
    ├── monomer_db.py       ← MonomerDB class
    ├── HELMCoreLibrary.json    ← Standard HELM monomers
    ├── custom_monomers.json    ← Custom additions
    └── cycpeptmpdb_monomers.json ← CycPeptMPDB non-proteinogenic residues
```
