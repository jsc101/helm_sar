# HELM SAR Toolkit

Two complementary tools for HELM-encoded peptide work:

- **SAR Report pipeline** — extract a paper into a HELM library, then generate an interactive HTML alignment + activity report
- **Peptide Design tool** — author new compounds in a flat Excel grid and generate valid HELM strings with correct sidechain connectivity

---

## Tool 1: SAR Report Pipeline

### Overview

```
Paper (PDF / URL)
      │
      ▼
[Step 1 — interactive LLM agent]
  Reads tables, figures, footnotes.
  Constructs HELM strings and extracts assay values.
  Output: library.csv  (Name, HELM, assay columns…)
      │
      ▼
[Step 2 — run_sar.py]
  HELM parsing → alignment → colours → scatter plots → HTML
      │
      ▼
  report.html  (interactive, self-contained)
```

Step 1 is agent-driven — papers encode SAR in too many formats for a rigid script.
An LLM agent reads all formats, builds valid HELM strings, and writes the CSV.
See `AGENTS.md` for full agent instructions.

### Quickstart

```bash
pip install -r requirements.txt

# Generate a report
python run_sar.py --input library.csv --activity IC50_nM

# From Excel or Apple Numbers
python run_sar.py --input library.xlsx --activity pIC50 --out report.html
python run_sar.py --input library.numbers
```

### Input format

| Column | Required | Notes |
|--------|----------|-------|
| `Name` | yes | Also accepts: `ID`, `Compound`, `Cmpd` |
| `HELM` | yes | HELM V2 string. Also accepts: `HELM_string`, `Sequence` |
| activity | no | Any name — pass with `--activity` |

Reference compound is auto-detected from the `Name` column: any name containing
`parent`, `ref`, `wt`, `wildtype`, or `reference` (case-insensitive).
If none matches, the first row is used.

### Report features

- **Sequence alignment** — coloured by Zappo chemistry class (aromatic, aliphatic, charged, polar…)
- **Charge/LogP map** — same alignment, coloured red (positive) / blue (negative) / yellow-orange (lipophilic)
- **Sidechain columns** — main chain positions as integers; sidechain sub-positions as `N.1`, `N.2`, … (green headers, reading outward from backbone)
- **Scatter plots** — activity vs lipidation position; activity vs sidechain composition
- **Activity column** — colour-coded potent (green) → inactive (red), sortable
- **MW column** — Σ(residue MW) − n_bonds × 18.015
- **% Identity** — fraction of positions matching the reference
- **Changes vs Ref** — explicit substitution list per compound
- **Conservation bars** — per-position residue frequency
- **Gap columns** — shorter/longer compounds get NW gap cells, not truncation
- **Sortable columns** — click any header; reference and consensus rows stay pinned

---

## Tool 2: Peptide Design (Excel → HELM)

### Overview

Design new analogs in a flat Excel grid — one row per compound, one column per
backbone position — then generate valid HELM V2 strings automatically.

The key challenge is sidechain connectivity: a plain monomer list loses which
R-groups form each bond, which matters as soon as a CHEM monomer (e.g. a click
triazole) or a non-standard bond appears. The design format captures this
explicitly with alternating `(bond, monomer)` pairs.

### Design format — one row per compound

```
| Name | 1 | 2 | … | N | Site | b1    | SC1  | b2    | SC2 | b3    | SC3 | b4    | SC4  |
|------|---|---|---|---|------|-------|------|-------|-----|-------|-----|-------|------|
| ref  | I | K | … | Y |      |       |      |       |     |       |     |       |      |
| A28  | I | K | … | K | 30   | R3-R1 | gGlu | R2-R1 | Ado | R2-R1 | Ado | R2-R1 | C18d |
```

- **Main chain columns** (`1`…`N`) — one residue symbol per HELM position
- **Site** — HELM position where the sidechain attaches (e.g. `30` for Lys30)
- **b`k`** — bond entering monomer `k`, written `proxRg-distRg` (e.g. `R3-R1`)
- **SC`k`** — monomer symbol, reading outward from backbone
- Up to 3 sidechains per compound: repeat `Site/b/SC` block with suffix `_2`, `_3`

Click-chemistry example (non-default bonds, mixed PEPTIDE+CHEM):

```
| Site | b1    | SC1 | b2    | SC2        | b3    | SC3 | b4    | SC4 | b5    | SC5  |
| 37   | R3-R1 | Hpg | R3-R1 | Triazole14 | R4-R3 | Aha | R2-R1 | A   | R2-R1 | C18d |
```

### CLI

```bash
# CSV or XLSX → adds HELM and HELM_status columns
python -m peptide_design.cli my_design.csv
python -m peptide_design.cli my_design.xlsx

# Print to terminal
python -m peptide_design.cli my_design.csv --stdout
```

### Excel deliverable

A ready-to-use `.xlsm` template with example PYY analogs pre-filled:

```bash
python -m peptide_design.make_template            # → peptide_design_template.xlsx
python -m peptide_design.make_template --out /path/to/file.xlsx
```

Open in Excel → VBE (Alt+F11) → File → Import → `helm_builder.bas` → assign
`BuildHELM` to a button on the Design sheet. The macro reads the Design sheet
and writes HELM strings + status to the HELM sheet, no Python required.

### HELM generation rule

| Branch type | Condition | Output |
|-------------|-----------|--------|
| **Collapse** | All monomers PEPTIDE + all bonds default head-to-tail (`R2-R1`, except b1 distal = `R1`) | Single `PEPTIDE{N}{[SC1].[SC2]…}` chain + one connection line |
| **Expand** | Any CHEM monomer, or any non-default bond | One chain per monomer + explicit bond lines |

The user writes the same format either way; the tool decides.

---

## Project Structure

```
helm_sar/
│
├── run_sar.py                    ← SAR report entry point
├── AGENTS.md                     ← Full agent instructions (Step 1 + Step 2 detail)
│
├── peptide_design/               ← Design → HELM module
│   ├── core.py                   parse_row() → Compound/Sidechain dataclasses
│   ├── generator.py              generate_helm() — collapse + expand paths
│   ├── validator.py              validate_helm() — HELMParser + monomer DB check
│   ├── cli.py                    CSV/XLSX in → HELM + status column out
│   └── make_template.py          generates peptide_design_template.xlsx
│
├── peptide_design_template.xlsx  ← Ready-to-use Excel (Design/HELM/MonomerDB/Example sheets)
├── helm_builder.bas              ← VBA macro (same algorithm, standalone in Excel)
│
├── scripts/                      ← SAR report engine
│   ├── helm_parser.py            HELM V2 string → HELMObject
│   ├── helm_engine.py            HELMObject — query methods (get_jpv_flat, get_lipidation_pos…)
│   ├── helm_alignment.py         Cosine rotation alignment + Needleman-Wunsch
│   ├── helm_compare.py           Structural identity
│   ├── residue_colors.py         Zappo + charge/LogP colour schemes
│   ├── rdkit_bridge.py           RDKit descriptor calculation
│   └── report.py                 HTML report generator (build_data, build_html)
│
├── monomer_db/                   ← Shared by both tools
│   ├── monomer_db.py             MonomerDB class — symbol/SMILES lookup
│   ├── HELMCoreLibrary.json      705 standard monomers (Pistoia Alliance)
│   ├── custom_monomers.json      Aib, aMePhe, Triazole14, gGlu, Ado, C12d–C20d…
│   └── cycpeptmpdb_monomers.json CycPeptMPDB non-proteinogenic residues
│
├── data/                         ← Step-1 extraction outputs
│   ├── pyy_lipidation_scan.csv   PYY3-36 lipidation scan (52 analogs)
│   ├── medi7219_sar.csv          GLP-1 / MEDI7219 series (5 compounds)
│   └── *.html                    Generated SAR reports
│
└── tests/peptide_design/         ← 30 tests (core, generator, validator, PYY integration)
```

## Custom monomers

Both tools share `monomer_db/custom_monomers.json`:

| Symbol | Description |
|--------|-------------|
| `gGlu` | γ-glutamic acid linker |
| `Ado` | 8-amino-3,6-dioxaoctanoic acid (PEG spacer) |
| `C12d`–`C20d` | C12/C14/C16/C18/C20 diacid terminal monomers |
| `Aib` | 2-Aminoisobutyric acid (α-methylalanine) |
| `aMePhe` | α-Methyl-L-phenylalanine |
| `aMeSer` | α-Methyl-L-serine |
| `aMeLys` | α-Methyl-L-lysine (ε-N for lipidation) |
| `Triazole14` | 1,4-disubstituted 1,2,3-triazole (CuAAC product) — CHEM type |
| `Hpg` / `Pra` | Homopropargylglycine / propargylglycine (alkyne handles) |
| `Aha` | L-Azidohomoalanine (azide handle) |
| `meR` | N-methyl-arginine |
| `AcI` | N-acetyl-isoleucine (N-terminal cap) |

## Dependencies

```
rdkit
numpy
openpyxl            # Excel read/write
numbers-parser      # Apple Numbers support
```
