# HELM SAR Toolkit

Generate interactive sequence-alignment + SAR reports for HELM-encoded peptide libraries.

## Pipeline Overview

```
Paper (PDF / URL)
      │
      ▼
[Step 1 — interactive LLM agent]
  Reads tables, figures, footnotes.
  Constructs HELM strings and extracts assay values.
  Output: library.csv
      │
      ▼
[Step 2 — run_sar.py]
  HELM parsing → alignment → colours → scatter plots → HTML
      │
      ▼
  report.html  (interactive, self-contained)
```

Step 1 is intentionally agent-driven — academic papers encode SAR data in too many formats
for a rigid script.  An LLM agent reads all formats, applies chemical reasoning to build
valid HELM strings, and writes the CSV.  See `AGENTS.md` for the agent's full instructions.

## Quickstart (Step 2 — report generation)

```bash
# Install dependencies
pip install -r requirements.txt

# Generate a report from a CSV file
python run_sar.py --input library.csv

# Include activity data
python run_sar.py --input library.csv --activity IC50_nM

# From Excel or Apple Numbers
python run_sar.py --input library.xlsx --activity pIC50 --out report.html
python run_sar.py --input library.numbers
```

## Input Format

Your file needs two required columns and any number of optional activity columns:

| Column | Description |
|--------|-------------|
| `Name` | Compound identifier (also accepts: ID, Compound, Cmpd) |
| `HELM` | HELM V2 string (also accepts: HELM_string, Sequence) |
| `IC50_nM` | Optional — any activity column name, pass with `--activity` |

**Reference compound** is auto-detected from the `Name` column: any name containing
"parent", "ref", "wt", "wildtype", or "reference" (case-insensitive) is used.
If none matches, the first row is the reference.

### Example CSV

```csv
Name,HELM,IC50_nM
Parent,"PEPTIDE1{H.F.R.W.[meA].K}$PEPTIDE1,PEPTIDE1,1:R1-6:R2|PEPTIDE1,PEPTIDE1,3:R3-5:R3$$$V2.0",45.2
Analog_1,"PEPTIDE1{H.F.R.W.[meG].K}$PEPTIDE1,PEPTIDE1,1:R1-6:R2|PEPTIDE1,PEPTIDE1,3:R3-5:R3$$$V2.0",12.1
Analog_2,"PEPTIDE1{H.Y.R.W.[meA].K}$PEPTIDE1,PEPTIDE1,1:R1-6:R2|PEPTIDE1,PEPTIDE1,3:R3-5:R3$$$V2.0",89.5
```

## Report Features

- **Sequence alignment** coloured by Zappo chemistry class (aromatic, aliphatic, charged, polar…)
- **Charge/LogP map** — same alignment, coloured red (positive) / blue (negative) / yellow-orange (lipophilic)
- **Sidechain position labels** — main chain positions numbered as integers; sidechain sub-positions as `N.1`, `N.2`, … with green column headers (reading outward from backbone)
- **Scatter plots** — embedded above the table: activity vs lipidation position, and activity vs sidechain composition string
- **Activity column** — colour-coded from potent (green) to inactive (red), sortable
- **MW column** — calculated from residue MW sum minus bond water losses
- **% Identity** — fraction of positions matching the reference, with bar
- **Alignment score** — cosine descriptor similarity after cyclic rotation
- **Changes vs Ref** — explicit list of substitutions per compound
- **Conservation bars** — per-position residue frequency, labelled by position (e.g. `28.1`)
- **Gap columns** — shorter/longer compounds get NW gap cells rather than being truncated
- **Sortable columns** — click any header; reference and consensus rows stay pinned

## How It Works

1. HELM strings are parsed to `HELMObject` (full structural graph with connectivity)
2. Cyclic same-length pairs are rotation-aligned by cosine similarity of 8 RDKit descriptors
3. Sequences are flattened via `get_jpv_flat()`: sidechain monomers are inserted after their backbone attachment position, reading outward from the backbone. Each token carries a position label (`"28"` for main chain, `"28.1"` for first sidechain monomer at position 28, etc.)
4. Flattened sequences are globally NW-aligned; column headers display the position labels
5. MW = Σ(residue MW) − n_bonds × 18.015 (exact peptide bond water losses)
6. Residue colours come from RDKit-calculated LogP and SMARTS-based charge classification

For the full technical detail — and instructions for AI agents — see `AGENTS.md`.

## Dependencies

```
rdkit
numpy
numbers-parser      # Apple Numbers support
openpyxl            # Excel support
```

## Project Structure

```
helm_sar/
├── run_sar.py          ← Entry point
├── AGENTS.md           ← AI agent instructions (full mental model)
├── scripts/
│   ├── helm_engine.py  ← HELMObject
│   ├── helm_parser.py  ← HELM V2 parser
│   ├── helm_alignment.py
│   ├── helm_compare.py
│   ├── residue_colors.py
│   ├── rdkit_bridge.py
│   └── report.py       ← HTML report generator
└── monomer_db/         ← Monomer libraries (standard + custom + CycPeptMPDB)
```
