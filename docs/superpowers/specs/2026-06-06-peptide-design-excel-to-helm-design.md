# Peptide Design: Excel → HELM

**Date:** 2026-06-06
**Status:** Design approved, pending spec review

## Problem

Medicinal chemists design lipidated/conjugated peptide analogs in a flat,
position-by-position "JPV" mindset — one column per main-chain position, the way
they already think about a sequence. But the flat view throws away exactly the
information HELM needs to be valid:

1. **Which** main-chain position a sidechain attaches to (the labeling site)
2. **Through which R-groups** each bond forms (head-to-tail R2→R1 for plain
   peptides, but R3-R1 / R4-R3 / etc. once a CHEM monomer like a click triazole
   is involved)
3. How to keep **multiple independent** sidechains separate

A single "parent HELM" template can't anchor this — it can't represent a *new*
sidechain the parent lacks, can't relocate an attachment site, and can't carry
multiple branches. The connectivity must live in the spreadsheet structure
itself.

The reference stress-test is the PYY library (`data/pyy_lipidation_scan.csv`,
52 analogs), which exercises every wrinkle: an attachment-site scan (pos 4→35),
an N-terminal attach via R1 (not R3), a protractor-composition scan
(C14/16/18/20 diacids, mono-acid, 0–6× Ado, ±gGlu), combinatorial backbones
that mutate several positions at once, and a monomer (`C16m`) not yet in the
monomer DB.

## Format: one row per compound

The user authors and *sees* the whole molecule on a single line, left→right.

### Columns

```
| main chain pos 1..N | Site | b1 | SC1 | b2 | SC2 | b3 | SC3 | ... |
```

- **Main-chain block** — one column per HELM residue position. A full,
  editable row per compound (fill-down is a convenience, NOT a constraint):
  the PYY group-C analogs mutate positions 1, 2, 15, 28 simultaneously, so the
  main chain is *not* globally constant.
- **Site** — the HELM residue *index* the sidechain attaches to. This is the
  index of the cell where the attaching residue sits, derived from the
  main-chain block — NOT the paper's own position numbering (which is frequently
  offset; "pos 4" in the PYY paper is HELM residue 2).
- **Alternating (bond, monomer) pairs** — `b1 SC1 b2 SC2 … bN SCN`.
  - `SC_k` = monomer code (proximal→distal, i.e. reading outward from backbone).
  - `b_k` = the bond *into* `SC_k`, written `proximalRg-distalRg` (two R-groups,
    always). For `b1` the proximal side is the backbone residue's R-group; for
    `b_k` (k>1) it is the prior monomer's R-group. The distal side is `SC_k`'s
    facing R-group.
  - N monomers → N bonds. The **final monomer has no trailing bond** — its
    distal terminus is free (e.g. C18d's free distal COOH).

A compound with no sidechain simply leaves Site and the SC/bond columns blank.
A compound with a second sidechain uses a second block of (Site, b/SC pairs) —
the row format repeats; nothing else changes.

### Worked examples

Click-chemistry branch (non-default bonds, mixed polymer types):

```
| Site | b1    | SC1 | b2    | SC2        | b3    | SC3 | b4    | SC4 | b5    | SC5  |
| 37   | R3-R1 | Hpg | R3-R1 | Triazole14 | R4-R3 | Aha | R2-R1 | A   | R2-R1 | C18d |
```

Reads: backbone pos37 R3 ─R3-R1─ Hpg ─R3-R1─ Triazole14 ─R4-R3─ Aha ─R2-R1─
Ala ─R2-R1─ C18d(free COOH).

Plain lipidation branch (all default head-to-tail):

```
| Site | b1    | SC1  | b2    | SC2 | b3    | SC3 | b4    | SC4  |
| 30   | R3-R1 | gGlu | R2-R1 | Ado | R2-R1 | Ado | R2-R1 | C18d |
```

## Generation rule: collapse vs expand

For each sidechain branch the generator inspects its bonds and monomer polymer
types:

- **Auto-collapse** — when *every* bond is the default head-to-tail
  (`proximal R2 - distal R1`, with the first bond being `backboneRg-R1`) AND all
  monomers share one polymer type → emit a single chain:
  `PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}` plus one connection line
  `PEPTIDE2,PEPTIDE1,1:R1-30:R3`.

- **Explicit-bond expansion** — when *any* bond is non-default OR the branch
  mixes polymer types (e.g. PEPTIDE Hpg/Aha + CHEM Triazole14) → emit each
  segment as its own chain and connect with explicit bond lines:

  ```
  PEPTIDE1{...}|PEPTIDE2{[Hpg]}|CHEM1{[Triazole14]}|PEPTIDE3{[Aha]}|PEPTIDE4{A.[C18d]}$
  PEPTIDE2,PEPTIDE1,1:R1-37:R3|CHEM1,PEPTIDE2,1:R1-1:R3|PEPTIDE3,CHEM1,1:R3-1:R4|PEPTIDE4,PEPTIDE3,1:R1-1:R2$$V2.0
  ```

The author writes the *same* (bond, monomer) pairs either way; the tool decides
the representation. Default head-to-tail runs within an expansion may still be
collapsed into a multi-monomer chain where polymer type is continuous.

## Two deliverables (shared core)

1. **Ready-to-go Excel** (`.xlsm`):
   - **Design** sheet — the one-row-per-compound layout above.
   - **HELM** sheet — generated output, one HELM per compound + validation status.
   - **Example** sheet — pre-filled with a known library (PYY and/or MEDI7219)
     so the format is self-documenting.
   - A **VBA macro** ("Build HELM" button) runs generation. A macro, not a pure
     cell formula — assembling multi-chain HELM with connection lines is past
     what a sane single formula can do.

2. **Standalone directory** (`peptide_design/`) — a small Python module that
   reads the *same* two-region layout from `.xlsx`/CSV and emits identical HELM,
   for deployment as CLI/script/service later. The generation logic is shared in
   spirit between the two front-ends (the VBA and Python implement the same
   collapse/expand algorithm against the same monomer DB).

## Validation: separate pass

Generation and validation are decoupled (explicit user preference). After HELM
strings are produced, a separate validator pass runs each through the existing
`HELMParser` + monomer-DB lookup and writes a status per compound:

- unknown monomer (e.g. PYY's `C16m`) → flagged "unknown monomer — add to DB",
  not silently emitted as broken HELM
- malformed connectivity / parse failure → flagged with the parse error
- valid → OK (optionally with computed MW for a sanity check)

The validator never edits the design; it only reports.

## Out of scope

- Round-tripping arbitrary existing HELM back into the design grid (this is a
  design/authoring tool, not a general HELM editor).
- Branches on branches (a sidechain off a sidechain). Current scope is linear
  branches off the main chain; revisit if a real library needs it.
- Automatic monomer creation. Unknown monomers are flagged for the user to add
  to `custom_monomers.json` by hand.

## Open items for the plan

- Exact column-naming / header convention on the Design sheet (how repeated
  second-sidechain blocks are headed).
- Whether the Python core and VBA share a serialized monomer DB or each read
  `monomer_db/` directly.
- Where `peptide_design/` lives relative to `helm_sar` (sibling dir vs subdir).
