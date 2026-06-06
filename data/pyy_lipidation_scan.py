#!/usr/bin/env python3
"""Generate PYY3-36 lipidation scan CSV from Ostergaard et al. 2021 (Sci Rep).

Multi-chain HELM encoding:
  PEPTIDE1 — backbone (PYY positions 3-36 = HELM positions 1-34)
  PEPTIDE2 — lipidation protractor (sidechain, pos 1 = closest to backbone)
  Connection — PEPTIDE2,PEPTIDE1,1:R1-{helm_pos}:{rg}

Sidechain monomer order (outward from backbone):
  gGlu -> Ado -> Ado -> C18d  (standard C18-gGlu-2xAdo)
"""

import csv
import math
from pathlib import Path
from typing import Optional

# PYY3-36 backbone: HELM positions 1-34 = PYY positions 3-36
BASE = list('IKPEAPGEDASPEELNRYYASLRHYLNLVTRQRY')
assert len(BASE) == 34

# PYY position → HELM position: helm_pos = pyy_pos - 2
def h(pyy_pos: int) -> int:
    return pyy_pos - 2


def make_helm(backbone_mods: dict, sidechain: Optional[str] = None,
              sc_attach_pos: Optional[int] = None, sc_rgroup: str = 'R3') -> str:
    """
    Build a multi-chain HELM V2 string.

    backbone_mods : {helm_pos: symbol}  — substitutions on PEPTIDE1
    sidechain     : HELM sequence string for PEPTIDE2, e.g. '[gGlu].[Ado].[Ado].[C18d]'
                    Position 1 of PEPTIDE2 is always the monomer closest to the backbone.
    sc_attach_pos : attachment position on PEPTIDE1 (1-based)
    sc_rgroup     : R-group on PEPTIDE1 at sc_attach_pos ('R3' for Lys epsilon-N, 'R1' for N-term)
    """
    seq = BASE.copy()
    for pos, sym in backbone_mods.items():
        seq[pos - 1] = sym
    p1 = 'PEPTIDE1{' + '.'.join(seq) + '}'

    if sidechain is None or sc_attach_pos is None:
        return p1 + '$$$$V2.0'

    p2 = f'PEPTIDE2{{{sidechain}}}'
    conn = f'PEPTIDE2,PEPTIDE1,1:R1-{sc_attach_pos}:{sc_rgroup}'
    return f'{p1}|{p2}${conn}$$V2.0'


def pec50(ec50_nm):
    if ec50_nm is None:
        return ''
    return f'{-math.log10(ec50_nm * 1e-9):.2f}'


rows = []

# Standard lipidation protractor: C18 diacid - gGlu - 2xAdo
# PEPTIDE2 read from backbone outward: gGlu(1) - Ado(2) - Ado(3) - C18d(4)
SC_STANDARD = '[gGlu].[Ado].[Ado].[C18d]'

# ── Reference ────────────────────────────────────────────────────────────────────
rows.append({
    'Name': 'PYY3-36-ref',
    'HELM': make_helm({}),
    'Y2R_EC50_nM': 0.60,
    'HalfLife_h': '',
    'pEC50_Y2R': pec50(0.60),
    'Group': 'ref',
    'PYY_pos': '',
    'Protractor': 'none',
})

# ── Lipidation position scan: analogues 1-33 ────────────────────────────────────
# All use C18 diacid-gGlu-2xAdo protractor
# Data from Table 1 (half-life) and Table 2 (EC50)
scan = [
    # (analogue_no, pyy_lipid_pos, ec50_nM, halflife_h, note)
    # pyy_lipid_pos='Na' = N-terminal (HELM pos 1 R1 attachment)
    (1,  'Na', 10,   17,   'N-alpha acylation'),
    (2,  4,    25,   14,   ''),
    (3,  5,    25,   8.8,  ''),
    (4,  6,    7.9,  7.2,  ''),
    (5,  7,    2.0,  11,   ''),
    (6,  8,    40,   5.9,  ''),
    (7,  9,    13,   12,   ''),
    (8,  10,   1.3,  8.4,  ''),
    (9,  11,   5.0,  11,   ''),
    (10, 12,   16,   15,   ''),
    (11, 13,   32,   13,   ''),
    (12, 14,   4.0,  11,   ''),
    (13, 15,   10,   8.8,  ''),
    (14, 16,   20,   17,   ''),
    (15, 17,   5.0,  39,   ''),
    (16, 18,   7.9,  22,   ''),
    (17, 19,   4.0,  29,   ''),
    (18, 20,   20,   33,   ''),
    (19, 21,   4.0,  34,   ''),
    (20, 22,   2.0,  36,   ''),
    (21, 23,   13,   19,   ''),
    (22, 24,   50,   66,   ''),
    (23, 25,   79,   55,   ''),
    (24, 26,   40,   41,   ''),
    (25, 27,   32,   52,   ''),
    (26, 28,   32,   62,   ''),
    (27, 29,   500,  39,   ''),
    (28, 30,   5.0,  76,   'best position'),
    (29, 31,   79,   75,   ''),
    (30, 32,   50,   49,   ''),
    (31, 33,   None, 56,   ''),
    (32, 34,   None, 30,   ''),
    (33, 35,   None, 67,   ''),
]

for (no, pos, ec50, hl, note) in scan:
    if pos == 'Na':
        # N-terminal: PEPTIDE2 pos 1 (gGlu) connects to PEPTIDE1 pos 1 R1
        helm = make_helm({}, SC_STANDARD, sc_attach_pos=1, sc_rgroup='R1')
        pyy_label = 'Na'
        helm_pos = 1
    else:
        helm_pos = h(pos)
        # Introduce Lys at lipidation position for sidechain attachment via R3
        # (positions that are naturally Lys: pos 4 = HELM 2; all others need K substitution)
        bk_mods = {helm_pos: 'K'} if BASE[helm_pos - 1] != 'K' else {}
        helm = make_helm(bk_mods, SC_STANDARD, sc_attach_pos=helm_pos, sc_rgroup='R3')
        pyy_label = str(pos)

    rows.append({
        'Name': f'Analog_{no:02d}-pos{pyy_label}',
        'HELM': helm,
        'Y2R_EC50_nM': ec50 if ec50 is not None else '',
        'HalfLife_h': hl,
        'pEC50_Y2R': pec50(ec50),
        'Group': 'A',
        'PYY_pos': pyy_label,
        'Protractor': 'C18d-gGlu-2xAdo',
    })

# ── Fatty acid variants at position 30 (analogues 34-43) ────────────────────────
# Lys at HELM pos 28 (PYY pos 30), varied sidechain
fa_variants = [
    # (no, sidechain_helm, fa_label, has_gGlu, ec50, halflife)
    (34, '[gGlu].[Ado].[Ado].[C14d]', 'C14d-gGlu-2xAdo', True,  None, 4),
    (35, '[gGlu].[Ado].[Ado].[C16d]', 'C16d-gGlu-2xAdo', True,  None, 28),
    # analogue 28 already in scan (C18d-gGlu-2xAdo)
    (36, '[gGlu].[Ado].[Ado].[C20d]', 'C20d-gGlu-2xAdo', True,  None, 99),
    (37, '[Ado].[Ado].[C14d]',         'C14d-2xAdo',      False, None, 2),
    (38, '[Ado].[Ado].[C16d]',         'C16d-2xAdo',      False, None, 4),
    (39, '[Ado].[Ado].[C18d]',         'C18d-2xAdo',      False, None, 13),
    (40, '[Ado].[Ado].[C20d]',         'C20d-2xAdo',      False, None, 20),
    (41, '[gGlu].[C18d]',              'C18d-gGlu-noAdo', True,  None, 97),
    (42, '[gGlu].[Ado].[Ado].[Ado].[Ado].[C18d]', 'C18d-gGlu-4xAdo', True, None, 75),
    (43, '[gGlu].[Ado].[Ado].[Ado].[Ado].[Ado].[Ado].[C18d]', 'C18d-gGlu-6xAdo', True, None, 78),
    (44, '[gGlu].[Ado].[Ado].[C16m]',  'C16mono-gGlu-2xAdo', True, None, 0.5),
]

for (no, sc, fa_label, has_gGlu, ec50, hl) in fa_variants:
    helm = make_helm({}, sc, sc_attach_pos=28, sc_rgroup='R3')
    rows.append({
        'Name': f'Analog_{no:02d}-{fa_label}',
        'HELM': helm,
        'Y2R_EC50_nM': ec50 if ec50 is not None else '',
        'HalfLife_h': hl,
        'pEC50_Y2R': pec50(ec50),
        'Group': 'B',
        'PYY_pos': '30',
        'Protractor': fa_label,
    })

# ── Backbone modifications at pos 30 (analogues 45-52) ──────────────────────────
# All have C18 diacid-gGlu-2xAdo at HELM pos 28 (PYY pos 30)
backbone_mods = [
    # (no, backbone_substitutions_helm_pos, sc_attach_pos, ec50, halflife, label)
    # analogue 45: MeArg35 + Lys4 lipidation (Lys4 = HELM 2, native K)
    (45, {h(35): 'meR'},               h(4),  None, 83,  'meArg35-Lys4'),
    # analogue 46: Ala4 mutation + Lys30 lipidation
    (46, {h(4): 'A'},                  28,    None, 84,  'Ala4-Lys30'),
    # analogue 47: Arg4 (K4→R) + Lys30 lipidation
    (47, {h(4): 'R'},                  28,    None, 62,  'Arg4-Lys30'),
    # analogue 48: Asp18 (N18→D) + Lys30 lipidation
    (48, {h(18): 'D'},                 28,    None, 104, 'Asp18-Lys30'),
    # analogue 49: AcI N-term + Ala4 + Lys30 lipidation
    (49, {1: 'AcI', h(4): 'A'},        28,    None, 113, 'AcI-Ala4-Lys30'),
    # analogue 50: AcI N-term + Arg4 + Lys30 lipidation
    (50, {1: 'AcI', h(4): 'R'},        28,    None, 120, 'AcI-Arg4-Lys30'),
    # analogue 51: AcI N-term + Ala4 + Asp18 + Lys30 lipidation
    (51, {1: 'AcI', h(4): 'A', h(18): 'D'}, 28, None, 114, 'AcI-Ala4-Asp18-Lys30'),
    # analogue 52: Arg4 + Gln18 + Lys30 lipidation
    (52, {h(4): 'R', h(18): 'Q'},      28,    4.0,  78,  'Arg4-Gln18-Lys30'),
]

for (no, bk_mods, sc_pos, ec50, hl, label) in backbone_mods:
    # Ensure Lys at attachment position unless it's already K
    backbone = dict(bk_mods)
    if BASE[sc_pos - 1] != 'K' and sc_pos not in backbone:
        backbone[sc_pos] = 'K'
    helm = make_helm(backbone, SC_STANDARD, sc_attach_pos=sc_pos, sc_rgroup='R3')
    rows.append({
        'Name': f'Analog_{no:02d}-{label}',
        'HELM': helm,
        'Y2R_EC50_nM': ec50 if ec50 is not None else '',
        'HalfLife_h': hl,
        'pEC50_Y2R': pec50(ec50),
        'Group': 'C',
        'PYY_pos': str(sc_pos + 2),  # convert back to PYY position
        'Protractor': 'C18d-gGlu-2xAdo',
    })

# ── Write CSV ────────────────────────────────────────────────────────────────────
out = Path(__file__).parent / 'pyy_lipidation_scan.csv'
fields = ['Name', 'HELM', 'Y2R_EC50_nM', 'HalfLife_h', 'pEC50_Y2R', 'Group', 'PYY_pos', 'Protractor']

with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print(f'Written {len(rows)} rows to {out}')
for r in rows[:3]:
    print(f"  {r['Name']}: {r['HELM'][:100]}")
