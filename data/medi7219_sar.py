#!/usr/bin/env python3
"""
Step 1 extraction — Pechenov et al. 2021, Sci Rep 11:22521
Development of an orally delivered GLP-1 receptor agonist (MEDI7219)

Compounds:
  GLP-1(7-36)  — native parent (30 AA, C-terminal amide in vivo; encoded as linear)
  Semaglutide  — Aib8, K26 mono-lipidated (gGlu-Ado-Ado-C18d)
  J211         — 12 backbone substitutions, no lipidation
  J229         — J211 + K26 mono-lipidated (gGlu-Ado-Ado-C18d), same protractor as semaglutide
  MEDI7219     — J211 variant, bis-lipidated at K19+K31 with Ado-Ado-C12d

GLP-1 numbering (positions 7-36 → HELM positions 1-30):
  HELM 1=H7, 2=A8, 3=E9, 4=G10, 5=T11, 6=F12, 7=T13, 8=S14, 9=D15, 10=V16,
  11=S17, 12=S18, 13=Y19, 14=L20, 15=E21, 16=G22, 17=Q23, 18=A24, 19=A25, 20=K26,
  21=E27, 22=F28, 23=I29, 24=A30, 25=W31, 26=L32, 27=V33, 28=K34, 29=G35, 30=R36

J211 substitutions (12 total, 8 alpha-methyl AAs):
  A8(2)→Aib, T11(5)→S, F12(6)→aMePhe, S17(11)→aMeSer, Y19(13)→aMePhe,
  Q23(17)→E, K26(20)→aMeLys, F28(22)→aMePhe, W31(25)→aMePhe, L32(26)→V,
  K34(28)→aMeLys, R36(30)→G

Activity data — Table 1 (CHO-hGLP-1R cAMP, EC50 pM):
  0.1% BSA: GLP-1=2.1, Sema=12, J229=132, MEDI7219=3.4, J211=nd
  4.4% HSA: GLP-1=3.3, Sema=2630, J229=nt, MEDI7219=398, J211=nt
"""

import csv
import math
import sys
sys.path.insert(0, '.')

from scripts.helm_parser import HELMParser

# ── Backbone sequences ────────────────────────────────────────────────────────

GLP1_BASE = [
    'H','A','E','G','T','F','T','S','D','V',
    'S','S','Y','L','E','G','Q','A','A','K',
    'E','F','I','A','W','L','V','K','G','R',
]
assert len(GLP1_BASE) == 30


def make_backbone(mods: dict) -> list:
    seq = list(GLP1_BASE)
    for pos, sym in mods.items():
        seq[pos - 1] = sym
    return seq


def fmt_sym(s):
    return f'[{s}]' if len(s) > 1 else s


def make_helm(backbone: list,
              sc1_seq: str = None, sc1_pos: int = None, sc1_rg: str = 'R3',
              sc2_seq: str = None, sc2_pos: int = None, sc2_rg: str = 'R3') -> str:
    p1 = 'PEPTIDE1{' + '.'.join(fmt_sym(s) for s in backbone) + '}'
    if sc1_seq is None:
        return p1 + '$$$$V2.0'

    chains = [p1, f'PEPTIDE2{{{sc1_seq}}}']
    conns  = [f'PEPTIDE2,PEPTIDE1,1:R1-{sc1_pos}:{sc1_rg}']

    if sc2_seq is not None:
        chains.append(f'PEPTIDE3{{{sc2_seq}}}')
        conns.append(f'PEPTIDE3,PEPTIDE1,1:R1-{sc2_pos}:{sc2_rg}')

    return '|'.join(chains) + '$' + '|'.join(conns) + '$$$V2.0'


# ── Lipidation protractors ────────────────────────────────────────────────────

# Semaglutide / J229: C18 diacid-gGlu-2xAdo (same convention as PYY paper)
# Reading outward from backbone: gGlu(1) → Ado(2) → Ado(3) → C18d(4)
SC_C18 = '[gGlu].[Ado].[Ado].[C18d]'

# MEDI7219: C12 diacid-2xAdo (bis-lipidation at K19 and K31)
# Reading outward: Ado(1) → Ado(2) → C12d(3)
SC_C12 = '[Ado].[Ado].[C12d]'


def pec50(ec50_pm):
    if ec50_pm is None:
        return ''
    return f'{-math.log10(ec50_pm * 1e-12):.2f}'


# ── J211 backbone mods ────────────────────────────────────────────────────────

J211_MODS = {
    2:  'Aib',    # A8  → alpha-methyl-Ala (DPP-IV protection)
    5:  'S',      # T11 → Ser (enables aMePhe12)
    6:  'aMePhe', # F12 → alpha-methyl-Phe (chymotrypsin/pepsin/neprilysin)
    11: 'aMeSer', # S17 → alpha-methyl-Ser (elastase)
    13: 'aMePhe', # Y19 → alpha-methyl-Phe (chymotrypsin/pepsin/neprilysin)
    17: 'E',      # Q23 → Glu (solubility, removes deamidation)
    20: 'aMeLys', # K26 → alpha-methyl-Lys (trypsin; lipidated in J229)
    22: 'aMePhe', # F28 → alpha-methyl-Phe (chymotrypsin/pepsin/neprilysin)
    25: 'aMePhe', # W31 → alpha-methyl-Phe (chymotrypsin/pepsin/neprilysin)
    26: 'V',      # L32 → Val (reduce hydrophobicity)
    28: 'aMeLys', # K34 → alpha-methyl-Lys (trypsin)
    30: 'G',      # R36 → Gly (trypsin protection)
}

# MEDI7219 differs from J211 at three positions:
# Y19(13) → K (lipidated with C12 instead of aMePhe)
# W31(25) → K (lipidated with C12 instead of aMePhe)
# K34(28) → E (reduce alpha-methyl count; aMeLys → Glu)
MEDI7219_MODS = dict(J211_MODS)
MEDI7219_MODS[13] = 'K'      # lipidation site 1
MEDI7219_MODS[25] = 'K'      # lipidation site 2
MEDI7219_MODS[28] = 'E'      # K34 → Glu

# ── Build HELM strings ────────────────────────────────────────────────────────

glp1_helm    = make_helm(GLP1_BASE)
sema_helm    = make_helm(make_backbone({2: 'Aib'}), SC_C18, sc1_pos=20)
j211_helm    = make_helm(make_backbone(J211_MODS))
j229_helm    = make_helm(make_backbone(J211_MODS), SC_C18, sc1_pos=20)
medi_helm    = make_helm(make_backbone(MEDI7219_MODS),
                         sc1_seq=SC_C12, sc1_pos=13,
                         sc2_seq=SC_C12, sc2_pos=25)

rows = [
    {'Name': 'GLP-1_ref',   'HELM': glp1_helm,
     'EC50_BSA_pM': 2.1,   'EC50_HSA_pM': 3.3,
     'pEC50_BSA': pec50(2.1), 'pEC50_HSA': pec50(3.3),
     'Lipidation': 'none', 'Alpha_methyl_count': 0},

    {'Name': 'Semaglutide', 'HELM': sema_helm,
     'EC50_BSA_pM': 12,    'EC50_HSA_pM': 2630,
     'pEC50_BSA': pec50(12), 'pEC50_HSA': pec50(2630),
     'Lipidation': 'mono-C18', 'Alpha_methyl_count': 1},

    {'Name': 'J211',        'HELM': j211_helm,
     'EC50_BSA_pM': '',    'EC50_HSA_pM': '',
     'pEC50_BSA': '',      'pEC50_HSA': '',
     'Lipidation': 'none', 'Alpha_methyl_count': 8},

    {'Name': 'J229',        'HELM': j229_helm,
     'EC50_BSA_pM': 132,   'EC50_HSA_pM': '',
     'pEC50_BSA': pec50(132), 'pEC50_HSA': '',
     'Lipidation': 'mono-C18', 'Alpha_methyl_count': 8},

    {'Name': 'MEDI7219',    'HELM': medi_helm,
     'EC50_BSA_pM': 3.4,   'EC50_HSA_pM': 398,
     'pEC50_BSA': pec50(3.4), 'pEC50_HSA': pec50(398),
     'Lipidation': 'bis-C12', 'Alpha_methyl_count': 5},
]

# ── Validate ──────────────────────────────────────────────────────────────────

print("=== Validating HELM strings ===")
for r in rows:
    try:
        obj = HELMParser.parse(r['HELM'])
        chains = obj.data['_chains']
        main   = obj.get_chain()
        lip    = obj.get_lipidation_pos()
        sc     = obj.get_sidechain_string()
        print(f"  {r['Name']:12s}  chains={len(chains)}  main_len={len(main['monomers'])}  "
              f"lip_pos={lip}  sc={sc or 'none'}")
    except Exception as e:
        print(f"  {r['Name']:12s}  ERROR: {e}")
        raise

# ── Write CSV ─────────────────────────────────────────────────────────────────

from pathlib import Path
out = Path(__file__).parent / 'medi7219_sar.csv'
fields = ['Name', 'HELM', 'EC50_BSA_pM', 'pEC50_BSA', 'EC50_HSA_pM', 'pEC50_HSA',
          'Lipidation', 'Alpha_methyl_count']

with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print(f"\nWritten {len(rows)} rows → {out}")
for r in rows:
    print(f"  {r['Name']:12s}  pEC50_BSA={r['pEC50_BSA'] or 'nd':6s}  "
          f"lip={r['Lipidation']}  αMe={r['Alpha_methyl_count']}")
