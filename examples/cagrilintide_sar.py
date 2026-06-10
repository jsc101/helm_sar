#!/usr/bin/env python3
"""Build the amylin/cagrilintide SAR library from Kruse et al. 2021 (J. Med. Chem. 64:11183).

Source tables:
  Table 2 — compound structures (mutations relative to h-amylin + N-terminal protractor)
  Table 3 — in vitro functional potency (EC50, pM) and receptor binding (IC50, pM)

Encoding decisions:
  - Backbone = human amylin (37 aa), disulfide Cys2-Cys7, C-terminal amide (amide implicit,
    as in the repo's PYY example).
  - All compounds are N-terminally lipidated: the protractor is PEPTIDE2 and connects to the
    backbone N-terminus via  PEPTIDE2,PEPTIDE1,1:R1-1:R1  (mirrors the PYY example).
  - Protractor chain is written proximal->distal and capped distally with the C20 diacid (C20d):
        "gGlu"                -> [gGlu].[C20d]
        "gGlu-His-His"        -> [gGlu].H.H.[C20d]
        "gGlu-2xOEG-Arg-Arg"  -> [gGlu].[Ado].[Ado].R.R.[C20d]   (OEG == Ado)
  - Pramlintide is the reference row (h-amylin + 25P,28P,29P, no lipidation).
  - Compound 20 (11Agb) is omitted: the Agb monomer is not in the DB.
  - s-calcitonin (32 aa) is omitted from the aligned table (different length).

Activity: pEC50 = 12 - log10(EC50 in pM); pIC50 = 12 - log10(IC50 in pM).
"""
import csv
import math
from pathlib import Path

# ── Human amylin backbone (1-based positions) ──────────────────────────────────
HAMYLIN = list('KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY')
assert len(HAMYLIN) == 37

PRAMLINTIDE_MODS = ['25P', '28P', '29P']  # h-amylin -> pramlintide

# ── Table 2: (mods relative to h-amylin, protractor linker) ────────────────────
# protractor linker is proximal->distal; C20 diacid is appended distally.
COMPOUNDS = {
    '01':            (['1H', '17H', '25P', '28P', '29P'],                              'gGlu-His-His'),
    '02':            (['1H', '14E', '17H', '18R', '25P', '28P', '29P'],               'gGlu'),
    '03':            (['1R', '17H', '25P', '28P', '29P', '35R'],                      'gGlu-2xOEG-Arg-Arg'),
    '04':            (['1H', '14E', '17H', '18R', '26P'],                             'gGlu'),
    '05':            (['1H', '14E', '17H', '18R', '25E', '29R'],                      'gGlu'),
    '06':            (['1H', '14E', '17H', '18R', '21S', '22S', '25P', '28P', '29P', '31D', '35D'], 'gGlu'),
    '07':            (['1H', '14E', '17H', '25P', '28P', '29P'],                      'gGlu-His-His'),
    '08':            (['1H', '14E', '17R', '25P', '28P', '29P'],                      'gGlu-His-His'),
    '09':            (['1H', '14R', '17H', '25P', '28P', '29P'],                      'gGlu-His-His'),
    '10':            (['1H', '14E', '18R', '25P', '28P', '29P'],                      'gGlu'),
    '11':            (['1H', '14E', '17H', '25P', '28P', '29P'],                      'gGlu-His'),
    '12':            (['1H', '3R', '14E', '17H', '25P', '28P', '29P'],                'gGlu-His-His'),
    '13':            (['14E', '17H', '25P', '28P', '29P'],                            'gGlu-Glu'),
    '14':            (['14E', '17H', '25P', '28P', '29P'],                            'gGlu'),
    '15':            (['des1', '14E', '17H', '25P', '28P', '29P'],                    'gGlu'),
    '16':            (['1E', '14E', '17H', '25P', '28P', '29P'],                      'gGlu'),
    '17':            (['14E', '17H', '18R', '25P', '28P', '29P'],                     'gGlu-2xOEG-D-Arg-D-Arg'),
    '18':            (['14E', '17H', '18R', '21A', '25P', '28P', '29P', '35S'],       'gGlu-2xOEG-D-Arg-D-Arg'),
    '19':            (['1H', '11hArg', '17H', '25P', '28P', '29P'],                   'gGlu-His-His'),
    # '20' omitted: 11Agb monomer unavailable
    '21':            (['14E', '17R', '25P', '28P', '29P'],                            'gGlu'),
    '22':            (['14E', '17H', '25P', '28P', '29P', '37P'],                     'gGlu'),
    '23-cagrilintide': (['14E', '17R', '25P', '28P', '29P', '37P'],                   'gGlu'),
    '24':            (['14E', '17H', '25P', '28P', '29P', '35H'],                     'gGlu'),
    '25':            (['14E', '17H', '25P', '28P', '29P', '35H'],                     'gGlu-Glu'),
    '26':            (['14E', '17R', '25P', '28P', '29P', '35H'],                     'gGlu'),
    '27':            (['1H', '14E', '17H', '25P', '28P', '29P', '37P'],              'gGlu'),
}

# ── Table 3: functional EC50 [pM] and binding IC50 [pM] ────────────────────────
EC50_hAMY3R = {'pram': 5, '01': 74, '02': 83, '03': 76, '04': 107, '05': 77, '06': 166,
               '07': 74, '08': 72, '09': 69, '10': 83, '11': 102, '12': 96, '13': 129,
               '14': 121, '15': 125, '16': 377, '17': 78, '18': 140, '19': 102, '21': 100,
               '22': 61, '23-cagrilintide': 49, '24': 145, '25': 297, '26': 111, '27': 91}
EC50_hCTR = {'pram': 70, '01': 136, '02': 128, '03': 1244, '04': 218, '05': 129, '06': 180,
             '07': 200, '08': 407, '09': 398, '10': 125, '11': 165, '12': 157, '13': 205,
             '14': 107, '15': 117, '16': 657, '17': 608, '18': 1081, '19': 391, '21': 447,
             '22': 52, '23-cagrilintide': 62, '24': 147, '25': 347, '26': 379, '27': 67}
IC50_hAMY3R = {'pram': 114, '01': 216, '02': 137, '03': 213, '05': 818, '07': 723, '10': 499,
               '11': 359, '14': 192, '15': 1348, '17': 52, '19': 1107, '21': 350, '22': 186,
               '23-cagrilintide': 170, '24': 1072, '27': 293}


def apply_mods(mods):
    seq = list(HAMYLIN)
    drop = []
    for m in mods:
        if m.startswith('des'):
            drop.append(int(m[3:]))
            continue
        pos = int(''.join(c for c in m if c.isdigit()))
        res = m[len(str(pos)):]
        seq[pos - 1] = res
    for pos in sorted(drop, reverse=True):
        del seq[pos - 1]
    return seq


def backbone_str(seq):
    return '.'.join(f'[{s}]' if len(s) > 1 else s for s in seq)


def protractor_chain(linker):
    """proximal->distal monomer list, distally capped with C20 diacid."""
    parts = []
    for tok in linker.replace('D-Arg', 'dArg').split('-'):
        if tok == '2xOEG':
            parts += ['Ado', 'Ado']
        elif tok == 'dArg':
            parts.append('dR')
        elif tok == 'Arg':
            parts.append('R')
        elif tok == 'His':
            parts.append('H')
        elif tok == 'Glu':
            parts.append('E')
        elif tok == 'gGlu':
            parts.append('gGlu')
        else:
            raise ValueError(f'unknown linker token: {tok}')
    parts.append('C20d')
    return '.'.join(f'[{p}]' if len(p) > 1 else p for p in parts)


def build_helm(seq, linker=None):
    bb = backbone_str(seq)
    n = len(seq)
    dis = f'PEPTIDE1,PEPTIDE1,2:R3-7:R3'  # Cys2-Cys7 disulfide
    if linker is None:
        return f'PEPTIDE1{{{bb}}}${dis}$$V2.0'
    prot = protractor_chain(linker)
    conns = f'PEPTIDE2,PEPTIDE1,1:R1-1:R1|{dis}'
    return f'PEPTIDE1{{{bb}}}|PEPTIDE2{{{prot}}}${conns}$$V2.0'


def pX(val_pM):
    return round(12 - math.log10(val_pM), 2) if val_pM else ''


# ── Assemble rows ──────────────────────────────────────────────────────────────
rows = []

# reference: pramlintide (no lipidation)
rows.append({
    'Name': 'Pramlintide-ref',
    'HELM': build_helm(apply_mods(PRAMLINTIDE_MODS)),
    'pEC50_hAMY3R': pX(EC50_hAMY3R['pram']),
    'pEC50_hCTR':   pX(EC50_hCTR['pram']),
    'pIC50_hAMY3R': pX(IC50_hAMY3R['pram']),
    'Protractor': 'none',
    'Mods': '+'.join(PRAMLINTIDE_MODS),
})

for cid, (mods, linker) in COMPOUNDS.items():
    seq = apply_mods(mods)
    rows.append({
        'Name': f'Cmpd_{cid}',
        'HELM': build_helm(seq, linker),
        'pEC50_hAMY3R': pX(EC50_hAMY3R.get(cid)),
        'pEC50_hCTR':   pX(EC50_hCTR.get(cid)),
        'pIC50_hAMY3R': pX(IC50_hAMY3R.get(cid)),
        'Protractor': f'{linker}-C20diacid',
        'Mods': '+'.join(mods),
    })

# ── Write CSV ──────────────────────────────────────────────────────────────────
out = Path(__file__).resolve().parent.parent / 'data' / 'cagrilintide_sar.csv'
fields = ['Name', 'HELM', 'pEC50_hAMY3R', 'pEC50_hCTR', 'pIC50_hAMY3R', 'Protractor', 'Mods']
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print(f'Written {len(rows)} rows -> {out}')
for r in rows[:3]:
    print(f"  {r['Name']:18} pEC50_hAMY3R={r['pEC50_hAMY3R']}  {r['HELM'][:60]}...")
