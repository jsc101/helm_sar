#!/usr/bin/env python3
"""
report.py — Sequence alignment + property report for a HELM peptide library.

Called by run_sar.py. Can also be used as a library:

    from scripts.report import build_data, build_html

Public API
----------
build_data(pairs, db)               → (ref_obj, rows)
build_html(ref_obj, rows, db,
           activity_map, act_label) → HTML string
"""

from __future__ import annotations

import sys, os, argparse, html
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import logging
logging.disable(logging.WARNING)

from scripts.helm_parser import HELMParser
from scripts.helm_alignment import HELMAlignment
from scripts.residue_colors import residue_color, text_color, LEGEND_CATEGORIES, category_of
from monomer_db.monomer_db import MonomerDB


# ── MW calculation ──────────────────────────────────────────────────────────────

def calc_mw(obj) -> float | None:
    desc = obj.all_monomer_descriptors()
    if not desc:
        return None
    if any(not v['in_db'] for v in desc.values()):
        return None
    total = sum(v['descriptors']['MW'] for v in desc.values())
    n = len(desc)
    n_explicit = len(obj.data.get('connectivity_graph', []))
    return total - ((n - 1) + n_explicit) * 18.015


# ── Charge/polarity color scheme ────────────────────────────────────────────────
#
# positive at pH 7  → red   (K, R, H, Dab, Orn, …)
# negative at pH 7  → blue  (D, E, gGlu, …)
# polar / neutral   → white (S, T, N, Q, Hyp, C, G, P, …)
# aliphatic/aromatic→ yellow → orange scaled by LogP

def _charge_logp_color(sym: str, desc: dict, db=None) -> str:
    # Classification is structural, not electrochemical:
    #   'negative' = residue has >1 COOH (sidechain acid + backbone)
    #   'positive' = residue has extra amine(s), guanidinium, or imidazole
    # His is marked positive even though its imidazole pKa≈6 means ~50% protonated at pH 7.
    # If you wanted actual pH-7 charge, you'd need RDKit's MoleculeEnumerator with Epik or a
    # Henderson-Hasselbalch step — significantly more complex. Keeping structural heuristic as-is.
    cat = category_of(sym, db)
    logp = desc.get('LogP', 0.0) if desc else 0.0

    if cat == 'positive':
        return '#e62020'
    if cat == 'negative':
        return '#1e50e6'
    if cat in ('polar', 'cysteine', 'conformational'):
        return '#f8f8f8'
    if cat in ('aliphatic', 'aromatic'):
        # yellow (logp≈0) → orange (logp≈4+), driven by RDKit Crippen LogP
        t = max(0.0, min(1.0, logp / 4.0))
        g = int(220 - t * 130)   # green channel 220→90
        return f'#ff{g:02x}00'
    return '#d8d8d8'


# ── Data helpers ────────────────────────────────────────────────────────────────

def _main_chain_syms(obj) -> list[str]:
    """Main chain symbols only — used for NW alignment (no sidechain monomers)."""
    if not hasattr(obj, 'get_jpv_flat'):
        return []
    return [t['symbol'] for t in obj.get_jpv_flat() if t['is_main']]


def _main_chain_labels(obj) -> list[str]:
    """Position labels for main chain only, parallel to _main_chain_syms."""
    if not hasattr(obj, 'get_jpv_flat'):
        return []
    return [t['label'] for t in obj.get_jpv_flat() if t['is_main']]


def _token_descs(obj) -> tuple[list[dict], dict]:
    """
    Return (main_descs, sc_desc_map) where:
      main_descs   — per-position descriptor dict for each main chain token (in order)
      sc_desc_map  — {main_pos: [desc, desc, ...]} sidechain descriptor lists (outward order)
    """
    if not hasattr(obj, 'get_jpv_flat'):
        return [], {}

    chains = {c['chain_id']: c for c in obj.data.get('_chains', [])}
    primary = obj.get_chain()
    if primary is None:
        return [], {}
    pid = primary['chain_id']

    all_chain_descs = {cid: obj.all_monomer_descriptors(cid) for cid in chains}
    primary_descs = all_chain_descs.get(pid, {})

    # Build sidechain position map (main_pos → list of (sc_chain_id, sc_helm_pos) outward)
    conn_graph = obj.data.get('connectivity_graph', [])
    sc_order: dict[int, list] = {}
    for conn in conn_graph:
        fc, fp = conn['from_chain'], conn['from_pos']
        tc, tp = conn['to_chain'],   conn['to_pos']
        if fc == pid and tc != pid:
            sc_id, main_attach, sc_attach = tc, fp, tp
        elif tc == pid and fc != pid:
            sc_id, main_attach, sc_attach = fc, tp, fp
        else:
            continue
        sc = chains.get(sc_id)
        if sc is None:
            continue
        n = len(sc['monomers'])
        positions = list(range(1, n + 1))
        if sc_attach == n:
            positions = list(reversed(positions))
        elif sc_attach != 1:
            idx = sc_attach - 1
            positions = positions[idx:] + list(reversed(positions[:idx]))
        sc_order[main_attach] = [(sc_id, p) for p in positions]

    main_descs = []
    sc_desc_map: dict[int, list] = {}
    for token in obj.get_jpv_flat():
        mp = token['main_pos']
        if token['is_main']:
            d = primary_descs.get(mp, {})
            main_descs.append(d.get('descriptors', {}))
        else:
            k = token['sc_pos'] - 1
            entry_list = sc_order.get(mp, [])
            if k < len(entry_list):
                sc_cid, sc_p = entry_list[k]
                d = all_chain_descs.get(sc_cid, {}).get(sc_p, {})
                sc_desc_map.setdefault(mp, []).append(d.get('descriptors', {}))
            else:
                sc_desc_map.setdefault(mp, []).append({})

    return main_descs, sc_desc_map


def _changes_vs_ref(syms: list[str], ref_syms: list[str]) -> str:
    changes = [
        f"P{i+1}:{ref_syms[i]}→{s}"
        for i, s in enumerate(syms)
        if i < len(ref_syms) and s != ref_syms[i]
    ]
    return ', '.join(changes) if changes else '—'


# ── Numbers reader ──────────────────────────────────────────────────────────────

def read_numbers_helms(path: str) -> list[tuple[str, str]]:
    from numbers_parser import Document
    doc = Document(path)
    sheets = {s.name: s for s in doc.sheets}
    if 'HELM_Builder' not in sheets:
        raise ValueError("No 'HELM_Builder' sheet found")
    table = sheets['HELM_Builder'].tables[0]
    pairs = []
    for i, row in enumerate(table.rows()):
        if i == 0:
            continue
        name = str(row[0].value).strip() if row[0].value else ''
        helm = str(row[1].value).strip() if len(row) > 1 and row[1].value else ''
        if helm and 'PEPTIDE' in helm:
            pairs.append((name, helm))
    return pairs


# ── Alignment pipeline ──────────────────────────────────────────────────────────

def build_data(pairs: list[tuple[str, str]], db: MonomerDB) -> tuple:
    ref_name, ref_helm = next(
        ((n, h) for n, h in pairs if 'parent' in n.lower()),
        pairs[0]
    )
    ref_obj = HELMParser.parse(ref_helm)
    ref_mw  = calc_mw(ref_obj)
    ref_row = {
        'id': ref_name, 'helm': ref_helm, 'obj': ref_obj,
        'mw': ref_mw, 'score': 1.0, 'is_ref': True,
    }

    aligner = HELMAlignment(ref_obj)
    rows = [ref_row]

    for name, helm in pairs:
        if name.lower() == ref_name.lower():
            continue
        try:
            obj = HELMParser.parse(helm)
            if len(obj.positions()) != len(ref_obj.positions()):
                rows.append({'id': name, 'helm': helm, 'obj': obj,
                             'mw': calc_mw(obj), 'score': None, 'is_ref': False})
                continue
            result = aligner.align(obj)
            rows.append({
                'id':    name,
                'helm':  helm,
                'obj':   result.rotated,
                'mw':    calc_mw(obj),
                'score': result.best_score,
                'is_ref': False,
            })
        except Exception as e:
            rows.append({'id': name, 'helm': helm, 'obj': None,
                         'mw': None, 'score': None, 'is_ref': False, 'err': str(e)})

    return ref_obj, rows


# ── Embedded scatter plots ──────────────────────────────────────────────────────

def _make_scatter_plots(rows: list, activity_map: dict | None, act_label: str) -> str:
    """
    Generate two embedded scatter plots as base64 PNG.
    Plot 1: Lipidation position (main chain 1-based, PYY numbering) vs activity.
    Plot 2: Sidechain string vs activity (x-axis rotated 90°).
    Returns HTML string with <img> tags, or '' if no activity data.
    """
    if not activity_map:
        return ''

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import base64, io
    except ImportError:
        return ''

    def _to_b64(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='white')
        buf.seek(0)
        data = base64.b64encode(buf.read()).decode()
        plt.close(fig)
        return data

    def _pyy_pos(helm_pos):
        return helm_pos + 2  # HELM 1-based → PYY 3-36 numbering

    lip_pts: list[tuple] = []
    sc_pts:  list[tuple] = []

    for r in rows:
        if r.get('obj') is None or r.get('is_ref'):
            continue
        name = r['id']
        val  = activity_map.get(name)
        if val is None or not isinstance(val, (int, float)):
            continue
        obj = r['obj']
        lp  = obj.get_lipidation_pos() if hasattr(obj, 'get_lipidation_pos') else None
        sc  = obj.get_sidechain_string() if hasattr(obj, 'get_sidechain_string') else ''
        if lp is not None:
            lip_pts.append((_pyy_pos(lp), float(val), name))
        if sc:
            sc_pts.append((sc, float(val), name))

    parts = []

    # ── Plot 1: position vs activity ──────────────────────────────────────
    if lip_pts:
        fig, ax = plt.subplots(figsize=(11, 4))
        xs = [p[0] for p in lip_pts]
        ys = [p[1] for p in lip_pts]
        ax.scatter(xs, ys, s=65, alpha=0.85, color='#2255bb',
                   edgecolors='white', linewidths=0.6, zorder=5)
        ax.set_xlabel('Lipidation position (PYY numbering)', fontsize=10)
        ax.set_ylabel(act_label, fontsize=10)
        ax.set_title('Activity vs. lipidation position', fontsize=11, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.yaxis.grid(True, alpha=0.3, linestyle=':')
        ax.set_axisbelow(True)
        b64 = _to_b64(fig)
        parts.append(
            f'<h2>Activity vs. Lipidation Position</h2>'
            f'<p class="subtitle">Each point = one analogue. X = lipidation site on PYY backbone.</p>'
            f'<img src="data:image/png;base64,{b64}" style="max-width:100%;margin:8px 0 20px">'
        )

    # ── Plot 2: sidechain string vs activity ──────────────────────────────
    if sc_pts:
        sc_strings = sorted(set(p[0] for p in sc_pts))
        sc_idx = {s: i for i, s in enumerate(sc_strings)}
        xs2 = [sc_idx[p[0]] for p in sc_pts]
        ys2 = [p[1] for p in sc_pts]

        fig_h = max(4.5, len(sc_strings) * 0.35)
        fig, ax = plt.subplots(figsize=(8, fig_h))
        ax.scatter(xs2, ys2, s=75, alpha=0.85, color='#2CA02C',
                   edgecolors='white', linewidths=0.6, zorder=5)
        ax.set_xticks(range(len(sc_strings)))
        ax.set_xticklabels(sc_strings, rotation=90, ha='center', fontsize=8)
        ax.set_ylabel(act_label, fontsize=10)
        ax.set_title('Activity vs. sidechain composition', fontsize=11, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.yaxis.grid(True, alpha=0.3, linestyle=':')
        ax.set_axisbelow(True)
        plt.tight_layout()
        b64 = _to_b64(fig)
        parts.append(
            f'<h2>Activity vs. Sidechain Composition</h2>'
            f'<p class="subtitle">Sidechain string = monomers from backbone outward (dash-separated). '
            f'Each point = one analogue.</p>'
            f'<img src="data:image/png;base64,{b64}" style="max-width:100%;margin:8px 0 20px">'
        )

    return '\n'.join(parts)


# ── MW color (for MW column badge) ─────────────────────────────────────────────

def _mw_color(val: float | None, vmin: float, vmax: float) -> str:
    if val is None:
        return '#e8e8e8'
    t = max(0.0, min(1.0, (val - vmin) / max(vmax - vmin, 1e-6)))
    if t < 0.5:
        s = t * 2
        g = int(255 - s * 100); b = int(255 - s * 255)
        return f'#ff{g:02x}{b:02x}'
    else:
        s = (t - 0.5) * 2
        g = int(155 - s * 155)
        return f'#ff{g:02x}00'


# ── HTML generation ─────────────────────────────────────────────────────────────

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Courier New', monospace; background: #f4f6f9; color: #222; padding: 20px; }
h1 { font-size: 1.3rem; font-weight: 700; margin-bottom: 4px; color: #1a2a4a; }
.subtitle { font-size: 0.8rem; color: #555; margin-bottom: 22px; }
h2 { font-size: 1.0rem; font-weight: 700; margin: 22px 0 8px; color: #2a3a5a;
     border-left: 4px solid #4477cc; padding-left: 8px; }

/* ── Shared table shell ───────────────────────── */
.aln-wrap { overflow-x: auto; }
table.aln { border-collapse: collapse; font-size: 0.72rem; table-layout: auto; }
table.aln th, table.aln td { border: 1px solid #ccc; text-align: center; white-space: nowrap; }

/* fixed left columns */
table.aln th.row-hdr { text-align: right; padding: 2px 6px; font-weight: 600; min-width: 120px; }
table.aln td.row-hdr { text-align: right; padding: 2px 6px; min-width: 120px;
                        font-size: 0.68rem; color: #333; }
table.aln td.row-hdr.ref-hdr { font-weight: 700; color: #003388; }
table.aln td.row-hdr.con-hdr { font-weight: 700; color: #006600; font-style: italic; }
table.aln .score-col { min-width: 74px; padding: 2px 4px; }
table.aln .id-col    { min-width: 68px; padding: 2px 4px; }
table.aln .mw-col    { min-width: 76px; padding: 2px 6px; font-weight: 700; }
table.aln .chg-col   { min-width: 160px; max-width: 260px; padding: 2px 6px;
                        text-align: left; font-size: 0.65rem; color: #444; white-space: normal; }
table.aln .act-col   { min-width: 80px; padding: 2px 6px; font-weight: 600;
                        text-align: center; font-size: 0.7rem; }

/* residue cells */
.cell { width: 58px; height: 26px; padding: 0; position: relative; cursor: default; }
.gap-cell { width: 58px; height: 26px; background: #ebebeb; border-color: #ddd !important; }
.cell span { display: block; width: 100%; height: 100%; line-height: 26px;
             font-size: 0.7rem; font-weight: 600; overflow: hidden; text-overflow: ellipsis; }
.cell:hover::after {
  content: attr(data-tip); position: absolute; bottom: 110%; left: 50%;
  transform: translateX(-50%); background: #222; color: #fff;
  padding: 2px 6px; border-radius: 4px; font-size: 0.65rem; white-space: nowrap;
  z-index: 100; pointer-events: none;
}
tr.ref-row td { border-bottom: 2px solid #336; }
tr.sep-row td { border-top: 2px dashed #555; }

/* sortable header */
th.sortable { cursor: pointer; user-select: none; }
th.sortable:hover { background: #d0d8ee; }
th.sort-asc::after  { content: ' ▲'; font-size: 0.6rem; }
th.sort-desc::after { content: ' ▼'; font-size: 0.6rem; }
.pos-hdr { background: #e8eef6; font-weight: 700; font-size: 0.75rem;
           padding: 3px 2px; color: #1a3a6a; }

/* score bar */
.score-bar-wrap { display: flex; align-items: center; gap: 4px; height: 22px; }
.score-bar { height: 12px; border-radius: 3px; }

/* conservation */
.cons-wrap { display: flex; gap: 1px; margin: 6px 0 14px; overflow-x: auto; }
.cons-col { display: flex; flex-direction: column; align-items: center; width: 58px; flex-shrink: 0; }
.cons-bar-track { width: 46px; height: 40px; background: #e0e4ea; border-radius: 3px;
                  display: flex; align-items: flex-end; overflow: hidden; border: 1px solid #ccc; }
.cons-bar-fill { width: 100%; background: #2255bb; border-radius: 2px 2px 0 0; }
.cons-label { font-size: 0.62rem; color: #444; margin-top: 2px; }
.cons-pos { font-size: 0.6rem; color: #888; }

/* legend */
.legend-wrap { display: flex; flex-wrap: wrap; gap: 8px; margin: 6px 0 16px; }
.legend-item { display: flex; align-items: center; gap: 4px; font-size: 0.68rem; }
.legend-swatch { width: 16px; height: 16px; border-radius: 3px; border: 1px solid #aaa; flex-shrink: 0; }
.charge-legend { display: flex; flex-wrap: wrap; gap: 12px; margin: 6px 0 14px; font-size: 0.68rem; }
.charge-legend-item { display: flex; align-items: center; gap: 5px; }
.charge-swatch { width: 18px; height: 18px; border-radius: 3px; border: 1px solid #aaa; }
"""

_JS = """
// ── Column hover highlighting ───────────────────────
function highlight(cls) {
  document.querySelectorAll('.' + cls).forEach(el => {
    el.style.outline = '2px solid #ff8800'; el.style.zIndex = '10';
  });
}
function unhighlight(cls) {
  document.querySelectorAll('.' + cls).forEach(el => {
    el.style.outline = ''; el.style.zIndex = '';
  });
}

// ── Table sorting ───────────────────────────────────
var _sortState = {};   // tableId → {col, dir}

function sortTable(tableId, colIdx, type) {
  var tbl   = document.getElementById(tableId);
  var tbody = tbl.tBodies[0];
  var rows  = Array.from(tbody.rows);

  // separate ref row (first) and consensus row (last, class sep-row)
  var pinned_top = rows.filter(r => r.classList.contains('ref-row'));
  var pinned_bot = rows.filter(r => r.classList.contains('sep-row'));
  var sortable   = rows.filter(r => !r.classList.contains('ref-row') && !r.classList.contains('sep-row'));

  var st  = _sortState[tableId] || {col: -1, dir: 1};
  var dir = (st.col === colIdx) ? -st.dir : 1;
  _sortState[tableId] = {col: colIdx, dir: dir};

  sortable.sort(function(a, b) {
    var av = _cellVal(a, colIdx, type);
    var bv = _cellVal(b, colIdx, type);
    if (av < bv) return -dir;
    if (av > bv) return  dir;
    return 0;
  });

  // Clear sort indicators on this table's headers
  tbl.querySelectorAll('th.sortable').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
  });
  var th = tbl.querySelectorAll('th.sortable')[
    Array.from(tbl.querySelectorAll('th')).indexOf(
      tbl.querySelectorAll('th')[colIdx])];
  // simpler: mark by data-col attribute
  tbl.querySelectorAll('th[data-col="' + colIdx + '"]').forEach(th => {
    th.classList.add(dir === 1 ? 'sort-asc' : 'sort-desc');
  });

  pinned_top.concat(sortable).concat(pinned_bot).forEach(r => tbody.appendChild(r));
}

function _cellVal(row, colIdx, type) {
  var td = row.cells[colIdx];
  if (!td) return '';
  var raw = (td.dataset.sort !== undefined) ? td.dataset.sort : td.innerText.trim();
  if (type === 'num') { var n = parseFloat(raw); return isNaN(n) ? -Infinity : n; }
  return raw.toLowerCase();
}
"""


# ── Needleman-Wunsch sequence alignment ────────────────────────────────────────

def _nw_align(ref_syms: list, query_syms: list,
              match: float = 2.0, mismatch: float = 0.0,
              gap: float = -1.5) -> tuple:
    """
    Global Needleman-Wunsch alignment of two monomer sequences.
    Returns (aligned_ref, aligned_query) — lists where None marks a gap.
    Scores: same symbol = match, different = mismatch, insertion/deletion = gap.
    """
    n, m = len(ref_syms), len(query_syms)
    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1): dp[i][0] = i * gap
    for j in range(m + 1): dp[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if ref_syms[i-1] == query_syms[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap)
    ar, aq = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if ref_syms[i-1] == query_syms[j-1] else mismatch
            if abs(dp[i][j] - (dp[i-1][j-1] + s)) < 1e-9:
                ar.append(ref_syms[i-1]); aq.append(query_syms[j-1])
                i -= 1; j -= 1; continue
        if i > 0 and (j == 0 or abs(dp[i][j] - (dp[i-1][j] + gap)) < 1e-9):
            ar.append(ref_syms[i-1]); aq.append(None); i -= 1
        else:
            ar.append(None); aq.append(query_syms[j-1]); j -= 1
    ar.reverse(); aq.reverse()
    return ar, aq


def _align_group_rows(grp_ref_syms: list, rows: list,
                      grp_ref_labels: list | None = None) -> tuple:
    """
    Align rows to the reference backbone, then insert sidechain columns structurally.

    NW operates on **main chain symbols only**. Sidechain monomers are inserted after
    their anchor main-chain position based on HELM connectivity — never repositioned by
    NW. This prevents sidechain tokens from bleeding into adjacent main-chain columns
    when the reference has no sidechain.

    Column types in master_ref / master_labels:
      master_ref[c] not None                    → main chain position
      master_ref[c] None, master_labels[c] '.' → sidechain slot (green header)
      master_ref[c] None, master_labels[c] None → NW gap column (grey '·' header)

    Returns (n_cols, master_ref, master_labels, aligned_syms_list, aligned_desc_list).
    """
    if grp_ref_labels is None:
        grp_ref_labels = [str(i + 1) for i in range(len(grp_ref_syms))]

    if not rows:
        return len(grp_ref_syms), list(grp_ref_syms), list(grp_ref_labels), [], []

    ref_n = len(grp_ref_syms)

    # ── NW on main chain only ──────────────────────────────────────────────
    pairwise = []
    for r in rows:
        q = _main_chain_syms(r['obj']) if r['obj'] is not None else []
        pairwise.append(_nw_align(grp_ref_syms, q) if q != grp_ref_syms
                        else (list(grp_ref_syms), list(q)))

    gaps_before = [0] * (ref_n + 1)
    for aref, _ in pairwise:
        p, run = 0, 0
        for sym in aref:
            if sym is None:
                run += 1
            else:
                gaps_before[p] = max(gaps_before[p], run)
                p += 1; run = 0
        gaps_before[ref_n] = max(gaps_before[ref_n], run)

    # ── Max sidechain depth at each ref label ──────────────────────────────
    sc_len: dict[str, int] = {}
    for r in rows:
        if r['obj'] is None or not hasattr(r['obj'], 'get_jpv_flat'):
            continue
        for tok in r['obj'].get_jpv_flat():
            if not tok['is_main']:
                lbl = str(tok['main_pos'])
                sc_len[lbl] = max(sc_len.get(lbl, 0), tok['sc_pos'])

    # ── Master column layout ───────────────────────────────────────────────
    master_ref, master_labels = [], []
    for p in range(ref_n):
        master_ref.extend([None] * gaps_before[p])
        master_labels.extend([None] * gaps_before[p])
        master_ref.append(grp_ref_syms[p])
        master_labels.append(grp_ref_labels[p])
        n_sc = sc_len.get(grp_ref_labels[p], 0)
        for k in range(1, n_sc + 1):
            master_ref.append(None)
            master_labels.append(f'{grp_ref_labels[p]}.{k}')
    master_ref.extend([None] * gaps_before[ref_n])
    master_labels.extend([None] * gaps_before[ref_n])
    n_cols = len(master_ref)

    # ── Per-row projection ─────────────────────────────────────────────────
    aligned_syms_list, aligned_desc_list = [], []
    for (aref, aqry_main), r in zip(pairwise, rows):
        main_descs, sc_desc_map = _token_descs(r['obj']) if r['obj'] is not None else ([], {})

        sc_syms: dict[int, list] = {}
        if r['obj'] is not None and hasattr(r['obj'], 'get_jpv_flat'):
            for tok in r['obj'].get_jpv_flat():
                if not tok['is_main']:
                    sc_syms.setdefault(tok['main_pos'], []).append(tok['symbol'])

        result, d_list = [], []
        ai = 0; main_qi = 0

        for p in range(ref_n):
            n_g, n_used = gaps_before[p], 0
            while ai < len(aref) and aref[ai] is None and n_used < n_g:
                sym = aqry_main[ai]
                result.append(sym)
                d_list.append(main_descs[main_qi] if sym is not None and main_qi < len(main_descs) else {})
                if sym is not None: main_qi += 1
                ai += 1; n_used += 1
            result.extend([None] * (n_g - n_used))
            d_list.extend([{}]   * (n_g - n_used))

            if ai < len(aref) and aref[ai] == grp_ref_syms[p]:
                sym = aqry_main[ai]
                result.append(sym)
                d_list.append(main_descs[main_qi] if main_qi < len(main_descs) else {})
                if sym is not None: main_qi += 1
                ai += 1
            else:
                result.append(None); d_list.append({})

            n_sc = sc_len.get(grp_ref_labels[p], 0)
            try:    mp = int(grp_ref_labels[p])
            except: mp = None
            syms_here  = sc_syms.get(mp, [])    if mp is not None else []
            descs_here = sc_desc_map.get(mp, []) if mp is not None else []
            for k in range(n_sc):
                result.append(syms_here[k]  if k < len(syms_here)  else None)
                d_list.append(descs_here[k] if k < len(descs_here) else {})

        n_t, n_u = gaps_before[ref_n], 0
        while ai < len(aref) and aref[ai] is None and n_u < n_t:
            sym = aqry_main[ai]
            result.append(sym)
            d_list.append(main_descs[main_qi] if sym is not None and main_qi < len(main_descs) else {})
            if sym is not None: main_qi += 1
            ai += 1; n_u += 1
        result.extend([None] * (n_t - n_u))
        d_list.extend([{}]   * (n_t - n_u))

        aligned_syms_list.append(result)
        aligned_desc_list.append(d_list)

    return n_cols, master_ref, master_labels, aligned_syms_list, aligned_desc_list


def build_html(ref_obj, rows: list, db=None,
               activity_map: dict | None = None,
               act_label: str = 'Activity',
               max_rows: int = 50,
               min_identity: float = 0.3) -> str:
    """
    activity_map  : {compound_name: value}  — optional; adds a sortable activity column.
    act_label     : column header label (e.g. 'IC50 (nM)', 'pIC50').
    max_rows      : maximum number of compounds shown in table (default 50).
    min_identity  : exclude compounds with identity < this fraction (default 0.3).
    """
    ref_n      = len(ref_obj.positions())
    ref_syms   = _main_chain_syms(ref_obj)
    ref_labels = _main_chain_labels(ref_obj)
    valid_rows = [r for r in rows if r.get('obj') is not None]

    all_mw = [r['mw'] for r in valid_rows if r['mw'] is not None]
    mw_min = min(all_mw) if all_mw else 1000
    mw_max = max(all_mw) if all_mw else 2000

    # Quick identity screen: align all to reference, filter, sort, cap
    if valid_rows:
        _, mr_full, _, als_full, _ = _align_group_rows(ref_syms, valid_rows)
        def _quick_id(syms_al):
            return sum(
                1 for c, ms in enumerate(mr_full)
                if ms is not None and c < len(syms_al) and syms_al[c] == ms
            ) / ref_n if ref_n else 0.0

        ref_rows   = [r for r, s in zip(valid_rows, als_full) if r.get('is_ref')]
        other      = [(r, _quick_id(s)) for r, s in zip(valid_rows, als_full)
                      if not r.get('is_ref')]
        passing    = sorted([(r, i) for r, i in other if i >= min_identity],
                            key=lambda x: -x[1])
        excluded   = [(r, i) for r, i in other if i < min_identity]
        cap = max(0, max_rows - len(ref_rows))
        if len(passing) > cap:
            excluded.extend(passing[cap:])
            passing = passing[:cap]
        display_rows = ref_rows + [r for r, _ in passing]
        excluded_info = sorted(excluded, key=lambda x: -x[1])
    else:
        display_rows = []
        excluded_info = []

    _GAP_CELL = '<td class="gap-cell" data-sort=""></td>'

    def _cell(sym, hov_cls, color_fn, desc=None):
        bg = color_fn(sym, desc, db)
        fg = text_color(bg)
        safe = html.escape(sym)
        return (
            f'<td class="cell {hov_cls}" '
            f'style="background:{bg};color:{fg}" '
            f'data-tip="{safe}" data-sort="{safe}" '
            f'onmouseenter="highlight(\'{hov_cls}\')" onmouseleave="unhighlight(\'{hov_cls}\')">'
            f'<span>{safe}</span></td>'
        )

    def _score_td(score, is_ref):
        if is_ref:
            return '<td class="score-col" style="font-style:italic;color:#336" data-sort="1.0">REF</td>'
        if score is None:
            return '<td class="score-col" data-sort="-1">—</td>'
        pct = int(score * 100)
        bar_w = int(score * 50)
        clr = '#22aa55' if score >= 0.85 else ('#aaaa00' if score >= 0.70 else '#cc4422')
        return (
            f'<td class="score-col" data-sort="{score:.4f}">'
            f'<div class="score-bar-wrap">'
            f'<div class="score-bar" style="width:{bar_w}px;background:{clr}"></div>'
            f'<span style="font-size:0.65rem">{pct}%</span>'
            f'</div></td>'
        )

    def _identity_td(identity, is_ref, grp_ref_n):
        if is_ref:
            return '<td class="id-col" style="font-style:italic;color:#336" data-sort="1.0">REF</td>'
        if identity is None:
            return '<td class="id-col" style="color:#aaa" data-sort="-1">—</td>'
        n_same = int(round(identity * grp_ref_n))
        bar_w  = int(identity * 50)
        clr = '#1a8a1a' if identity >= 0.93 else ('#cc8800' if identity >= 0.79 else '#cc2222')
        return (
            f'<td class="id-col" data-sort="{identity:.4f}">'
            f'<div class="score-bar-wrap">'
            f'<div class="score-bar" style="width:{bar_w}px;background:{clr}"></div>'
            f'<span style="font-size:0.65rem">{n_same}/{grp_ref_n}</span>'
            f'</div></td>'
        )

    def _mw_td(mw_val):
        if mw_val is None:
            return '<td class="mw-col" style="color:#999" data-sort="-1">N/A</td>'
        bg = _mw_color(mw_val, mw_min, mw_max)
        fg = text_color(bg)
        return (
            f'<td class="mw-col" style="background:{bg};color:{fg}" '
            f'data-sort="{mw_val:.1f}">{mw_val:.1f}</td>'
        )

    # Pre-compute activity range for colour-coding (lower = more potent if IC50-like)
    _act_vals = [v for v in (activity_map or {}).values()
                 if v is not None and isinstance(v, (int, float))]
    _act_min  = min(_act_vals) if _act_vals else 0.0
    _act_max  = max(_act_vals) if _act_vals else 1.0

    def _act_td(cid, is_ref):
        if activity_map is None:
            return ''
        if is_ref:
            return f'<td class="act-col" style="font-style:italic;color:#336" data-sort="-999">REF</td>'
        val = activity_map.get(cid)
        if val is None:
            return '<td class="act-col" style="color:#aaa" data-sort="-1">—</td>'
        try:
            fval = float(val)
            span = max(_act_max - _act_min, 1e-6)
            # potency: low value = more potent → green; high = red
            t = max(0.0, min(1.0, (fval - _act_min) / span))
            r = int(60 + t * 195); g = int(200 - t * 160)
            bg = f'#{r:02x}{g:02x}40'
            return (f'<td class="act-col" data-sort="{fval}" '
                    f'style="background:{bg}">{val}</td>')
        except (TypeError, ValueError):
            return f'<td class="act-col" data-sort="-1">{html.escape(str(val))}</td>'

    _charge_legend_items = [
        ('#e62020', 'Positive at pH 7',     'K, R, H, Dab, Orn'),
        ('#1e50e6', 'Negative at pH 7',     'D, E, gGlu'),
        ('#f8f8f8', 'Polar / Neutral',      'S, T, N, Q, G, P, C'),
        ('#ffdc00', 'Aliphatic (low LogP)', 'A, V, L, I — yellow'),
        ('#ff5500', 'Aliphatic (high LogP)','meI, meL, Pip — orange'),
        ('#ffaa00', 'Aromatic',              'F, W, Y, 1Nal — scaled'),
    ]

    def _group_section(title, subtitle, grp_rows, grp_ref_syms, grp_ref_labels, pfx):
        if not grp_rows:
            return ''
        grp_ref_n = len(grp_ref_syms)
        n_cols, master_ref, master_labels, al_syms, al_desc = _align_group_rows(
            grp_ref_syms, grp_rows, grp_ref_labels)

        # Per-row identity and changes vs group reference
        identities, changes = [], []
        for r, syms_al in zip(grp_rows, al_syms):
            is_ref = r.get('is_ref', False)
            if is_ref:
                identities.append(1.0); changes.append('—'); continue
            ref_p, matches, chg_list = 0, 0, []
            for c, ms in enumerate(master_ref):
                if ms is None:
                    continue
                qs = syms_al[c] if c < len(syms_al) else None
                if qs is None:
                    chg_list.append(f'P{ref_p+1}:{ms}→gap')
                elif qs != ms:
                    chg_list.append(f'P{ref_p+1}:{ms}→{qs}')
                else:
                    matches += 1
                ref_p += 1
            identities.append(matches / grp_ref_n if grp_ref_n else 0.0)
            changes.append(', '.join(chg_list) if chg_list else '—')

        # Consensus and conservation (only over ref positions in master_ref)
        ref_col_indices = [c for c, s in enumerate(master_ref) if s is not None]
        ref_col_labels  = [master_labels[c] for c in ref_col_indices]
        cons_syms, conservation = [], []
        for c in ref_col_indices:
            col_syms = [al_syms[i][c] for i in range(len(grp_rows))
                        if c < len(al_syms[i]) and al_syms[i][c] is not None]
            if not col_syms:
                cons_syms.append('?'); conservation.append(0.0); continue
            ctr = Counter(col_syms)
            top, cnt = ctr.most_common(1)[0]
            cons_syms.append(top)
            conservation.append(cnt / len(col_syms))

        _POS_OFF = 2
        _ANN_OFF = _POS_OFF + n_cols + 1  # +1 for right spacer column

        # Position headers:
        #   ms not None  → main chain position (blue header)
        #   ms None, ml has '.' → sidechain slot (green header)
        #   ms None, ml None → NW gap column (grey '·')
        pos_ths = ''
        for c, (ms, ml) in enumerate(zip(master_ref, master_labels)):
            hov_cls = f'hov-{pfx}-{c}'
            col_idx = _POS_OFF + c
            if ms is not None:
                # main chain position
                lbl = html.escape(str(ml)) if ml is not None else str(c)
                pos_ths += (
                    f'<th class="pos-hdr sortable" data-col="{col_idx}" '
                    f'style="background:#e8eef6;color:#1a3a6a" '
                    f'onclick="sortTable(\'{{TBL_ID}}\',{col_idx},\'str\')" '
                    f'onmouseenter="highlight(\'{hov_cls}\')" onmouseleave="unhighlight(\'{hov_cls}\')">'
                    f'{lbl}</th>'
                )
            elif ml is not None:
                # sidechain slot (ms=None, ml="N.K")
                lbl = html.escape(str(ml))
                pos_ths += (
                    f'<th class="pos-hdr sortable" data-col="{col_idx}" '
                    f'style="background:#d4edda;color:#155724" '
                    f'onclick="sortTable(\'{{TBL_ID}}\',{col_idx},\'str\')" '
                    f'onmouseenter="highlight(\'{hov_cls}\')" onmouseleave="unhighlight(\'{hov_cls}\')">'
                    f'{lbl}</th>'
                )
            else:
                # NW gap column
                pos_ths += (
                    f'<th class="pos-hdr" style="color:#aaa;background:#f0f0f8" '
                    f'onmouseenter="highlight(\'{hov_cls}\')" onmouseleave="unhighlight(\'{hov_cls}\')">'
                    f'&middot;</th>'
                )

        _act_th = (
            f'<th class="act-col sortable" data-col="{_ANN_OFF+4}" '
            f'onclick="sortTable(\'{{TBL_ID}}\',{_ANN_OFF+4},\'num\')">'
            f'{html.escape(act_label)}</th>'
        ) if activity_map is not None else ''

        ann_ths = (
            f'<th data-col="{_ANN_OFF-1}" style="min-width:4px;padding:0;background:#dee4f0"></th>'
            f'<th class="score-col sortable" data-col="{_ANN_OFF}" '
            f'onclick="sortTable(\'{{TBL_ID}}\',{_ANN_OFF},\'num\')">Align score</th>'
            f'<th class="id-col sortable" data-col="{_ANN_OFF+1}" '
            f'onclick="sortTable(\'{{TBL_ID}}\',{_ANN_OFF+1},\'num\')">% Identity</th>'
            f'<th class="mw-col sortable" data-col="{_ANN_OFF+2}" '
            f'onclick="sortTable(\'{{TBL_ID}}\',{_ANN_OFF+2},\'num\')">MW (Da)</th>'
            f'<th class="sortable" data-col="{_ANN_OFF+3}" '
            f'onclick="sortTable(\'{{TBL_ID}}\',{_ANN_OFF+3},\'str\')">Changes vs Ref</th>'
            f'{_act_th}'
        )
        fixed_ths = (
            f'<th class="row-hdr sortable" data-col="0" '
            f'onclick="sortTable(\'{{TBL_ID}}\',0,\'str\')">ID</th>'
            f'<th data-col="1" style="min-width:4px;padding:0;background:#dee4f0"></th>'
        )
        header_tmpl = fixed_ths + pos_ths + ann_ths

        def _build_rows(scheme):
            html_out = ''
            for r, syms_al, desc_al, identity, chg in zip(
                    grp_rows, al_syms, al_desc, identities, changes):
                is_ref = r.get('is_ref', False)
                hdr_cls = 'ref-hdr' if is_ref else ''
                row_cls = 'ref-row' if is_ref else ''
                cells = ''
                for c, (qs, qd) in enumerate(zip(syms_al, desc_al)):
                    hov_cls = f'hov-{pfx}-{c}'
                    if qs is None:
                        cells += _GAP_CELL
                    elif scheme == 'zappo':
                        cells += _cell(qs, hov_cls,
                                       lambda s, d, db_: residue_color(s, db_), desc={})
                    else:
                        cells += _cell(qs, hov_cls, _charge_logp_color, desc=qd)
                cid = r['id']
                html_out += (
                    f'<tr class="{row_cls}">'
                    f'<td class="row-hdr {hdr_cls}" data-sort="{html.escape(cid)}">'
                    f'{html.escape(cid)}</td>'
                    f'<td style="width:4px;padding:0;background:#dee4f0"></td>'
                    f'{cells}'
                    f'<td style="width:4px;padding:0;background:#dee4f0"></td>'
                    f'{_score_td(r.get("score"), is_ref)}'
                    f'{_identity_td(identity, is_ref, grp_ref_n)}'
                    f'{_mw_td(r.get("mw"))}'
                    f'<td class="chg-col">{html.escape(chg)}</td>'
                    f'{_act_td(r["id"], is_ref)}'
                    f'</tr>\n'
                )
            # Consensus row — ref positions filled, gap columns show gap cell
            ref_p = 0
            cons_cells = ''
            for c, ms in enumerate(master_ref):
                hov_cls = f'hov-{pfx}-{c}'
                if ms is None:
                    cons_cells += _GAP_CELL
                else:
                    cs = cons_syms[ref_p] if ref_p < len(cons_syms) else '?'
                    if scheme == 'zappo':
                        cons_cells += _cell(cs, hov_cls,
                                            lambda s, d, db_: residue_color(s, db_), desc={})
                    else:
                        cons_cells += _cell(cs, hov_cls, _charge_logp_color, desc={})
                    ref_p += 1
            html_out += (
                f'<tr class="sep-row">'
                f'<td class="row-hdr con-hdr" data-sort="~consensus">Consensus</td>'
                f'<td style="width:4px;padding:0;background:#dee4f0"></td>'
                f'{cons_cells}'
                f'<td style="width:4px;padding:0;background:#dee4f0"></td>'
                f'<td class="score-col"></td><td class="id-col"></td>'
                f'<td class="mw-col"></td><td class="chg-col"></td>'
                + ('<td class="act-col"></td>' if activity_map is not None else '')
                + '</tr>\n'
            )
            return html_out

        tid_z = f'tbl_{pfx}_z'
        tid_c = f'tbl_{pfx}_c'
        hdr_z = header_tmpl.replace('{TBL_ID}', tid_z)
        hdr_c = header_tmpl.replace('{TBL_ID}', tid_c)
        zappo_rows_html  = _build_rows('zappo')
        charge_rows_html = _build_rows('charge')

        cons_bars = ''.join(
            f'<div class="cons-col">'
            f'<div class="cons-bar-track">'
            f'<div class="cons-bar-fill" style="height:{int(conservation[i]*38)}px"></div>'
            f'</div>'
            f'<div class="cons-label">{int(conservation[i]*100)}%</div>'
            f'<div class="cons-pos">{html.escape(str(ref_col_labels[i]))}</div>'
            f'</div>'
            for i in range(len(ref_col_indices))
        )
        legend_zappo = ''.join(
            f'<div class="legend-item">'
            f'<div class="legend-swatch" style="background:{hex_}"></div>'
            f'<span><b>{label}</b> <span style="color:#777">({ex})</span></span>'
            f'</div>'
            for hex_, label, ex in LEGEND_CATEGORIES
        )
        legend_charge = ''.join(
            f'<div class="charge-legend-item">'
            f'<div class="charge-swatch" style="background:{clr};border-color:#999"></div>'
            f'<span><b>{lbl}</b> <span style="color:#777">({ex})</span></span>'
            f'</div>'
            for clr, lbl, ex in _charge_legend_items
        )

        return f"""
<h2>{html.escape(title)}</h2>
<p class="subtitle">{html.escape(subtitle)}</p>

<h3 style="font-size:0.85rem;margin:14px 0 6px;color:#2a3a5a">Zappo Chemistry Colours</h3>
<p class="subtitle">Residues coloured by chemistry class. Hover a position number or cell to highlight the column.</p>
<div class="aln-wrap">
<table class="aln" id="{tid_z}">
<thead><tr>{hdr_z}</tr></thead>
<tbody>{zappo_rows_html}</tbody>
</table>
</div>

<h3 style="font-size:0.85rem;margin:14px 0 4px;color:#2a3a5a">Conservation</h3>
<div class="cons-wrap">{cons_bars}</div>

<h3 style="font-size:0.85rem;margin:10px 0 4px;color:#2a3a5a">Colour Legend — Zappo</h3>
<div class="legend-wrap">{legend_zappo}</div>

<h3 style="font-size:0.85rem;margin:20px 0 6px;color:#2a3a5a">Charge / Lipophilicity Map</h3>
<p class="subtitle">positive (pH 7) → red &nbsp;|&nbsp; negative → blue &nbsp;|&nbsp; polar/neutral → white &nbsp;|&nbsp; aliphatic/aromatic → yellow–orange by LogP</p>
<div class="aln-wrap">
<table class="aln" id="{tid_c}">
<thead><tr>{hdr_c}</tr></thead>
<tbody>{charge_rows_html}</tbody>
</table>
</div>

<h3 style="font-size:0.85rem;margin:10px 0 4px;color:#2a3a5a">Colour Legend — Charge / LogP</h3>
<div class="charge-legend">{legend_charge}</div>
"""

    act_note = f' &nbsp;|&nbsp; {html.escape(act_label)} shown' if activity_map is not None else ''

    scatter_html = _make_scatter_plots(valid_rows, activity_map, act_label)

    sec = _group_section(
        'Sequence Alignment',
        (f'NW-aligned to reference. Gap columns (·) show insertions/deletions. '
         f'Main chain positions numbered by HELM index; sidechain sub-positions as N.1, N.2, … '
         f'(green headers, reading outward from backbone). '
         f'Showing {len(display_rows)} of {len(valid_rows)} compounds '
         f'(identity ≥ {min_identity:.0%}, max {max_rows} rows).'),
        display_rows, ref_syms, ref_labels, 'main',
    )

    # Excluded compounds note
    excl_html = ''
    if excluded_info:
        items = ', '.join(
            f'{html.escape(r["id"])} ({i:.0%})'
            for r, i in excluded_info[:15]
        )
        more = f' … and {len(excluded_info)-15} more' if len(excluded_info) > 15 else ''
        excl_html = (
            f'<p class="subtitle" style="color:#999;margin-top:8px">'
            f'Excluded {len(excluded_info)} compound(s) '
            f'(identity &lt; {min_identity:.0%} or display limit reached): '
            f'{items}{more}</p>'
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Peptide Library SAR Report</title>
<style>{_CSS}</style>
</head>
<body>
<h1>Peptide Library — Sequence Alignment &amp; SAR Report</h1>
<p class="subtitle">{len(valid_rows)} compounds &nbsp;|&nbsp; {ref_n}-mer reference &nbsp;|&nbsp;
MW = Σ residue MW − bonds × 18.015 &nbsp;|&nbsp; Click any column header to sort{act_note}</p>

{scatter_html}
{sec}
{excl_html}

<script>{_JS}</script>
</body>
</html>
"""
