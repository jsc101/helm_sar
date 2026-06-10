"""HTML report generation — the sequence-alignment SAR table with two colour schemes,
conservation bars, scatter plots, and sortable annotation columns."""
from __future__ import annotations

import html
from collections import Counter

from scripts.residue_colors import residue_color, text_color, LEGEND_CATEGORIES
from sar_report.assets import CSS, JS
from sar_report.colors import charge_logp_color, mw_color, CHARGE_LEGEND_ITEMS
from sar_report.columns import align_group_rows
from sar_report.flatten import main_chain_syms, main_chain_labels
from sar_report.plots import make_scatter_plots


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
    ref_syms   = main_chain_syms(ref_obj)
    ref_labels = main_chain_labels(ref_obj)
    valid_rows = [r for r in rows if r.get('obj') is not None]

    all_mw = [r['mw'] for r in valid_rows if r['mw'] is not None]
    mw_min = min(all_mw) if all_mw else 1000
    mw_max = max(all_mw) if all_mw else 2000

    # Quick identity screen: align all to reference, filter, sort, cap
    if valid_rows:
        _, mr_full, _, als_full, _ = align_group_rows(ref_syms, valid_rows)
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
        bg = mw_color(mw_val, mw_min, mw_max)
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

    def _group_section(title, subtitle, grp_rows, grp_ref_syms, grp_ref_labels, pfx):
        if not grp_rows:
            return ''
        grp_ref_n = len(grp_ref_syms)
        n_cols, master_ref, master_labels, al_syms, al_desc = align_group_rows(
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
                        cells += _cell(qs, hov_cls, charge_logp_color, desc=qd)
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
                        cons_cells += _cell(cs, hov_cls, charge_logp_color, desc={})
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
            for clr, lbl, ex in CHARGE_LEGEND_ITEMS
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
<p class="subtitle">positive (pH 7) → red &nbsp;|&nbsp; negative → blue &nbsp;|&nbsp; polar/neutral → white &nbsp;|&nbsp; aliphatic/aromatic → yellow–orange by LogP</p>
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

    scatter_html = make_scatter_plots(valid_rows, activity_map, act_label)

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
<style>{CSS}</style>
</head>
<body>
<h1>Peptide Library — Sequence Alignment &amp; SAR Report</h1>
<p class="subtitle">{len(valid_rows)} compounds &nbsp;|&nbsp; {ref_n}-mer reference &nbsp;|&nbsp;
MW = Σ residue MW − bonds × 18.015 &nbsp;|&nbsp; Click any column header to sort{act_note}</p>

{scatter_html}
{sec}
{excl_html}

<script>{JS}</script>
</body>
</html>
"""
