"""Project NW-aligned rows onto a shared master column layout for the report table.

NW operates on main-chain symbols only; sidechain monomers are inserted after their
anchor main-chain position from HELM connectivity, never repositioned by NW.
"""
from __future__ import annotations

from scripts.helm_alignment import nw_align
from sar_report.flatten import main_chain_syms, token_descs


def align_group_rows(grp_ref_syms: list, rows: list,
                     grp_ref_labels: list | None = None) -> tuple:
    """
    Align rows to the reference backbone, then insert sidechain columns structurally.

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
        q = main_chain_syms(r['obj']) if r['obj'] is not None else []
        pairwise.append(nw_align(grp_ref_syms, q) if q != grp_ref_syms
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
        main_descs, sc_desc_map = token_descs(r['obj']) if r['obj'] is not None else ([], {})

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
