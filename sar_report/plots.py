"""Embedded scatter plots for the report (base64 PNG)."""
from __future__ import annotations


def make_scatter_plots(rows: list, activity_map: dict | None, act_label: str) -> str:
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
