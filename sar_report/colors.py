"""Colour schemes for the report: charge/LogP residue colouring and the MW badge gradient."""
from __future__ import annotations

from scripts.residue_colors import category_of

# Charge/LogP legend rows: (hex, label, examples)
CHARGE_LEGEND_ITEMS = [
    ('#e62020', 'Positive at pH 7',      'K, R, H, Dab, Orn'),
    ('#1e50e6', 'Negative at pH 7',      'D, E, gGlu'),
    ('#f8f8f8', 'Polar / Neutral',       'S, T, N, Q, G, P, C'),
    ('#ffdc00', 'Aliphatic (low LogP)',  'A, V, L, I — yellow'),
    ('#ff5500', 'Aliphatic (high LogP)', 'meI, meL, Pip — orange'),
    ('#ffaa00', 'Aromatic',              'F, W, Y, 1Nal — scaled'),
]


def charge_logp_color(sym: str, desc: dict, db=None) -> str:
    """Structural (not electrochemical) charge classification:
      'negative' = residue has >1 COOH (sidechain acid + backbone)
      'positive' = residue has extra amine(s), guanidinium, or imidazole
    His is marked positive even though its imidazole pKa≈6 means ~50% protonated at pH 7.
    True pH-7 charge would need RDKit's MoleculeEnumerator with Epik or a
    Henderson-Hasselbalch step — keeping the structural heuristic as-is.
    """
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


def mw_color(val: float | None, vmin: float, vmax: float) -> str:
    """White→yellow→orange gradient for the MW column badge."""
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
