#!/usr/bin/env python3
"""
Recreate Fig 2a and Fig 2b from Østergaard et al. 2021 (Sci Rep 11:21179)
as dot/scatter plots with categorical x-axis.

Fig 2a: half-life vs lipidation position in PYY3-36 (C18 diacid-gGlu-2xAdo)
Fig 2b: half-life vs fatty acid length ± γGlu, fixed at position 30
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 150,
})

# ── DATA ─────────────────────────────────────────────────────────────────────────

# Fig 2a: position scan, all C18 diacid-gGlu-2xAdo (Table 1, analogues 1-33)
# (pyy_pos, half_life_h)
position_scan = [
    ('Nα',  17),
    ('4',   14),
    ('5',   8.8),
    ('6',   7.2),
    ('7',   11),
    ('8',   5.9),
    ('9',   12),
    ('10',  8.4),
    ('11',  11),
    ('12',  15),
    ('13',  13),
    ('14',  11),
    ('15',  8.8),
    ('16',  17),
    ('17',  39),
    ('18',  22),
    ('19',  29),
    ('20',  33),
    ('21',  34),
    ('22',  36),
    ('23',  19),
    ('24',  66),
    ('25',  55),
    ('26',  41),
    ('27',  52),
    ('28',  62),
    ('29',  39),
    ('30',  76),   # best position — highlighted green
    ('31',  75),
    ('32',  49),
    ('33',  56),
    ('34',  30),
    ('35',  67),
]

# Fig 2b: fatty acid length ± γGlu at position 30 (Table 1, analogues 34-43)
# Groups: without γGlu (analogues 37-40), with γGlu (analogues 34-36, 41-43 Ado variants)
# Format: (label, half_life_h, fatty_acid, has_gGlu)
fa_scan = [
    # without γGlu
    ('C14\n−γGlu',   2,  'C14', False),
    ('C16\n−γGlu',   4,  'C16', False),
    ('C18\n−γGlu',  13,  'C18', False),
    ('C20\n−γGlu',  20,  'C20', False),
    # with γGlu (2xAdo unless noted)
    ('C14\n+γGlu',   4,  'C14', True),
    ('C16\n+γGlu',  28,  'C16', True),
    ('C18\n+γGlu',  76,  'C18', True),   # analogue 28, 2xAdo
    ('C20\n+γGlu',  99,  'C20', True),
    # Ado spacer variants at C18-γGlu (all similar half-life)
    ('C18\nno Ado',  97,  'C18', True),  # analogue 41
    ('C18\n4xAdo',   75,  'C18', True),  # analogue 42
    ('C18\n6xAdo',   78,  'C18', True),  # analogue 43
]

# Color map for fatty acid chain length
FA_COLORS = {
    'C14': '#4472C4',  # blue
    'C16': '#70AD47',  # green
    'C18': '#FFC000',  # amber
    'C20': '#FF4B4B',  # red
}

# ── FIG 2a ───────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

ax = axes[0]

positions = [r[0] for r in position_scan]
halflives = [r[1] for r in position_scan]
x = np.arange(len(positions))

# Color: highlight pos 30 in green, others by half-life (grey-blue gradient)
colors_a = []
for p, hl in position_scan:
    if p == '30':
        colors_a.append('#2CA02C')   # paper green
    elif hl >= 50:
        colors_a.append('#1F77B4')   # blue for long half-life
    else:
        colors_a.append('#6B6B6B')   # grey for short

# Draw vertical lines from 0 to dot (lollipop style)
for xi, hl, col in zip(x, halflives, colors_a):
    ax.vlines(xi, 0, hl, colors=col, linewidth=1.5, alpha=0.7)

# Draw dots
scatter = ax.scatter(x, halflives, c=colors_a, s=70, zorder=5,
                     edgecolors='white', linewidths=0.5)

# Annotate the best position (30)
best_idx = positions.index('30')
ax.annotate('Pos 30\n(76 h)', xy=(best_idx, 76), xytext=(best_idx + 2, 82),
            fontsize=8, color='#2CA02C', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#2CA02C', lw=1.2))

# Horizontal reference line at 46 h (semaglutide half-life in minipig, cited in paper)
ax.axhline(46, linestyle='--', color='#FF7F0E', linewidth=1.2, alpha=0.8)
ax.text(len(x) - 0.5, 47, 'semaglutide\n(~46 h)', fontsize=7.5,
        color='#FF7F0E', ha='right', va='bottom')

ax.set_xticks(x)
ax.set_xticklabels(positions, rotation=45, ha='right')
ax.set_xlabel('Position in PYY₃₋₃₆')
ax.set_ylabel('Half-life in minipigs (hours)')
ax.set_title('(a)  Half-life vs. lipidation position\nC18 diacid-γGlu-2xAdo')
ax.set_ylim(0, 100)
ax.yaxis.grid(True, alpha=0.3, linestyle=':')
ax.set_axisbelow(True)

# Legend
long_patch = mpatches.Patch(color='#1F77B4', label='≥50 h')
short_patch = mpatches.Patch(color='#6B6B6B', label='<50 h')
best_patch  = mpatches.Patch(color='#2CA02C', label='pos 30 (best)')
ax.legend(handles=[best_patch, long_patch, short_patch],
          loc='upper left', fontsize=8, framealpha=0.8)

# ── FIG 2b ───────────────────────────────────────────────────────────────────────

ax2 = axes[1]

labels  = [r[0] for r in fa_scan]
hls     = [r[1] for r in fa_scan]
fas     = [r[2] for r in fa_scan]
gglus   = [r[3] for r in fa_scan]
x2      = np.arange(len(labels))

colors_b = [FA_COLORS[fa] for fa in fas]

# Background shading: left half = no γGlu, right half = +γGlu
n_no = sum(1 for g in gglus if not g)
ax2.axvspan(-0.5, n_no - 0.5, alpha=0.05, color='grey', label='')
ax2.axvspan(n_no - 0.5, len(labels) - 0.5, alpha=0.05, color='#2CA02C')

# Lollipop
for xi, hl, col in zip(x2, hls, colors_b):
    ax2.vlines(xi, 0, hl, colors=col, linewidth=1.8, alpha=0.8)
ax2.scatter(x2, hls, c=colors_b, s=90, zorder=5,
            edgecolors='white', linewidths=0.5)

# Value labels on dots
for xi, hl in zip(x2, hls):
    ax2.text(xi, hl + 2.5, str(hl), ha='center', va='bottom', fontsize=7.5, color='#333333')

# Group divider
ax2.axvline(n_no - 0.5, color='#999999', linestyle=':', linewidth=1)
ax2.text((n_no / 2) - 0.5, 128, '− γGlu', ha='center', fontsize=9,
         color='#555555', style='italic')
ax2.text(n_no + (len(labels) - n_no) / 2 - 0.5, 128, '+ γGlu', ha='center',
         fontsize=9, color='#2CA02C', style='italic', fontweight='bold')

ax2.set_xticks(x2)
ax2.set_xticklabels(labels, rotation=0, ha='center', fontsize=8.5)
ax2.set_xlabel('Fatty acid / linker (position 30)')
ax2.set_ylabel('Half-life in minipigs (hours)')
ax2.set_title('(b)  Half-life vs. fatty acid length ± γGlu\nall analogues at position 30')
ax2.set_ylim(0, 135)
ax2.yaxis.grid(True, alpha=0.3, linestyle=':')
ax2.set_axisbelow(True)

# FA color legend
legend_patches = [mpatches.Patch(color=FA_COLORS[k], label=k)
                  for k in ['C14', 'C16', 'C18', 'C20']]
ax2.legend(handles=legend_patches, title='Fatty acid', loc='upper left',
           fontsize=8, title_fontsize=8, framealpha=0.8)

# ── SAVE ─────────────────────────────────────────────────────────────────────────

_data_dir = Path(__file__).resolve().parent.parent / 'data'
out_png = str(_data_dir / 'pyy_fig2_halflives.png')
out_pdf = str(_data_dir / 'pyy_fig2_halflives.pdf')

fig.suptitle('Østergaard et al. 2021 — PYY₃₋₃₆ lipidation SAR (recreated)',
             fontsize=11, y=1.01, color='#444444')

fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
print(f'Saved:\n  {out_png}\n  {out_pdf}')
plt.show()
