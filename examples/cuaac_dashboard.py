"""CuAAC triazole → HELM pipeline dashboard.

Runs a representative set of CuAAC-stapled peptides through smiles_to_helm()
and emits a self-contained HTML report showing:

  - Fragment breakdown (colour-coded PEPTIDE / CHEM residue boxes)
  - Syntax-highlighted HELM string
  - Connection topology
  - Coverage statistics (100% = all residues resolved)

Usage
-----
    python -m examples.cuaac_dashboard              # → data/cuaac_dashboard.html
    python -m examples.cuaac_dashboard --out /tmp/out.html
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

# Allow running from repo root without pip install
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from monomer_db.monomer_db import MonomerDB
from scripts.smiles_to_helm import fragment_smiles, smiles_to_helm

# ---------------------------------------------------------------------------
# Test compounds
# ---------------------------------------------------------------------------

COMPOUNDS = [
    {
        "name": "Aha–Hpg macrolactam",
        "note": "Two-residue head-to-tail macrolactam; minimal CuAAC staple",
        "smiles": "N[C@@H]1CCn2cc(nn2)CCC[C@H](C(=O)O)NC1=O",
        "cyclic": False,
    },
    {
        "name": "Aha–Pra macrolactam",
        "note": "Shorter alkyne handle (Pra, 1-carbon) vs Hpg (3-carbon)",
        "smiles": "N[C@@H]1CCn2cc(C[C@H](NC1=O)C(=O)O)nn2",
        "cyclic": False,
    },
    {
        "name": "Ala–Aha⌁Hpg–Ala linear staple",
        "note": "Four-residue linear peptide; Aha at pos 2, Hpg at pos 3",
        "smiles": (
            "C[C@@H](NC(=O)[C@@H](CCn1cc(CCC[C@@H](C(=O)N[C@@H](C)C(=O)O)"
            "N)nn1)N)C(=O)O"
        ),
        "cyclic": False,
    },
    {
        "name": "Gly–Aha⌁Hpg–Gly linear staple",
        "note": "Glycine flanks; tests achiral backbone with click staple",
        "smiles": (
            "NCC(=O)N[C@@H](CCn1cc(CCC[C@@H](C(=O)NCC(=O)O)N)nn1)C(=O)O"
        ),
        "cyclic": False,
    },
    {
        "name": "Aha–Ala–Hpg extended staple",
        "note": "Non-adjacent Aha/Hpg with intervening Ala; exposes ordering limitation",
        "smiles": (
            "N[C@@H](CCn1cc(CCC[C@@H](C(=O)N[C@@H](C)C(=O)O)N)nn1)"
            "C(=O)N[C@@H](C)C(=O)O"
        ),
        "cyclic": False,
    },
]

# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

def run_pipeline(compounds: list[dict]) -> list[dict]:
    """Return one result dict per compound."""
    db = MonomerDB()
    results = []

    for c in compounds:
        smiles = c["smiles"]
        helm, new_monomers = smiles_to_helm(smiles, db=db)
        frags = fragment_smiles(smiles, db=db)

        total = len(frags)
        known = sum(1 for f in frags if f.symbol and not f.symbol.startswith("UNK_"))
        coverage = known / total if total else 0.0

        chain_info = _build_chain_info(frags)

        results.append({
            **c,
            "helm":         helm,
            "fragments":    frags,
            "coverage":     coverage,
            "known":        known,
            "total":        total,
            "new_monomers": new_monomers,
            "chain_info":   chain_info,
        })

    return results


def _build_chain_info(frags) -> dict:
    """Return dicts describing PEPTIDE and CHEM chains for display."""
    peptide, chem = [], []
    for i, f in enumerate(frags):
        is_chem = bool(f.entry and f.entry.get("polymerType") == "CHEM")
        info = {
            "idx":     i,
            "symbol":  f.symbol or "???",
            "is_unk":  not f.symbol or f.symbol.startswith("UNK_"),
            "is_chem": is_chem,
            "smiles":  f.smiles,
            "cross_bonds": list(f.cross_bonds),
            "r3_partners": list(f.r3_partners),
        }
        if is_chem:
            chem.append(info)
        else:
            peptide.append(info)
    return {"peptide": peptide, "chem": chem}


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

_CSS = """
:root {
  --bg: #0f1117;
  --surface: #1a1d27;
  --surface2: #23263a;
  --border: #2e3248;
  --text: #e2e4f0;
  --muted: #7c84a0;
  --peptide: #3b6cf6;
  --chem: #c46ef5;
  --ok: #2eb87b;
  --warn: #e8a733;
  --err: #e85555;
  --code-bg: #0c0e15;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: 'SF Pro Text', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  background: var(--bg); color: var(--text);
  padding: 32px 24px; max-width: 1100px; margin: 0 auto; line-height: 1.5;
}
h1 { font-size: 1.65rem; font-weight: 700; margin-bottom: 4px; }
.subtitle { color: var(--muted); font-size: 0.9rem; margin-bottom: 32px; }

/* Summary bar */
.summary-bar {
  display: flex; gap: 20px; flex-wrap: wrap;
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 10px; padding: 16px 20px; margin-bottom: 36px;
}
.stat { display: flex; flex-direction: column; }
.stat-val { font-size: 1.7rem; font-weight: 700; line-height: 1; }
.stat-label { font-size: 0.75rem; color: var(--muted); margin-top: 3px; text-transform: uppercase; letter-spacing: .05em; }

/* Cards */
.card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 12px; margin-bottom: 28px; overflow: hidden;
}
.card-header {
  display: flex; align-items: center; justify-content: space-between;
  padding: 14px 20px; border-bottom: 1px solid var(--border);
  background: var(--surface2);
}
.card-title { font-size: 1.05rem; font-weight: 600; }
.card-note { font-size: 0.8rem; color: var(--muted); max-width: 55ch; }
.card-body { padding: 18px 20px; }

/* Coverage badge */
.cov-badge {
  display: inline-flex; align-items: center; gap: 6px;
  border-radius: 20px; padding: 3px 12px; font-size: 0.82rem; font-weight: 600;
  white-space: nowrap;
}
.cov-100 { background: rgba(46,184,123,.18); color: var(--ok); }
.cov-part { background: rgba(232,167,51,.18); color: var(--warn); }
.cov-fail { background: rgba(232,85,85,.18); color: var(--err); }

/* Chain strips */
.chain-strip { display: flex; flex-wrap: wrap; gap: 6px; align-items: center; margin-bottom: 14px; }
.chain-label {
  font-size: 0.7rem; font-weight: 700; text-transform: uppercase;
  letter-spacing: .08em; color: var(--muted); margin-right: 4px; min-width: 62px;
}
.residue {
  display: inline-flex; flex-direction: column; align-items: center;
  border-radius: 7px; padding: 5px 9px;
  font-size: 0.8rem; font-weight: 600; cursor: default;
  transition: transform .1s;
  border: 1.5px solid transparent;
}
.residue:hover { transform: translateY(-2px); }
.residue .sym { font-size: 0.85rem; }
.residue .type-tag { font-size: 0.6rem; font-weight: 400; margin-top: 2px; opacity: .75; }
.res-peptide { background: rgba(59,108,246,.18); border-color: rgba(59,108,246,.4); color: #8baeff; }
.res-chem    { background: rgba(196,110,245,.18); border-color: rgba(196,110,245,.4); color: #d9a0fc; }
.res-unk     { background: rgba(232,85,85,.15);  border-color: rgba(232,85,85,.4);  color: #ff8a8a; }

/* Connection lines */
.connections { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 14px; }
.conn-pill {
  display: inline-flex; align-items: center; gap: 4px;
  border-radius: 20px; padding: 3px 10px; font-size: 0.75rem; font-family: monospace;
  background: rgba(196,110,245,.1); border: 1px solid rgba(196,110,245,.3); color: #d9a0fc;
}
.conn-arrow { color: var(--muted); }

/* HELM code block */
.helm-block {
  background: var(--code-bg); border-radius: 8px; padding: 14px 16px;
  font-family: 'JetBrains Mono', 'Fira Code', 'Cascadia Code', monospace;
  font-size: 0.78rem; overflow-x: auto; white-space: pre; line-height: 1.55;
  border: 1px solid var(--border);
}
.hl-chain  { color: #61c7fa; }
.hl-mono   { color: #d0a0fc; }
.hl-chem   { color: #f0c040; }
.hl-conn   { color: #7ee8a2; }
.hl-sep    { color: #5c6380; }
.hl-ver    { color: #5c6380; }
.copy-btn {
  float: right; background: var(--surface2); color: var(--muted);
  border: 1px solid var(--border); border-radius: 6px; padding: 3px 10px;
  font-size: 0.72rem; cursor: pointer; margin-left: 10px;
}
.copy-btn:hover { color: var(--text); border-color: var(--peptide); }

/* Legend */
.legend {
  display: flex; gap: 18px; flex-wrap: wrap;
  font-size: 0.78rem; color: var(--muted); margin-bottom: 24px;
}
.legend-item { display: flex; align-items: center; gap: 6px; }
.legend-dot {
  width: 10px; height: 10px; border-radius: 3px; flex-shrink: 0;
}
.ld-peptide { background: rgba(59,108,246,.5); }
.ld-chem    { background: rgba(196,110,245,.5); }
.ld-unk     { background: rgba(232,85,85,.5); }

/* Limitation banner */
.limitation {
  font-size: 0.8rem; color: var(--warn);
  background: rgba(232,167,51,.08); border: 1px solid rgba(232,167,51,.25);
  border-radius: 7px; padding: 8px 14px; margin-bottom: 14px;
}
"""

_JS = """
function copyHELM(id) {
  const el = document.getElementById(id);
  if (!el) return;
  navigator.clipboard.writeText(el.dataset.raw || el.innerText)
    .then(() => {
      const btn = el.previousElementSibling;
      if (btn) { btn.innerText = 'Copied!'; setTimeout(() => btn.innerText = 'Copy', 1500); }
    });
}
"""


def _cov_class(cov: float) -> str:
    if cov >= 1.0:
        return "cov-100"
    if cov >= 0.5:
        return "cov-part"
    return "cov-fail"


def _res_class(info: dict) -> str:
    if info["is_unk"]:
        return "res-unk"
    if info["is_chem"]:
        return "res-chem"
    return "res-peptide"


def _res_label(info: dict) -> str:
    sym = info["symbol"]
    if info["is_chem"]:
        return "CHEM"
    if info["is_unk"]:
        return "UNK"
    return "AA"


def _highlight_helm(helm: str) -> str:
    """Return HTML-highlighted HELM string."""
    if not helm:
        return '<span style="color:var(--err)">None</span>'

    # Split at $ separators (HELM sections: chains / connections / groups / version)
    sep = '<span class="hl-sep">$</span>'
    parts = helm.split("$")
    highlighted_parts = []

    for k, part in enumerate(parts):
        if k == 0:
            # Chain declarations: PEPTIDE1{...}|CHEM1{...}
            chain_html = []
            for chain in part.split("|"):
                if "{" in chain:
                    name, rest = chain.split("{", 1)
                    inner = rest.rstrip("}")
                    # Highlight each monomer token
                    tokens = inner.split(".")
                    token_html = []
                    for t in tokens:
                        if t.startswith("[") and t.endswith("]"):
                            token_html.append(
                                f'<span class="hl-mono">{t}</span>'
                            )
                        else:
                            token_html.append(
                                f'<span class="hl-mono">{t}</span>'
                            )
                    inner_html = '<span class="hl-sep">.</span>'.join(token_html)
                    if name.startswith("CHEM"):
                        chain_html.append(
                            f'<span class="hl-chem">{name}</span>{{'
                            f'{inner_html}}}'
                        )
                    else:
                        chain_html.append(
                            f'<span class="hl-chain">{name}</span>{{'
                            f'{inner_html}}}'
                        )
                else:
                    chain_html.append(part)
            highlighted_parts.append(
                '<span class="hl-sep">|</span>'.join(chain_html)
            )
        elif k == 1:
            # Connection lines
            if part:
                conn_html = []
                for conn in part.split("|"):
                    conn_html.append(f'<span class="hl-conn">{conn}</span>')
                highlighted_parts.append('<span class="hl-sep">|</span>'.join(conn_html))
            else:
                highlighted_parts.append("")
        elif k == len(parts) - 1:
            # Version
            highlighted_parts.append(f'<span class="hl-ver">{part}</span>')
        else:
            highlighted_parts.append(part)

    return sep.join(highlighted_parts)


def _card_html(result: dict, idx: int) -> str:
    cov = result["coverage"]
    cov_pct = f"{cov * 100:.0f}%"
    cov_cls = _cov_class(cov)
    has_non_adjacent = "⌁" in result["name"]  # flag in name

    chain_info = result["chain_info"]
    helm = result["helm"] or ""

    # Build PEPTIDE strip
    peptide_residues = ""
    for info in chain_info["peptide"]:
        rc = _res_class(info)
        label = _res_label(info)
        sym = info["symbol"]
        title = f'SMILES: {info["smiles"]}'
        peptide_residues += (
            f'<div class="residue {rc}" title="{title}">'
            f'<span class="sym">{sym}</span>'
            f'<span class="type-tag">{label}</span>'
            f"</div>"
        )

    # Build CHEM strip
    chem_residues = ""
    for info in chain_info["chem"]:
        rc = _res_class(info)
        label = _res_label(info)
        sym = info["symbol"]
        title = f'SMILES: {info["smiles"]}'
        chem_residues += (
            f'<div class="residue {rc}" title="{title}">'
            f'<span class="sym">{sym}</span>'
            f'<span class="type-tag">{label}</span>'
            f"</div>"
        )

    # Build connection pills from cross_bonds
    conn_pills = ""
    seen: set[frozenset] = set()
    for info in chain_info["peptide"] + chain_info["chem"]:
        for cb in info["cross_bonds"]:
            j = cb.partner_idx
            key: frozenset = frozenset([(info["idx"], cb.my_rgroup), (j, cb.partner_rgroup)])
            if key in seen:
                continue
            seen.add(key)
            src_sym = result["fragments"][info["idx"]].symbol or "?"
            dst_sym = result["fragments"][j].symbol or "?"
            conn_pills += (
                f'<div class="conn-pill">'
                f'{src_sym}:{cb.my_rgroup}'
                f'<span class="conn-arrow">→</span>'
                f'{dst_sym}:{cb.partner_rgroup}'
                f"</div>"
            )
        for p in info["r3_partners"]:
            key = frozenset([info["idx"], p])
            if key in seen:
                continue
            seen.add(key)
            src_sym = result["fragments"][info["idx"]].symbol or "?"
            dst_sym = result["fragments"][p].symbol or "?"
            conn_pills += (
                f'<div class="conn-pill">'
                f'{src_sym}:R3'
                f'<span class="conn-arrow">⟷</span>'
                f'{dst_sym}:R3'
                f"</div>"
            )

    helm_id = f"helm-{idx}"
    limitation_html = ""
    if has_non_adjacent:
        limitation_html = (
            '<div class="limitation">'
            "⚠ Non-adjacent Aha/Hpg: backbone ordering in PEPTIDE1 may not reflect "
            "linear sequence — intervening residues shift to end of fragment list (known limitation)."
            "</div>"
        )

    unk_warn = ""
    if result["new_monomers"]:
        syms = ", ".join(m["symbol"] for m in result["new_monomers"])
        unk_warn = (
            f'<div class="limitation">'
            f"Unknown fragments assigned UNK symbols: {syms}"
            f"</div>"
        )

    return f"""
<div class="card">
  <div class="card-header">
    <div>
      <div class="card-title">{result['name']}</div>
      <div class="card-note">{result['note']}</div>
    </div>
    <div class="cov-badge {cov_cls}">{result['known']}/{result['total']} &nbsp;{cov_pct}</div>
  </div>
  <div class="card-body">
    {limitation_html}{unk_warn}
    <div class="chain-strip">
      <span class="chain-label">PEPTIDE1</span>
      {peptide_residues}
    </div>
    {"<div class='chain-strip'><span class='chain-label'>CHEM</span>" + chem_residues + "</div>" if chem_residues else ""}
    {"<div class='connections'>" + conn_pills + "</div>" if conn_pills else ""}
    <div style="position:relative">
      <button class="copy-btn" onclick="copyHELM('{helm_id}')">Copy</button>
      <div class="helm-block" id="{helm_id}" data-raw="{helm}">{_highlight_helm(helm)}</div>
    </div>
  </div>
</div>
"""


def build_html(results: list[dict]) -> str:
    total_compounds = len(results)
    total_100 = sum(1 for r in results if r["coverage"] >= 1.0)
    total_frags = sum(r["total"] for r in results)
    total_known = sum(r["known"] for r in results)
    overall_cov = total_known / total_frags if total_frags else 0

    cards_html = "".join(_card_html(r, i) for i, r in enumerate(results))

    return textwrap.dedent(f"""\
        <!DOCTYPE html>
        <html lang="en">
        <head>
          <meta charset="UTF-8">
          <meta name="viewport" content="width=device-width, initial-scale=1.0">
          <title>CuAAC → HELM Dashboard</title>
          <style>{_CSS}</style>
        </head>
        <body>
          <h1>CuAAC Triazole → HELM Pipeline</h1>
          <div class="subtitle">
            Click chemistry fragmentation dashboard &mdash;
            azide (Aha/Pra) + alkyne (Hpg/Pra) staples through the SMILES→HELM engine
          </div>

          <div class="summary-bar">
            <div class="stat">
              <span class="stat-val">{total_compounds}</span>
              <span class="stat-label">Compounds</span>
            </div>
            <div class="stat">
              <span class="stat-val">{total_100} / {total_compounds}</span>
              <span class="stat-label">100% coverage</span>
            </div>
            <div class="stat">
              <span class="stat-val">{overall_cov * 100:.0f}%</span>
              <span class="stat-label">Overall fragment coverage</span>
            </div>
            <div class="stat">
              <span class="stat-val">{total_frags}</span>
              <span class="stat-label">Total fragments</span>
            </div>
          </div>

          <div class="legend">
            <span class="legend-item"><span class="legend-dot ld-peptide"></span> PEPTIDE monomer</span>
            <span class="legend-item"><span class="legend-dot ld-chem"></span> CHEM monomer (triazole)</span>
            <span class="legend-item"><span class="legend-dot ld-unk"></span> Unknown (UNK_*)</span>
          </div>

          {cards_html}

          <script>{_JS}</script>
        </body>
        </html>
    """)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate CuAAC → HELM HTML dashboard")
    parser.add_argument("--out", default=None, help="Output HTML path (default: data/cuaac_dashboard.html)")
    args = parser.parse_args()

    out_path = Path(args.out) if args.out else _ROOT / "data" / "cuaac_dashboard.html"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Running pipeline…")
    results = run_pipeline(COMPOUNDS)
    for r in results:
        cov_pct = f"{r['coverage'] * 100:.0f}%"
        status = "OK" if r["coverage"] >= 1.0 else "partial"
        print(f"  {r['name']:<40} {r['known']}/{r['total']} frags  {cov_pct}  [{status}]")

    html = build_html(results)
    out_path.write_text(html, encoding="utf-8")
    print(f"\nDashboard written → {out_path}")


if __name__ == "__main__":
    main()
