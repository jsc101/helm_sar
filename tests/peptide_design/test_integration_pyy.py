import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import parse_row
from peptide_design.generator import generate_helm
from peptide_design.validator import validate_helm

_PYY_SEQ = list("IKPEAPGEDASPEE") + list("LNRYYASLRHYLNL") + list("VTRQRY")
assert len(_PYY_SEQ) == 34

def _base_row(name, backbone=None):
    row = {"Name": name}
    seq = backbone or _PYY_SEQ
    for i, aa in enumerate(seq, 1):
        row[str(i)] = aa
    return row

def _add_sc(row, site, bonds, monomers, suffix=""):
    row[f"Site{suffix}"] = str(site)
    for k, (b, m) in enumerate(zip(bonds, monomers), 1):
        row[f"b{k}{suffix}"] = b
        row[f"SC{k}{suffix}"] = m

def _run(row):
    c = parse_row(row)
    helm = generate_helm(c)
    ok, msg = validate_helm(helm)
    return helm, ok, msg

def test_pyy_ref_no_sidechain():
    row = _base_row("ref")
    helm, ok, msg = _run(row)
    assert ok, msg
    assert helm.startswith("PEPTIDE1{I.K.P.E.A.P.G.E.D.A.S.P.E.E.L.N.R.Y.Y.A.S.L.R.H.Y.L.N.L.V.T.R.Q.R.Y}")

def test_pyy_analog_pos30_c18():
    """Analog_28: C18 gGlu-2xAdo at HELM pos 30 (K30)."""
    backbone = list(_PYY_SEQ); backbone[29] = "K"
    row = _base_row("Analog_28", backbone)
    _add_sc(row, 30, ["R3-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}" in helm
    assert "PEPTIDE2,PEPTIDE1,1:R1-30:R3" in helm

def test_pyy_analog_n_terminal_r1():
    """Analog_01-posNa: sidechain attaches to N-terminus via R1."""
    row = _base_row("Analog_01-posNa")
    _add_sc(row, 1, ["R1-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "PEPTIDE2,PEPTIDE1,1:R1-1:R1" in helm

def test_pyy_analog_c14d_no_gglu():
    """Analog_37: 2xAdo-C14d (no gGlu) at pos 30."""
    backbone = list(_PYY_SEQ); backbone[29] = "K"
    row = _base_row("Analog_37", backbone)
    _add_sc(row, 30, ["R3-R1","R2-R1","R2-R1"], ["Ado","Ado","C14d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "PEPTIDE2{[Ado].[Ado].[C14d]}" in helm

def test_pyy_combinatorial_backbone_plus_lipid():
    """Analog_52: Arg2 + Gln16 + K30 + C18 lipid."""
    backbone = list(_PYY_SEQ)
    backbone[1]  = "R"
    backbone[15] = "Q"
    backbone[29] = "K"
    row = _base_row("Analog_52", backbone)
    _add_sc(row, 30, ["R3-R1","R2-R1","R2-R1","R2-R1"], ["gGlu","Ado","Ado","C18d"])
    helm, ok, msg = _run(row)
    assert ok, msg
    assert "R." in helm
    assert ".Q." in helm
