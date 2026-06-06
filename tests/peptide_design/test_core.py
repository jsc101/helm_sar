import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import Compound, Sidechain, parse_row

def _pyy_row():
    """Analog_28-pos30 from PYY library: C18 at HELM pos 30."""
    row = {"Name": "Analog_28-pos30"}
    seq = list("IKPEAPGEDASPEE") + list("LNRYYASLRHYLNL") + list("VTRQRY")
    assert len(seq) == 34
    for i, aa in enumerate(seq, 1):
        row[str(i)] = aa
    row["30"] = "K"
    row["Site"] = "30"
    row["b1"] = "R3-R1"
    row["SC1"] = "gGlu"
    row["b2"] = "R2-R1"
    row["SC2"] = "Ado"
    row["b3"] = "R2-R1"
    row["SC3"] = "Ado"
    row["b4"] = "R2-R1"
    row["SC4"] = "C18d"
    return row

def test_parse_row_name():
    c = parse_row(_pyy_row())
    assert c.name == "Analog_28-pos30"

def test_parse_row_main_chain_length():
    c = parse_row(_pyy_row())
    assert len(c.main_chain) == 34

def test_parse_row_main_chain_k_at_30():
    c = parse_row(_pyy_row())
    assert c.main_chain[29] == "K"  # 0-indexed

def test_parse_row_one_sidechain():
    c = parse_row(_pyy_row())
    assert len(c.sidechains) == 1

def test_parse_row_sidechain_site():
    c = parse_row(_pyy_row())
    assert c.sidechains[0].site == 30

def test_parse_row_sidechain_bonds():
    c = parse_row(_pyy_row())
    sc = c.sidechains[0]
    assert sc.bonds == ["R3-R1", "R2-R1", "R2-R1", "R2-R1"]

def test_parse_row_sidechain_monomers():
    c = parse_row(_pyy_row())
    sc = c.sidechains[0]
    assert sc.monomers == ["gGlu", "Ado", "Ado", "C18d"]

def test_parse_row_no_sidechain():
    row = {"Name": "ref"}
    for i, aa in enumerate(list("IKPEAPGE"), 1):
        row[str(i)] = aa
    c = parse_row(row)
    assert c.sidechains == []

def test_parse_row_two_sidechains():
    row = {"Name": "bis"}
    for i, aa in enumerate(list("IKPEAPGE"), 1):
        row[str(i)] = aa
    row.update({"Site": "3", "b1": "R3-R1", "SC1": "Ado", "b2": "R2-R1", "SC2": "C18d"})
    row.update({"Site_2": "6", "b1_2": "R3-R1", "SC1_2": "Ado", "b2_2": "R2-R1", "SC2_2": "C12d"})
    c = parse_row(row)
    assert len(c.sidechains) == 2
    assert c.sidechains[1].site == 6
    assert c.sidechains[1].monomers == ["Ado", "C12d"]

def test_parse_row_sidechain_gap_skipped():
    """Site absent but Site_2 present — should still find Site_2."""
    row = {"Name": "gap"}
    for i, aa in enumerate(list("IKPE"), 1):
        row[str(i)] = aa
    row.update({"Site_2": "3", "b1_2": "R3-R1", "SC1_2": "Ado"})
    c = parse_row(row)
    assert len(c.sidechains) == 1
    assert c.sidechains[0].site == 3
