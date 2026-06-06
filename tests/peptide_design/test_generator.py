import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import Compound, Sidechain
from peptide_design.generator import generate_helm

def _simple_compound(backbone, site, bonds, monomers):
    sc = Sidechain(site=site, bonds=bonds, monomers=monomers)
    return Compound(name="test", main_chain=list(backbone), sidechains=[sc])

def test_no_sidechain():
    c = Compound(name="ref", main_chain=["I","K","P","E"])
    helm = generate_helm(c)
    assert helm == "PEPTIDE1{I.K.P.E}$$$$V2.0"

def test_collapse_c18_at_pos4():
    c = _simple_compound(
        ["I","K","P","K"],
        site=4,
        bonds=["R3-R1","R2-R1","R2-R1","R2-R1"],
        monomers=["gGlu","Ado","Ado","C18d"],
    )
    helm = generate_helm(c)
    assert "PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}" in helm
    assert "PEPTIDE2,PEPTIDE1,1:R1-4:R3" in helm
    chain_part = helm.split("$")[0]
    assert chain_part.count("PEPTIDE") == 2

def test_collapse_n_terminal_r1():
    c = _simple_compound(
        ["I","K","P","E"],
        site=1,
        bonds=["R1-R1","R2-R1","R2-R1","R2-R1"],
        monomers=["gGlu","Ado","Ado","C18d"],
    )
    helm = generate_helm(c)
    assert "PEPTIDE2,PEPTIDE1,1:R1-1:R1" in helm

def test_collapse_single_monomer():
    c = _simple_compound(["A","K"], site=2,
                         bonds=["R3-R1"], monomers=["C18d"])
    helm = generate_helm(c)
    assert "PEPTIDE2{[C18d]}" in helm
    assert "PEPTIDE2,PEPTIDE1,1:R1-2:R3" in helm

def test_fmt_multi_char_monomer_gets_brackets():
    c = _simple_compound(["A","K"], site=2,
                         bonds=["R3-R1","R2-R1"], monomers=["Ado","C18d"])
    helm = generate_helm(c)
    assert "[Ado]" in helm
    assert "[C18d]" in helm

def test_fmt_single_char_monomer_no_brackets():
    c = _simple_compound(["A","K","G"], site=3,
                         bonds=["R3-R1"], monomers=["G"])
    helm = generate_helm(c)
    assert "PEPTIDE2{G}" in helm

# ── Expand tests ─────────────────────────────────────────────────────────────

def test_expand_triggered_by_chem_monomer():
    """Triazole14 is CHEM — forces expand path."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    assert "CHEM1{[Triazole14]}" in helm

def test_expand_backbone_connection_uses_b1():
    """PEPTIDE2 (Hpg) connects to PEPTIDE1 at site 37 via b1=R3-R1."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    assert "PEPTIDE2,PEPTIDE1,1:R1-37:R3" in helm

def test_expand_triazole_connects_to_hpg():
    """CHEM1 (Triazole14) connects to PEPTIDE2 (Hpg): b2=R3-R1."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    assert "CHEM1,PEPTIDE2,1:R1-1:R3" in helm

def test_expand_aha_connects_to_triazole():
    """PEPTIDE3 (Aha) connects to CHEM1 (Triazole14): b3=R4-R3."""
    c = _simple_compound(
        ["A"] * 37,
        site=37,
        bonds=["R3-R1", "R3-R1", "R4-R3", "R2-R1", "R2-R1"],
        monomers=["Hpg", "Triazole14", "Aha", "A", "C18d"],
    )
    helm = generate_helm(c)
    assert "PEPTIDE3,CHEM1,1:R3-1:R4" in helm

def test_expand_non_default_bond_without_chem():
    """Non-default bond (R3-R1 on b2) forces expand even for all-PEPTIDE."""
    c = _simple_compound(
        ["A","K"],
        site=2,
        bonds=["R3-R1", "R3-R1"],
        monomers=["Hpg", "Aha"],
    )
    helm = generate_helm(c)
    assert "PEPTIDE2{[Hpg]}" in helm
    assert "PEPTIDE3{[Aha]}" in helm
    assert "PEPTIDE3,PEPTIDE2,1:R1-1:R3" in helm
