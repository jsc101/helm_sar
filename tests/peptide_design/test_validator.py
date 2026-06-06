import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))

from peptide_design.core import Compound, Sidechain
from peptide_design.generator import generate_helm
from peptide_design.validator import validate_helm

def _make_helm(backbone, site=None, bonds=None, monomers=None):
    scs = []
    if site is not None:
        scs = [Sidechain(site=site, bonds=bonds, monomers=monomers)]
    c = Compound(name="t", main_chain=list(backbone), sidechains=scs)
    return generate_helm(c)

def test_valid_no_sidechain():
    helm = _make_helm(["A","K","P","E"])
    ok, msg = validate_helm(helm)
    assert ok, msg

def test_valid_with_collapse_sidechain():
    helm = _make_helm(["A","K"], site=2,
                      bonds=["R3-R1","R2-R1"], monomers=["Ado","C18d"])
    ok, msg = validate_helm(helm)
    assert ok, msg

def test_invalid_helm_string():
    ok, msg = validate_helm("NOT_VALID_HELM")
    assert not ok
    assert msg

def test_unknown_monomer_flagged():
    c = Compound(name="t", main_chain=["A","K"],
                 sidechains=[Sidechain(site=2, bonds=["R3-R1"], monomers=["C16m"])])
    helm = generate_helm(c)
    ok, msg = validate_helm(helm)
    assert not ok
    assert "C16m" in msg or "unknown" in msg.lower()
