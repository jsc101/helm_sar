"""Characterization tests for the SAR report engine (scripts/ + report.py).

These lock in current behavior of the previously-untested half of the repo so
the report.py decomposition can be refactored against a safety net. They assert
structural invariants (parse, lipidation detection, MW, alignment, HTML output)
rather than exact pixel/markup, so they survive cosmetic changes but catch
real regressions.
"""
import logging

logging.disable(logging.WARNING)

import pytest

from monomer_db.monomer_db import MonomerDB
from scripts.helm_parser import HELMParser
from sar_report import build_data, build_html, calc_mw

REF = (
    "PEPTIDE1{I.K.P.E.A.P.G.E.D.A.S.P.E.E.L.N.R.Y.Y.A.S.L.R.H.Y.L.N.L.V.T.R.Q.R.Y}"
    "$$$$V2.0"
)
POS30 = (
    "PEPTIDE1{I.K.P.E.A.P.G.E.D.A.S.P.E.E.L.N.R.Y.Y.A.S.L.R.H.Y.L.N.K.V.T.R.Q.R.Y}"
    "|PEPTIDE2{[gGlu].[Ado].[Ado].[C18d]}$PEPTIDE2,PEPTIDE1,1:R1-30:R3$$V2.0"
)


@pytest.fixture(scope="module")
def db():
    return MonomerDB(extra_sources=["monomer_db/custom_monomers.json"])


# ── HELMParser + HELMObject queries ──────────────────────────────────────────

def test_ref_parses_to_34_residues():
    obj = HELMParser.parse(REF)
    assert len(obj.get_chain()["monomers"]) == 34


def test_ref_jpv_linear_no_sidechain():
    obj = HELMParser.parse(REF)
    assert obj.get_jpv().startswith("I-K-P-E-A-P-G-E-D-A")
    assert obj.get_lipidation_pos() is None
    assert obj.get_sidechain_string() == ""


def test_lipidated_analog_detects_position_and_sidechain():
    obj = HELMParser.parse(POS30)
    assert obj.get_lipidation_pos() == 30
    assert obj.get_sidechain_string() == "gGlu-Ado-Ado-C18d"


def test_main_chain_flat_token_count():
    obj = HELMParser.parse(REF)
    main = [t for t in obj.get_jpv_flat() if t["is_main"]]
    assert len(main) == 34


# ── MW calculation ───────────────────────────────────────────────────────────

def test_ref_mw_within_tolerance():
    mw = calc_mw(HELMParser.parse(REF))
    assert mw == pytest.approx(4050.455, abs=0.01)


def test_lipidated_mw_within_tolerance():
    mw = calc_mw(HELMParser.parse(POS30))
    assert mw == pytest.approx(4047.46, abs=0.01)


# ── build_data pipeline ──────────────────────────────────────────────────────

def test_build_data_row_schema(db):
    ref_obj, rows = build_data([("PYY-ref", REF), ("A1", POS30)], db)
    assert ref_obj is not None
    assert len(rows) == 2
    assert set(rows[0].keys()) >= {"id", "helm", "obj", "mw", "score", "is_ref"}


def test_build_data_reference_autodetected(db):
    _, rows = build_data([("PYY-ref", REF), ("A1", POS30)], db)
    is_ref = {r["id"]: r["is_ref"] for r in rows}
    assert is_ref["PYY-ref"] is True
    assert is_ref["A1"] is False


# ── build_html output ────────────────────────────────────────────────────────

def test_build_html_produces_valid_document(db):
    ref_obj, rows = build_data([("PYY-ref", REF), ("A1", POS30)], db)
    html = build_html(
        ref_obj, rows, db=db,
        activity_map={"A1": 8.0}, act_label="pEC50",
    )
    assert "<html" in html.lower()
    assert "<table" in html
    assert len(html) > 10000
