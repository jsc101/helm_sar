"""Characterization tests for the unified tabular reader."""
from __future__ import annotations

import pytest

from scripts.table_io import read_table


def _write(path, text):
    path.write_text(text, encoding='utf-8')
    return str(path)


def test_read_csv_headers_and_rows(tmp_path):
    p = _write(tmp_path / 'in.csv',
               "Name,HELM,pEC50\n A1 , PEPTIDE1{A}$$$$V2.0 , 7.5\nA2,PEPTIDE1{C}$$$$V2.0,\n")
    headers, rows = read_table(p)
    assert headers == ['Name', 'HELM', 'pEC50']
    assert len(rows) == 2
    # values are stripped
    assert rows[0] == {'Name': 'A1', 'HELM': 'PEPTIDE1{A}$$$$V2.0', 'pEC50': '7.5'}
    # empty cell -> ''
    assert rows[1]['pEC50'] == ''


def test_read_tsv_uses_tab_delimiter(tmp_path):
    p = _write(tmp_path / 'in.tsv', "Name\tHELM\nA1\tPEPTIDE1{A}$$$$V2.0\n")
    headers, rows = read_table(p)
    assert headers == ['Name', 'HELM']
    assert rows[0]['HELM'] == 'PEPTIDE1{A}$$$$V2.0'


def test_read_csv_skips_utf8_bom(tmp_path):
    p = tmp_path / 'bom.csv'
    p.write_bytes('﻿Name,HELM\nA1,PEPTIDE1{A}$$$$V2.0\n'.encode('utf-8'))
    headers, _ = read_table(str(p))
    assert headers[0] == 'Name'


def test_read_xlsx_roundtrip(tmp_path):
    openpyxl = pytest.importorskip('openpyxl')
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(['Name', 'HELM', 'pEC50'])
    ws.append(['A1', 'PEPTIDE1{A}$$$$V2.0', 7.5])
    ws.append([None, None, None])  # blank row should be dropped
    p = tmp_path / 'in.xlsx'
    wb.save(p)
    headers, rows = read_table(str(p))
    assert headers == ['Name', 'HELM', 'pEC50']
    assert len(rows) == 1
    assert rows[0]['Name'] == 'A1'
    assert rows[0]['pEC50'] == '7.5'
