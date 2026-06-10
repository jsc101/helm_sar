"""Static CSS and JS for the SAR report HTML."""

CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Courier New', monospace; background: #f4f6f9; color: #222; padding: 20px; }
h1 { font-size: 1.3rem; font-weight: 700; margin-bottom: 4px; color: #1a2a4a; }
.subtitle { font-size: 0.8rem; color: #555; margin-bottom: 22px; }
h2 { font-size: 1.0rem; font-weight: 700; margin: 22px 0 8px; color: #2a3a5a;
     border-left: 4px solid #4477cc; padding-left: 8px; }

/* ── Shared table shell ───────────────────────── */
.aln-wrap { overflow-x: auto; }
table.aln { border-collapse: collapse; font-size: 0.72rem; table-layout: auto; }
table.aln th, table.aln td { border: 1px solid #ccc; text-align: center; white-space: nowrap; }

/* fixed left columns */
table.aln th.row-hdr { text-align: right; padding: 2px 6px; font-weight: 600; min-width: 120px; }
table.aln td.row-hdr { text-align: right; padding: 2px 6px; min-width: 120px;
                        font-size: 0.68rem; color: #333; }
table.aln td.row-hdr.ref-hdr { font-weight: 700; color: #003388; }
table.aln td.row-hdr.con-hdr { font-weight: 700; color: #006600; font-style: italic; }
table.aln .score-col { min-width: 74px; padding: 2px 4px; }
table.aln .id-col    { min-width: 68px; padding: 2px 4px; }
table.aln .mw-col    { min-width: 76px; padding: 2px 6px; font-weight: 700; }
table.aln .chg-col   { min-width: 160px; max-width: 260px; padding: 2px 6px;
                        text-align: left; font-size: 0.65rem; color: #444; white-space: normal; }
table.aln .act-col   { min-width: 80px; padding: 2px 6px; font-weight: 600;
                        text-align: center; font-size: 0.7rem; }

/* residue cells */
.cell { width: 58px; height: 26px; padding: 0; position: relative; cursor: default; }
.gap-cell { width: 58px; height: 26px; background: #ebebeb; border-color: #ddd !important; }
.cell span { display: block; width: 100%; height: 100%; line-height: 26px;
             font-size: 0.7rem; font-weight: 600; overflow: hidden; text-overflow: ellipsis; }
.cell:hover::after {
  content: attr(data-tip); position: absolute; bottom: 110%; left: 50%;
  transform: translateX(-50%); background: #222; color: #fff;
  padding: 2px 6px; border-radius: 4px; font-size: 0.65rem; white-space: nowrap;
  z-index: 100; pointer-events: none;
}
tr.ref-row td { border-bottom: 2px solid #336; }
tr.sep-row td { border-top: 2px dashed #555; }

/* sortable header */
th.sortable { cursor: pointer; user-select: none; }
th.sortable:hover { background: #d0d8ee; }
th.sort-asc::after  { content: ' ▲'; font-size: 0.6rem; }
th.sort-desc::after { content: ' ▼'; font-size: 0.6rem; }
.pos-hdr { background: #e8eef6; font-weight: 700; font-size: 0.75rem;
           padding: 3px 2px; color: #1a3a6a; }

/* score bar */
.score-bar-wrap { display: flex; align-items: center; gap: 4px; height: 22px; }
.score-bar { height: 12px; border-radius: 3px; }

/* conservation */
.cons-wrap { display: flex; gap: 1px; margin: 6px 0 14px; overflow-x: auto; }
.cons-col { display: flex; flex-direction: column; align-items: center; width: 58px; flex-shrink: 0; }
.cons-bar-track { width: 46px; height: 40px; background: #e0e4ea; border-radius: 3px;
                  display: flex; align-items: flex-end; overflow: hidden; border: 1px solid #ccc; }
.cons-bar-fill { width: 100%; background: #2255bb; border-radius: 2px 2px 0 0; }
.cons-label { font-size: 0.62rem; color: #444; margin-top: 2px; }
.cons-pos { font-size: 0.6rem; color: #888; }

/* legend */
.legend-wrap { display: flex; flex-wrap: wrap; gap: 8px; margin: 6px 0 16px; }
.legend-item { display: flex; align-items: center; gap: 4px; font-size: 0.68rem; }
.legend-swatch { width: 16px; height: 16px; border-radius: 3px; border: 1px solid #aaa; flex-shrink: 0; }
.charge-legend { display: flex; flex-wrap: wrap; gap: 12px; margin: 6px 0 14px; font-size: 0.68rem; }
.charge-legend-item { display: flex; align-items: center; gap: 5px; }
.charge-swatch { width: 18px; height: 18px; border-radius: 3px; border: 1px solid #aaa; }
"""

JS = """
// ── Column hover highlighting ───────────────────────
function highlight(cls) {
  document.querySelectorAll('.' + cls).forEach(el => {
    el.style.outline = '2px solid #ff8800'; el.style.zIndex = '10';
  });
}
function unhighlight(cls) {
  document.querySelectorAll('.' + cls).forEach(el => {
    el.style.outline = ''; el.style.zIndex = '';
  });
}

// ── Table sorting ───────────────────────────────────
var _sortState = {};   // tableId → {col, dir}

function sortTable(tableId, colIdx, type) {
  var tbl   = document.getElementById(tableId);
  var tbody = tbl.tBodies[0];
  var rows  = Array.from(tbody.rows);

  // separate ref row (first) and consensus row (last, class sep-row)
  var pinned_top = rows.filter(r => r.classList.contains('ref-row'));
  var pinned_bot = rows.filter(r => r.classList.contains('sep-row'));
  var sortable   = rows.filter(r => !r.classList.contains('ref-row') && !r.classList.contains('sep-row'));

  var st  = _sortState[tableId] || {col: -1, dir: 1};
  var dir = (st.col === colIdx) ? -st.dir : 1;
  _sortState[tableId] = {col: colIdx, dir: dir};

  sortable.sort(function(a, b) {
    var av = _cellVal(a, colIdx, type);
    var bv = _cellVal(b, colIdx, type);
    if (av < bv) return -dir;
    if (av > bv) return  dir;
    return 0;
  });

  // Clear sort indicators on this table's headers
  tbl.querySelectorAll('th.sortable').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
  });
  var th = tbl.querySelectorAll('th.sortable')[
    Array.from(tbl.querySelectorAll('th')).indexOf(
      tbl.querySelectorAll('th')[colIdx])];
  // simpler: mark by data-col attribute
  tbl.querySelectorAll('th[data-col="' + colIdx + '"]').forEach(th => {
    th.classList.add(dir === 1 ? 'sort-asc' : 'sort-desc');
  });

  pinned_top.concat(sortable).concat(pinned_bot).forEach(r => tbody.appendChild(r));
}

function _cellVal(row, colIdx, type) {
  var td = row.cells[colIdx];
  if (!td) return '';
  var raw = (td.dataset.sort !== undefined) ? td.dataset.sort : td.innerText.trim();
  if (type === 'num') { var n = parseFloat(raw); return isNaN(n) ? -Infinity : n; }
  return raw.toLowerCase();
}
"""
