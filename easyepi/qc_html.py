"""A single-file HTML QC report: figures and numbers in one file.

Why one file
------------
A QC result is something a person sends to a supervisor, attaches to a ticket, or
reopens months later. A directory of JSON plus a directory of PDFs is none of
those: the PDFs arrive detached from the numbers, and the reviewer cannot tell
which figure belongs to which sample. Figures are inlined as base64 so the report
has no external references -- it opens from an email attachment, a network share,
or a ``file://`` URL with no assumptions.

Image sizing, and why every <img> carries explicit pixel dimensions
-------------------------------------------------------------------
CSS alone is not enough here. ``max-width: 100%`` resolves against the containing
block, and a containing block whose width is not yet established (a preview pane,
an email client, a print layout, a flex child without a basis) resolves it to
zero -- the image then renders at 0x0 and the report looks broken while the file
itself is fine. The intrinsic pixel size is therefore read from each image header
and written as ``width``/``height`` attributes, with ``width:100%; height:auto``
in CSS to let it shrink. An explicit non-zero width cannot collapse to zero.

What this report deliberately does not do
-----------------------------------------
* **It does not invent a pass/fail.** The badge shows the verdict the threshold
  layer already produced; nothing here re-judges a metric.
* **It does not turn "not measured" into a zero.** A metric that was never
  evaluated renders as a dash, never as ``0``.
* **It does not print a sentinel as a measurement.** ``assay_metrics`` writes
  ``tss_enrichment = 0.0`` when the ratio was undefined, so a floor threshold can
  still fire; the table shows absence instead of that sentinel.
* **It does not hide the caveats.** Anything the payload recorded as a problem is
  printed on the sample card rather than left in the JSON.

No new dependency: string formatting, ``html.escape`` and a small image-header
reader. jinja2 would be a third-party install for one module.
"""

from __future__ import annotations

import base64
import html
import struct
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

__all__ = ["render_qc_html", "render_merged_html", "fmt_value", "safe_text",
           "image_size"]

#: Shown where a value exists in the schema but was never computed. Never "0".
_NOT_MEASURED = "\u2014"

_VERDICT_STYLE = {
    "pass": ("pass", "ok"),
    "warn": ("warn", "warn"),
    "not_evaluated": ("not evaluated", "muted"),
    "incomplete": ("incomplete", "warn"),
}

#: Metrics worth putting in the headline table, in reading order, with a label
#: and the unit a reader needs to interpret the number.
_HEADLINE = (
    ("mapping_rate", "mapping rate", "\u00d7"),
    ("duplicate_rate", "duplicate rate", "\u00d7"),
    ("mitochondrial_fraction", "mitochondrial fraction", "\u00d7"),
    ("frip", "FRiP", "\u00d7"),
    ("signal_enrichment", "peak / flank signal", "\u00d7"),
    ("tss_enrichment", "TSS enrichment", "\u00d7"),
    ("tss_enrichment_bigwig", "TSS enrichment (BigWig)", "\u00d7"),
    ("median_fragment_length", "median fragment", "bp"),
    ("fraction_nfr", "fraction NFR", "\u00d7"),
    ("fraction_mononucleosome", "fraction mono-nucleosome", "\u00d7"),
    ("fraction_peaks_zero", "peaks with zero reads", "\u00d7"),
    ("blacklist_fraction", "blacklist fraction", "\u00d7"),
    ("mean_beta", "mean beta", "\u00d7"),
    ("nrl_bp", "nucleosome repeat length", "bp"),
)

#: Metrics whose 0.0 may be a stand-in for "undefined" rather than a measurement.
_SYNTHESISED_ZERO = frozenset({"tss_enrichment"})

#: String-valued entries that are genuinely per-sample state worth a table row.
_TEXT_METRICS = frozenset({"tss_enrichment_status", "tss_enrichment_source",
                           "interpretable", "duplicate_source"})

_STATE_HINTS = {
    "background_zero": "flanking pileup was zero, so the ratio is undefined",
    "no_signal": "no pileup at the TSS either",
    "no_reads": "no usable reads near any TSS",
    "no_tss": "no usable TSS positions",
    "error": "the reader could not read its input",
}

_CSS = """
*,*::before,*::after{box-sizing:border-box}
:root{
  --bg:#f6f7f9; --card:#fff; --ink:#14161a; --ink-2:#5b6270; --ink-3:#8a919e;
  --line:#e4e7ec; --line-strong:#cfd4dc;
  --ok-bg:#e7f6ec; --ok-ink:#136c33; --warn-bg:#fdf1dc; --warn-ink:#8a5200;
  --muted-bg:#eef0f3; --muted-ink:#5b6270; --bad-bg:#fdeceb; --bad-ink:#a1251b;
  --accent:#2f6fd0; --radius:10px;
}
html{-webkit-text-size-adjust:100%}
body{margin:0;padding:28px 20px 64px;background:var(--bg);color:var(--ink);
  font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",
  Arial,"Noto Sans",sans-serif;}
.wrap{max-width:1080px;margin:0 auto}
header.top{margin-bottom:20px}
h1{font-size:24px;line-height:1.25;margin:0 0 6px;letter-spacing:-.01em}
.sub{color:var(--ink-2);font-size:13px}
.sub .dot{margin:0 7px;color:var(--ink-3)}
h2{font-size:15px;margin:28px 0 10px;color:var(--ink);letter-spacing:.01em}
h2 .count{color:var(--ink-3);font-weight:400;margin-left:6px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(112px,1fr));
  gap:10px;margin:16px 0 4px}
.tile{background:var(--card);border:1px solid var(--line);border-radius:var(--radius);
  padding:12px 14px}
.tile .n{display:block;font-size:24px;font-weight:650;line-height:1.15;
  font-variant-numeric:tabular-nums}
.tile .k{display:block;font-size:11px;color:var(--ink-2);text-transform:uppercase;
  letter-spacing:.06em;margin-top:3px}
.tile.ok .n{color:var(--ok-ink)} .tile.warn .n{color:var(--warn-ink)}
.card{background:var(--card);border:1px solid var(--line);border-radius:var(--radius);
  padding:18px 20px;margin:12px 0}
.card-head{display:flex;align-items:center;gap:10px;flex-wrap:wrap;
  padding-bottom:12px;margin-bottom:4px;border-bottom:1px solid var(--line)}
.card-head h3{font-size:15px;margin:0;font-weight:650}
.badge{font-size:10.5px;font-weight:700;padding:3px 9px;border-radius:999px;
  text-transform:uppercase;letter-spacing:.06em;white-space:nowrap}
.badge.ok{background:var(--ok-bg);color:var(--ok-ink)}
.badge.warn{background:var(--warn-bg);color:var(--warn-ink)}
.badge.muted{background:var(--muted-bg);color:var(--muted-ink)}
.assay{margin-left:auto;color:var(--ink-3);font-size:12px}
table{border-collapse:collapse;width:100%;margin:6px 0 0;font-size:13px}
th,td{text-align:left;padding:7px 12px;border-bottom:1px solid var(--line)}
tr:last-child td{border-bottom:0}
th{font-size:10.5px;text-transform:uppercase;letter-spacing:.06em;
  color:var(--ink-2);font-weight:650;white-space:nowrap;background:#fafbfc}
th.sortable{cursor:pointer;user-select:none}
th.sortable:hover{color:var(--ink)}
th.sortable::after{content:"";display:inline-block;width:0}
td.num{text-align:right;font-variant-numeric:tabular-nums}
td.na{color:var(--ink-3);font-style:italic}
td.metric{white-space:nowrap}
td.unit{color:var(--ink-3);font-size:11.5px}
tr.breach td{background:var(--bad-bg)}
tr.breach td.metric{font-weight:650;color:var(--bad-ink)}
.note{border-left:3px solid var(--warn-ink);background:#fffdf7;padding:9px 13px;
  margin:10px 0;font-size:12.5px;border-radius:0 6px 6px 0;color:#4a3a12}
.note.bad{border-left-color:var(--bad-ink);background:#fff8f7;color:#5c1f19}
.note code{background:rgba(0,0,0,.06);padding:1px 5px;border-radius:4px;font-size:12px}
.figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));
  gap:14px;margin-top:14px}
figure{margin:0}
figure img{display:block;width:100%;height:auto;border:1px solid var(--line);
  border-radius:8px;background:#fff}
figcaption{font-size:11.5px;color:var(--ink-2);margin-top:6px}
.problems{border-left:3px solid var(--bad-ink);background:#fff8f7;padding:10px 14px;
  margin:12px 0;font-size:12.5px;border-radius:0 6px 6px 0}
.problems ul{margin:6px 0 0;padding-left:18px}
.empty{color:var(--ink-3);font-style:italic}
footer{margin-top:36px;padding-top:14px;border-top:1px solid var(--line);
  font-size:11.5px;color:var(--ink-3)}
@media print{body{background:#fff;padding:0}.card,.tile{break-inside:avoid;
  border-color:#d8dce3}}
"""

_JS = """
document.querySelectorAll('table.sortable').forEach(function (table) {
  var heads = table.tHead ? table.tHead.rows[0].cells : [];
  Array.prototype.forEach.call(heads, function (th, i) {
    if (!th.dataset.sortable) { th.dataset.sortable = '1'; }
    th.classList.add('sortable');
    th.addEventListener('click', function () {
      var body = table.tBodies[0];
      var rows = Array.prototype.slice.call(body.rows);
      var asc = th.dataset.asc !== '1';
      th.dataset.asc = asc ? '1' : '0';
      rows.sort(function (a, b) {
        var x = a.cells[i] ? (a.cells[i].dataset.v !== undefined
                 ? a.cells[i].dataset.v : a.cells[i].textContent) : '';
        var y = b.cells[i] ? (b.cells[i].dataset.v !== undefined
                 ? b.cells[i].dataset.v : b.cells[i].textContent) : '';
        var nx = parseFloat(x), ny = parseFloat(y);
        if (!isNaN(nx) && !isNaN(ny)) { return asc ? nx - ny : ny - nx; }
        return asc ? String(x).localeCompare(String(y))
                   : String(y).localeCompare(String(x));
      });
      rows.forEach(function (r) { body.appendChild(r); });
    });
  });
});
"""


def safe_text(value: Any) -> str:
    """Escape for HTML, keeping the distinction between None and empty string."""
    if value is None:
        return ""
    return html.escape(str(value), quote=True)


def image_size(blob: bytes) -> tuple[int, int] | None:
    """Intrinsic pixel size of a PNG/GIF/JPEG, or ``None`` if not recognised.

    Needed because ``width:100%`` alone can resolve to zero in a viewer whose
    containing block is not laid out yet, which renders the figure at 0x0 and
    makes a fine report look broken.
    """
    if len(blob) >= 24 and blob[:8] == b"\x89PNG\r\n\x1a\n" and blob[12:16] == b"IHDR":
        w, h = struct.unpack(">II", blob[16:24])
        return int(w), int(h)
    if blob[:6] in (b"GIF87a", b"GIF89a") and len(blob) >= 10:
        w, h = struct.unpack("<HH", blob[6:10])
        return int(w), int(h)
    if blob[:2] == b"\xff\xd8":  # JPEG: walk to a Start-Of-Frame marker
        i = 2
        while i + 9 < len(blob):
            if blob[i] != 0xFF:
                i += 1
                continue
            marker = blob[i + 1]
            if marker in (0xC0, 0xC1, 0xC2, 0xC3, 0xC5, 0xC6, 0xC7,
                          0xC9, 0xCA, 0xCB, 0xCD, 0xCE, 0xCF):
                h, w = struct.unpack(">HH", blob[i + 5:i + 9])
                return int(w), int(h)
            if marker in (0xD8, 0xD9) or 0xD0 <= marker <= 0xD7:
                i += 2
                continue
            seg = struct.unpack(">H", blob[i + 2:i + 4])[0]
            i += 2 + seg
        return None
    return None


def fmt_value(value: Any, *, digits: int = 4) -> tuple[str, str]:
    """Format one metric for a table cell; return ``(text, css_class)``.

    A missing value is rendered as a dash with a muted class, never as ``0``:
    the report is skimmed, and a skimming reader cannot tell a real zero from an
    absent measurement once both look like ``0``.
    """
    if value is None:
        return _NOT_MEASURED, "na"
    if isinstance(value, bool):
        return ("yes" if value else "no"), ""
    if isinstance(value, float):
        if value != value:  # NaN
            return _NOT_MEASURED, "na"
        if value.is_integer():
            return str(int(value)), "num"
        return f"{value:.{digits}g}", "num"
    if isinstance(value, int):
        return str(value), "num"
    return safe_text(value), ""


def _data_uri(path: Path) -> tuple[str, tuple[int, int] | None] | None:
    suffix = path.suffix.lower().lstrip(".")
    mime = {"png": "image/png", "svg": "image/svg+xml", "jpg": "image/jpeg",
            "jpeg": "image/jpeg", "gif": "image/gif"}.get(suffix)
    if mime is None:
        return None
    try:
        blob = path.read_bytes()
    except OSError:
        return None
    if not blob:
        return None
    return (f"data:{mime};base64,{base64.b64encode(blob).decode('ascii')}",
            image_size(blob))


def _collect_figures(out_dir: Path | None) -> dict[str, list[tuple[str, str, tuple[int, int] | None]]]:
    """Group figures under ``figures/`` by filename stem.

    PDF is skipped: browsers embed nothing useful from it without a plugin, and the
    point is a self-contained page.
    """
    found: dict[str, list[tuple[str, str, tuple[int, int] | None]]] = {}
    if out_dir is None:
        return found
    fig_dir = Path(out_dir) / "figures"
    if not fig_dir.is_dir():
        return found
    for path in sorted(fig_dir.iterdir()):
        if not path.is_file():
            continue
        got = _data_uri(path)
        if got is None:
            continue
        uri, size = got
        found.setdefault(path.stem, []).append((path.suffix.lstrip(".").upper(), uri, size))
    return found


def _collect_from_paths(paths) -> dict[str, list[tuple[str, str, tuple[int, int] | None]]]:
    """Group an explicit list of figure files by filename stem.

    Preferred over scanning `figures/` because a shared output directory holds the
    plots of every run that wrote into it: a MNase report would otherwise embed an
    ATAC run's charts. The caller passes exactly what this run produced.
    """
    found: dict[str, list[tuple[str, str, tuple[int, int] | None]]] = {}
    for raw_path in paths or ():
        path = Path(raw_path)
        if not path.is_file():
            continue
        got = _data_uri(path)
        if got is None:
            continue
        uri, size = got
        found.setdefault(path.stem, []).append((path.suffix.lstrip(".").upper(), uri, size))
    return found


def _verdict_badge(verdict: Mapping[str, Any] | None) -> str:
    overall = (verdict or {}).get("overall") or "not_evaluated"
    label, css = _VERDICT_STYLE.get(overall, (str(overall), "muted"))
    return f'<span class="badge {css}">{safe_text(label)}</span>'


def _metric_rows(flat: Mapping[str, Any]) -> list[tuple[str, str, Any]]:
    """``(key, label, value)`` for the headline metrics, then anything else numeric.

    A headline metric present but ``None`` is kept and printed as a dash: dropping
    it would silently shorten the table, and a reader comparing two reports could
    not tell "this sample has no TSS enrichment" from "this report has no TSS row".
    """
    units = {key: unit for key, _label, unit in _HEADLINE}
    labels = {key: label for key, label, _unit in _HEADLINE}
    rows: list[tuple[str, str, Any]] = []
    shown: set[str] = set()
    for key, label, _unit in _HEADLINE:
        if key in flat:
            rows.append((key, label, flat[key]))
            shown.add(key)
    for key in sorted(flat):
        if key in shown or key.startswith("_") or key in _TEXT_METRICS:
            continue
        value = flat[key]
        if isinstance(value, (int, float, bool)) or value is None:
            rows.append((key, labels.get(key, key.replace("_", " ")), value))
    return [(k, lab, v) for k, lab, v in rows], units


def _sample_card(sample_id: str, payload: Mapping[str, Any]) -> str:
    assay_qc = payload.get("assay_qc") or {}
    bam_qc = payload.get("bam_qc") or {}
    verdict = payload.get("qc_verdict") or {}

    flat: dict[str, Any] = {}
    for source in (bam_qc, assay_qc):
        for key, value in source.items():
            if isinstance(value, (int, float, bool)) or value is None:
                flat.setdefault(key, value)
            elif isinstance(value, str) and key in _TEXT_METRICS:
                flat.setdefault(key, value)

    breaches = {c.get("metric") for c in verdict.get("checks", [])
                if isinstance(c, dict) and c.get("status") == "breached"}

    rows_pairs, units = _metric_rows(flat)
    body_rows = []
    for key, label, value in rows_pairs:
        # A 0.0 that `assay_metrics` synthesised for a ratio it could not form is
        # not a measurement; the status sits right there saying so.
        if (value == 0 and key in _SYNTHESISED_ZERO
                and flat.get(f"{key}_status") in ("background_zero", "no_signal")):
            value = None
        text, align = fmt_value(value)
        cls = "breach" if key in breaches else ""
        data = "" if value is None else f' data-v="{safe_text(value)}"'
        unit = units.get(key, "")
        body_rows.append(
            f'<tr class="{cls}"><td class="metric">{safe_text(label)}</td>'
            f'<td class="{align}"{data}>{text}</td>'
            f'<td class="unit">{safe_text(unit)}</td></tr>')

    table = (f'<table class="sortable"><thead><tr><th>metric</th>'
             f'<th>value</th><th>unit</th></tr></thead><tbody>'
             + "".join(body_rows) + "</tbody></table>") if body_rows else (
        '<p class="empty">No assay metrics were computed for this sample: '
        'the inputs needed for them were not supplied.</p>')

    notes = []
    warned = verdict.get("warnings") or []
    if warned:
        notes.append(f'<div class="note bad"><strong>Breached:</strong> '
                     f'{safe_text(", ".join(str(w) for w in warned))}</div>')
    status = assay_qc.get("tss_enrichment_status")
    if status and status != "measured":
        hint = _STATE_HINTS.get(str(status), "")
        notes.append(f'<div class="note">TSS enrichment not reported: '
                     f'<code>{safe_text(status)}</code>'
                     + (f' &mdash; {safe_text(hint)}' if hint else "") + "</div>")
    disagreement = assay_qc.get("tss_enrichment_disagreement")
    if isinstance(disagreement, Mapping):
        notes.append(
            f'<div class="note">Two TSS definitions disagree by '
            f'{safe_text(disagreement.get("ratio"))}&times; '
            f'(BAM {safe_text(disagreement.get("bam"))}, '
            f'BigWig {safe_text(disagreement.get("bigwig"))}) &mdash; the windows '
            f'differ; <code>tss_enrichment</code> is the BAM definition.</div>')
    missing_core = verdict.get("missing_core_metrics") or []
    if missing_core:
        notes.append(f'<div class="note">Core metrics not evaluated: '
                     f'{safe_text(", ".join(str(m) for m in missing_core))}</div>')
    coverage = verdict.get("core_coverage")
    if coverage is not None:
        notes.append(f'<div class="note">Core-metric coverage: '
                     f'{safe_text(coverage)}</div>')

    assay_label = payload.get("assay")
    return f"""
<section class="card">
  <div class="card-head">
    <h3>{safe_text(sample_id)}</h3>
    {_verdict_badge(verdict)}
    {f'<span class="assay">{safe_text(assay_label)}</span>' if assay_label else ''}
  </div>
  {''.join(notes)}
  {table}
</section>"""


def _figures_section(figures: Mapping[str, list]) -> str:
    """One run-level figure grid.

    `qc` writes its figures to a single `figures/` directory for the whole run, not
    one directory per sample, so repeating them inside every sample card would show
    the same plot N times and imply it belonged to that sample.
    """
    if not figures:
        return ""
    items_html = []
    for stem, items in figures.items():
        label = stem.replace("_", " ").strip()
        for _ext, uri, size in items:
            # Explicit pixel dimensions: a percentage width can resolve to zero in
            # a pane whose containing block is not laid out yet, which renders the
            # figure at 0x0 and makes an intact report look broken.
            dims = (f' width="{size[0]}" height="{size[1]}"' if size
                    else ' width="960" height="640"')
            items_html.append(
                f'<figure><img alt="{safe_text(label)}" src="{uri}"{dims}>'
                f'<figcaption>{safe_text(label)}</figcaption></figure>')
    return ('<h2>Figures<span class="count">run level</span></h2>'
            f'<div class="figs">{"".join(items_html)}</div>')


def _summary_tiles(summary: Mapping[str, Any]) -> str:
    counts = summary.get("samples") or {}
    tiles = []
    for key, css in (("pass", "ok"), ("warn", "warn"),
                     ("not_evaluated", ""), ("incomplete", "warn")):
        if key in counts:
            tiles.append(f'<div class="tile {css}"><span class="n">'
                         f'{safe_text(counts[key])}</span><span class="k">'
                         f'{safe_text(key.replace("_", " "))}</span></div>')
    warned = summary.get("warning_metrics") or {}
    if warned:
        tiles.append(f'<div class="tile"><span class="n">{len(warned)}</span>'
                     f'<span class="k">metrics warned</span></div>')
    return f'<div class="tiles">{"".join(tiles)}</div>' if tiles else ""


def _page(title: str, subtitle: str, body: str) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{safe_text(title)}</title>
<style>{_CSS}</style>
</head><body><div class="wrap">
<header class="top">
  <h1>{safe_text(title)}</h1>
  <div class="sub">{subtitle}</div>
</header>
{body}
<footer>
  Generated by easyepi on {safe_text(stamp)}. Figures are embedded, so this file
  is self-contained. A dash means the metric was not measured, not that it was
  zero. Verdicts come from the threshold table; this page does not re-judge them.
</footer>
</div><script>{_JS}</script></body></html>
"""


def render_qc_html(payload: Mapping[str, Any], out_path, *,
                   figure_paths=None, figure_dir=None,
                   title: str = "easyepi QC report") -> Path:
    """Write the QC report for one run and return its path.

    ``figure_paths`` is the exact list of files this run produced and is the
    preferred input. ``figure_dir`` remains for callers that have a directory of
    figures and no manifest; it scans that directory, which is only safe when the
    directory belongs to this run alone.
    """
    samples = payload.get("samples") or {}
    summary = payload.get("verdict_summary") or {}

    if figure_paths is not None:
        figures = _collect_from_paths(figure_paths)
        suite_dir = None
        figure_dir = None
    else:
        suite_dir = figure_dir
    if figure_paths is None and isinstance(figure_dir, Mapping):
        # Per-sample directories: merge them under a "<sample> — <stem>" label so
        # one grid still shows every figure exactly once.
        figures: dict[str, list] = {}
        for sample_id, directory in figure_dir.items():
            for stem, items in _collect_figures(
                    Path(directory) if directory else None).items():
                figures[f"{sample_id} {stem}".replace("_", " ")] = items
    elif figure_paths is None:
        figures = _collect_figures(Path(suite_dir) if suite_dir else None)

    cards = [_sample_card(str(sid), dict(s) if isinstance(s, Mapping) else {})
             for sid, s in samples.items()]

    meta = [f"{len(samples)} sample" + ("" if len(samples) == 1 else "s")]
    if payload.get("assay"):
        meta.append(f"assay {payload['assay']}")
    if payload.get("reports"):
        meta.append(f"{len(payload['reports'])} source report(s)")
    subtitle = '<span class="dot">\u00b7</span>'.join(safe_text(m) for m in meta)

    body = (_summary_tiles(summary)
            + f'<h2>Samples<span class="count">{len(samples)}</span></h2>'
            + (''.join(cards) or '<p class="empty">No samples in this report.</p>')
            + _figures_section(figures))
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(_page(title, subtitle, body), encoding="utf-8")
    return out_path


def render_merged_html(merged: Mapping[str, Any], out_path, *,
                       title: str = "easyepi batch QC report") -> Path:
    """Batch-level report: every sample from every merged report in one table."""
    samples = merged.get("samples") or {}
    summary = merged.get("verdict_summary") or {}

    rows = []
    for sample_id, sample in sorted(samples.items()):
        sample = sample if isinstance(sample, Mapping) else {}
        verdict = sample.get("qc_verdict") or {}
        warned = verdict.get("warnings") or []
        overall = verdict.get("overall") or "not_evaluated"
        label, css = _VERDICT_STYLE.get(overall, (str(overall), "muted"))
        cls = "breach" if warned else ""
        rows.append(
            f'<tr class="{cls}">'
            f'<td data-v="{safe_text(sample_id)}">{safe_text(sample_id)}</td>'
            f'<td data-v="{safe_text(sample.get("batch"))}">'
            f'{safe_text(sample.get("batch"))}</td>'
            f'<td>{safe_text(sample.get("assay"))}</td>'
            f'<td><span class="badge {css}">{safe_text(label)}</span></td>'
            f'<td data-v="{len(warned)}">'
            f'{safe_text(", ".join(str(w) for w in warned))}</td></tr>')

    problems = merged.get("problems") or []
    problem_html = ""
    if problems:
        items = "".join(f"<li>{safe_text(p)}</li>" for p in problems)
        problem_html = (f'<div class="problems"><strong>{len(problems)} report(s) '
                        f'skipped</strong><ul>{items}</ul></div>')

    table = f"""<table class="sortable">
<thead><tr><th>sample</th><th>batch</th><th>assay</th><th>verdict</th>
<th>warnings</th></tr></thead>
<tbody>{''.join(rows) or '<tr><td colspan="5" class="empty">no samples</td></tr>'}
</tbody></table>"""

    subtitle = (f'{safe_text(merged.get("n_reports"))} report(s)'
                f'<span class="dot">\u00b7</span>'
                f'{safe_text(merged.get("n_samples"))} sample(s)')
    body = (_summary_tiles(summary) + problem_html
            + '<h2>Samples</h2>' + table)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(_page(title, subtitle, body), encoding="utf-8")
    return out_path
