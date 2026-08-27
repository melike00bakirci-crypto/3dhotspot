"""Self-contained offline REVIEW_PACK dashboard — Output Contract IX.2.

Single file: inline CSS, inline vanilla JS only, figures embedded as base64 PNG
data URIs, decision-trace tables embedded inline. NO CDN, no external fonts, no
fetch()/XHR — those all fail under file:// on Windows. Expandable sections use
native <details>. Opens by double-click from an unzipped folder with no server
and no internet.

MAJOR and BLOCKING warnings are NEVER inside a collapsed <details> (IX.7).
Negative results render in the same prominent style as positive ones (IX.8).
"""
from __future__ import annotations

import base64
import html
from pathlib import Path

_CSS = """
:root{--bg:#ffffff;--fg:#16181d;--muted:#5b6270;--line:#e2e5ea;--card:#f7f8fa;
--ok:#1a7f4b;--warn:#8a6100;--bad:#a32020;--info:#1f5f9e;--accent:#2a3f5f}
@media (prefers-color-scheme:dark){:root{--bg:#14161a;--fg:#e8eaed;--muted:#9aa2b1;
--line:#2b303a;--card:#1c1f26;--ok:#4ec98a;--warn:#e0b350;--bad:#f08585;
--info:#79b8f3;--accent:#a9c2e8}}
*{box-sizing:border-box}
body{margin:0;padding:0 0 4rem;background:var(--bg);color:var(--fg);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1100px;margin:0 auto;padding:0 1.25rem}
header{border-bottom:1px solid var(--line);margin-bottom:1.5rem}
h1{font-size:1.6rem;margin:1.5rem 0 .25rem;letter-spacing:-.01em}
h2{font-size:1.15rem;margin:2rem 0 .75rem;padding-bottom:.35rem;
border-bottom:1px solid var(--line)}
h3{font-size:1rem;margin:1.25rem 0 .5rem;color:var(--muted)}
.sub{color:var(--muted);margin:0 0 1.25rem}
.banner{padding:.9rem 1.1rem;border-radius:8px;font-weight:600;margin:1rem 0;
border-left:5px solid}
.banner.COMPLETED{background:color-mix(in srgb,var(--ok) 12%,transparent);
border-color:var(--ok);color:var(--ok)}
.banner.COMPLETED_WITH_WARNINGS{background:color-mix(in srgb,var(--warn) 14%,transparent);
border-color:var(--warn);color:var(--warn)}
.banner.COMPLETED_NEGATIVE{background:color-mix(in srgb,var(--info) 12%,transparent);
border-color:var(--info);color:var(--info)}
.banner.PARTIAL{background:color-mix(in srgb,var(--warn) 14%,transparent);
border-color:var(--warn);color:var(--warn)}
.banner.FAILED,.banner.BLOCKED{background:color-mix(in srgb,var(--bad) 12%,transparent);
border-color:var(--bad);color:var(--bad)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:.75rem}
.tile{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:.75rem .9rem}
.tile .k{font-size:.72rem;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
.tile .v{font-size:1.3rem;font-weight:650;margin-top:.15rem;
font-variant-numeric:tabular-nums}
.tile .u{font-size:.75rem;color:var(--muted)}
.flags{display:flex;flex-wrap:wrap;gap:.4rem;margin:.75rem 0}
.flag{font-size:.75rem;padding:.25rem .55rem;border-radius:999px;font-weight:600;
border:1px solid}
.flag.BLOCKING{background:color-mix(in srgb,var(--bad) 15%,transparent);
border-color:var(--bad);color:var(--bad)}
.flag.MAJOR{background:color-mix(in srgb,var(--warn) 15%,transparent);
border-color:var(--warn);color:var(--warn)}
.flag.ADVISORY{background:color-mix(in srgb,var(--info) 12%,transparent);
border-color:var(--info);color:var(--info)}
.flag.INFO{background:var(--card);border-color:var(--line);color:var(--muted)}
table{border-collapse:collapse;width:100%;font-size:.82rem;
font-variant-numeric:tabular-nums}
th,td{text-align:left;padding:.35rem .55rem;border-bottom:1px solid var(--line);
white-space:nowrap}
th{position:sticky;top:0;background:var(--bg);font-weight:650;font-size:.75rem;
text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
tr.sel td{background:color-mix(in srgb,var(--ok) 14%,transparent);font-weight:650}
.scroll{overflow-x:auto;max-height:460px;overflow-y:auto;border:1px solid var(--line);
border-radius:8px}
details{background:var(--card);border:1px solid var(--line);border-radius:8px;
padding:.6rem .9rem;margin:.6rem 0}
summary{cursor:pointer;font-weight:600}
figure{margin:1rem 0}
figure img{width:100%;height:auto;border:1px solid var(--line);border-radius:8px;
background:#fff}
figcaption{color:var(--muted);font-size:.8rem;margin-top:.35rem}
code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.82rem}
pre{background:var(--card);border:1px solid var(--line);border-radius:8px;
padding:.75rem;overflow-x:auto}
.neg{border-left:5px solid var(--info);background:color-mix(in srgb,var(--info) 8%,transparent);
padding:.9rem 1.1rem;border-radius:8px;margin:1rem 0}
.muted{color:var(--muted)}
.nowrap{white-space:nowrap}
footer{margin-top:3rem;padding-top:1rem;border-top:1px solid var(--line);
color:var(--muted);font-size:.8rem}
"""

_JS = """
// Vanilla only. No fetch/XHR — everything needed is already in this file.
document.addEventListener('click', function(e){
  var t = e.target.closest('[data-filter]');
  if(!t) return;
  var table = document.getElementById(t.dataset.target);
  if(!table) return;
  var want = t.dataset.filter;
  table.querySelectorAll('tbody tr').forEach(function(tr){
    tr.style.display = (want === 'all' || tr.dataset.qc === want) ? '' : 'none';
  });
});
"""


def _b64_png(path: Path) -> str | None:
    if not path.is_file():
        return None
    return base64.b64encode(path.read_bytes()).decode("ascii")


def _tile(key: str, value, unit: str = "") -> str:
    shown = "— not computed" if value is None else html.escape(str(value))
    unit_html = f'<div class="u">{html.escape(unit)}</div>' if unit else ""
    return (f'<div class="tile"><div class="k">{html.escape(key)}</div>'
            f'<div class="v">{shown}</div>{unit_html}</div>')


def _table(rows: list[dict], columns: list[str], table_id: str,
           selected_key: str | None = None, max_rows: int = 400) -> str:
    if not rows:
        return '<p class="muted">No rows — table not produced for this run.</p>'
    head = "".join(f"<th>{html.escape(c)}</th>" for c in columns)
    body = []
    for r in rows[:max_rows]:
        is_sel = selected_key and str(r.get(selected_key)).upper() == "TRUE"
        qc = html.escape(str(r.get("qc_status", "")))
        cells = "".join(f"<td>{html.escape(str(r.get(c, 'NA')))}</td>" for c in columns)
        body.append(f'<tr class="{"sel" if is_sel else ""}" data-qc="{qc}">{cells}</tr>')
    more = ""
    if len(rows) > max_rows:
        more = (f'<p class="muted">Showing first {max_rows} of {len(rows)} rows — '
                f'the complete table is in FULL_RESULTS.</p>')
    return (f'<div class="scroll"><table id="{table_id}"><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>{more}')


def _figure(name: str, caption: str, b64: str | None) -> str:
    if b64 is None:
        return (f'<figure><figcaption class="muted">Figure <b>{html.escape(name)}</b> '
                f'NOT CREATED — see manifest.tsv for the reason.</figcaption></figure>')
    return (f'<figure><img alt="{html.escape(caption)}" src="data:image/png;base64,{b64}">'
            f'<figcaption>{html.escape(caption)}</figcaption></figure>')


def render_index_html(*, identity: dict, run_status: str, banner: str,
                      outcome_type: str, metrics: dict, flags: list[dict],
                      sections: dict, figures: list[tuple[str, str, Path]],
                      negative_result: dict | None,
                      warnings_rows: list[dict]) -> str:
    """Build the complete offline dashboard as one HTML string."""
    ident_tiles = "".join(_tile(k, v) for k, v in identity.items())

    flag_html = "".join(
        f'<span class="flag {html.escape(str(f["severity"]))}">'
        f'{html.escape(str(f["warning_code"]))}</span>' for f in flags
    ) or '<span class="flag INFO">no flags raised</span>'

    major = [w for w in warnings_rows
             if str(w.get("severity")) in ("MAJOR", "BLOCKING")]
    other = [w for w in warnings_rows
             if str(w.get("severity")) not in ("MAJOR", "BLOCKING")]

    # MAJOR/BLOCKING are rendered OPEN, never inside a collapsed <details>.
    major_html = "".join(
        f'<div class="banner {"FAILED" if w.get("severity")=="BLOCKING" else "PARTIAL"}">'
        f'[{html.escape(str(w.get("severity")))}] {html.escape(str(w.get("warning_code")))}'
        f'<div style="font-weight:400;margin-top:.3rem">'
        f'{html.escape(str(w.get("message")))}<br>'
        f'<span class="muted">Consequence: {html.escape(str(w.get("potential_consequence")))} '
        f'&middot; Action: {html.escape(str(w.get("recommended_action")))}</span></div></div>'
        for w in major
    ) or '<p class="muted">No MAJOR or BLOCKING warnings.</p>'

    other_html = _table(other, ["severity", "warning_code", "stage", "message",
                                "recommended_action"], "tbl-warn-other") if other else ""

    neg_html = ""
    if negative_result:
        neg_html = (
            f'<div class="neg"><b>NEGATIVE RESULT — {html.escape(str(negative_result.get("condition","")))}'
            f'</b><p>{html.escape(str(negative_result.get("detail","")))}</p>'
            f'<p class="muted">This is a complete, valid scientific outcome, not a '
            f'technical failure. It licenses no method modification.</p></div>')

    fig_html = "".join(_figure(n, c, _b64_png(p)) for n, c, p in figures)

    section_html = []
    for title, body in sections.items():
        section_html.append(f"<h2>{html.escape(title)}</h2>{body}")

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>hotspot3d — {html.escape(str(identity.get('Gene','run')))}</title>
<style>{_CSS}</style></head><body><div class="wrap">
<header>
  <h1>3D Hotspot &amp; Footprint Analysis — {html.escape(str(identity.get('Gene','')))}</h1>
  <p class="sub">Offline review pack. Everything on this page is embedded; no
  internet connection is used or required.</p>
</header>

<div class="banner {html.escape(run_status)}">{html.escape(banner)}
<div style="font-weight:400;margin-top:.25rem">outcome_type = {html.escape(outcome_type)}</div></div>

{neg_html}

<h2>Run identity</h2>
<div class="grid">{ident_tiles}</div>

<h2>Flags</h2>
<div class="flags">{flag_html}</div>

<h2>Headline metrics</h2>
<div class="grid">{''.join(_tile(k, v.get('value'), v.get('unit','')) for k, v in metrics.items())}</div>

{''.join(section_html)}

<h2>Figures</h2>
{fig_html}

<h2>QC / Warnings</h2>
{major_html}
{f'<details><summary>Other warnings ({len(other)})</summary>{other_html}</details>' if other else ''}

<footer>
  <p>Generated by hotspot3d. REVIEW_PACK is derived, never recomputed — every number
  here is copied verbatim from a FULL_RESULTS artifact and checksum-verified.</p>
  <p>Full audit trail: <code>FULL_RESULTS/</code> &middot; Reproducibility:
  <code>FULL_RESULTS/12_REPRODUCIBILITY/</code></p>
</footer>
</div><script>{_JS}</script></body></html>
"""
