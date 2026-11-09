"""Self-contained, plain-language HTML report."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from jinja2 import Environment

from . import plots
from .explain import GLOSSARY, explain

TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>SeqGuide report – {{ r.sample }}</title>
<style>
:root{--ink:#0b0b0b;--ink2:#52514e;--line:#e4e3df;--bg:#fcfcfb;--card:#ffffff;--ok:#0b7a0b;--okbg:#e9f6e9;
--warn:#8a5a00;--warnbg:#fdf4dc;--bad:#a82828;--badbg:#fbe9e9;--info:#1c5cab;--infobg:#eaf2fc}
*{box-sizing:border-box}body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--ink);
background:var(--bg);margin:0;line-height:1.55}main{max-width:980px;margin:0 auto;padding:28px 20px 60px}
h1{font-size:24px;font-weight:600;margin:0 0 4px}h2{font-size:18px;font-weight:600;margin:34px 0 10px;padding-top:10px;border-top:1px solid var(--line)}
h3{font-size:15px;margin:18px 0 6px}.muted{color:var(--ink2)}.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;margin:10px 0}
.banner{border-radius:10px;padding:14px 16px;font-size:16px;margin:16px 0}.banner b{display:block;font-size:12px;letter-spacing:.04em;text-transform:uppercase;margin-bottom:2px}
.ok{background:var(--okbg);color:var(--ok)}.warn{background:var(--warnbg);color:var(--warn)}.block,.fail,.bad{background:var(--badbg);color:var(--bad)}.info{background:var(--infobg);color:var(--info)}
.msg{padding:8px 12px;border-radius:8px;margin:6px 0;font-size:14px}.msg .ico{font-weight:700;margin-right:6px}
table{border-collapse:collapse;width:100%;font-size:13.5px}th,td{text-align:left;padding:7px 8px;border-bottom:1px solid var(--line);vertical-align:top}
th{color:var(--ink2);font-weight:600;font-size:12.5px}td.num{text-align:right;font-variant-numeric:tabular-nums}
.pill{display:inline-block;padding:1px 8px;border-radius:999px;font-size:12px;font-weight:600;white-space:nowrap}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px}.stat{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.stat .v{font-size:20px;font-weight:600}.stat .k{font-size:12.5px;color:var(--ink2)}
details{margin:8px 0}summary{cursor:pointer;color:var(--info)}code,pre{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px}
pre{white-space:pre-wrap;background:#f4f4f1;padding:8px;border-radius:6px}img{max-width:100%}.flag{color:var(--warn);font-size:12.5px;display:block;margin-top:3px}
.why{font-size:13px;color:var(--ink2)}dl.gl dt{font-weight:600;margin-top:8px}dl.gl dd{margin:0 0 4px 0}
@media print{details{display:block}summary{display:none}}
</style></head><body><main>
<h1>Sequencing report: {{ r.sample }}</h1>
<div class="muted">{{ wf_title }} · run {{ r.started }} · SeqGuide {{ version }}{% if r.params.engine == 'builtin' %} · built-in engine{% endif %}</div>

{% if r.error %}
<div class="banner block"><b>The analysis stopped</b>{{ r.error }}</div>
<p class="why">What to do: check the warnings below. Too few reads, very low quality, or the wrong reference are the usual causes. The technical log is in <code>run.log</code> in the results folder.</p>
{% elif verdict %}
<div class="banner {{ verdict.level }}"><b>Bottom line</b>{{ verdict.text }}</div>
{% endif %}
{% for w in r.warnings %}<div class="msg warn"><span class="ico">!</span>{{ w }}</div>{% endfor %}

<h2>1 · Your data</h2>
<p class="why">This is what SeqGuide found in the files you gave it.</p>
{% set s = r.inspection %}
{% if s %}
<div class="grid">
<div class="stat"><div class="v">{{ layout_title }}</div><div class="k">Data type</div></div>
{% for f in s.files %}
<div class="stat"><div class="v">{{ "{:,}".format(f.estimated_total_reads or f.n_reads) }}{% if f.sampled %}*{% endif %}</div><div class="k">Reads in {{ f.name }}</div></div>
{% endfor %}
{% set f0 = s.files[0] %}
<div class="stat"><div class="v">{{ "%.0f"|format(f0.mean_len) }} bp</div><div class="k">Average read length</div></div>
{% if f0.mean_q is not none %}<div class="stat"><div class="v">Q{{ "%.0f"|format(f0.mean_q) }}</div><div class="k">Average quality</div></div>{% endif %}
{% if f0.plasmid_size_estimate %}<div class="stat"><div class="v">{{ "{:,}".format(f0.plasmid_size_estimate) }} bp</div><div class="k">Plasmid size from reads</div></div>{% endif %}
</div>
{% if s.files|selectattr("sampled")|list %}<p class="muted" style="font-size:12.5px">* estimated from the first part of the file.</p>{% endif %}
{% if s.pair_check %}<div class="msg {{ 'ok' if s.pair_check.ok else 'bad' }}">{{ s.pair_check.message }}</div>{% endif %}
{% if plot_len %}<img alt="Read length distribution" src="data:image/png;base64,{{ plot_len }}">{% endif %}
{% endif %}

<h2>2 · What we did and why</h2>
{% set rec = r.recommendation %}
{% if rec %}
<p><b>{{ rec.title }}.</b></p>
{% for m in rec.reasons %}<div class="msg {{ m.level }}"><span class="ico">✓</span>{{ m.text }}</div>{% endfor %}
{% for m in rec.messages %}{% if m.level != 'ok' %}<div class="msg {{ 'bad' if m.level=='block' else m.level }}"><span class="ico">{{ '!' if m.level in ('warn','block') else 'i' }}</span>{{ m.text }}</div>{% endif %}{% endfor %}
{% endif %}
<div class="card"><table><tr><th style="width:28%">Step</th><th>What it does</th><th style="width:9%">Result</th></tr>
{% for st in r.steps %}<tr><td>{{ st.step }}</td><td class="why">{{ st.explanation }}</td>
<td><span class="pill {{ 'ok' if st.status=='ok' else 'bad' }}">{{ 'done' if st.status=='ok' else 'failed' }}</span></td></tr>{% endfor %}
</table></div>

{% if r.contigs %}
<h2>3 · Assembly</h2>
<p class="why">The assembler builds sequences ("contigs") from the reads. For a plasmid, the ideal result is <b>one circular contig</b> of the expected size.</p>
<div class="card"><table><tr><th>Contig</th><th>Length</th><th>Circular?</th><th>Read depth</th><th>Note</th></tr>
{% for c in r.contigs[:15] %}<tr><td>{{ c.name }}</td><td class="num">{{ "{:,}".format(c.length) }} bp</td>
<td>{% if c.circular is none %}unknown{% elif c.circular %}<span class="pill ok">yes</span>{% else %}<span class="pill warn">no</span>{% endif %}</td>
<td class="num">{{ ("%.0f×"|format(c.mean_depth)) if c.mean_depth else "–" }}</td>
<td class="why">{% if c.name == r.main_contig %}{{ "Matches your reference" if r.comparison else "Main contig" }}{% if c.note %} · {{ c.note }}{% endif %}{% elif c.length < 1000 %}Small fragment; usually host DNA or a sequencing artefact (e.g. poly-G reads){% else %}Extra sequence: host (E. coli) DNA, or a second plasmid{% endif %}</td></tr>{% endfor %}
</table>
{% if r.contigs|length > 15 %}<p class="muted">… and {{ r.contigs|length - 15 }} more small contigs.</p>{% endif %}</div>
{% if r.contigs|length > 1 %}<p class="why">Several contigs are normal for Illumina plasmid preps: some DNA from the host bacterium is always sequenced too. Only the contig matching your plasmid matters.</p>{% endif %}
{% endif %}

{% if r.comparison %}
{% set c = r.comparison %}
<h2>4 · Comparison with your reference ({{ c.ref_name }})</h2>
<div class="grid">
<div class="stat"><div class="v">{{ "{:,}".format(c.qry_len) }} bp</div><div class="k">Your plasmid (assembled)</div></div>
<div class="stat"><div class="v">{{ "{:,}".format(c.ref_len) }} bp</div><div class="k">Reference</div></div>
<div class="stat"><div class="v">{{ "%.2f"|format(c.identity) }}%</div><div class="k">Identity where they align</div></div>
<div class="stat"><div class="v">{{ "%.1f"|format(c.ref_covered) }}%</div><div class="k">Reference found in your plasmid</div></div>
<div class="stat"><div class="v">{{ c.variants|length }}</div><div class="k">Differences</div></div>
</div>
{% for n in c.notes %}<div class="msg warn"><span class="ico">!</span>{{ n }}</div>{% endfor %}
{% if r.read_check %}<div class="msg ok"><span class="ico">✓</span>{{ "{:,}".format(r.read_check.reads) }} reads checked against the map: {{ r.read_check.pct_mapped }}% placed{% set cov = (r.coverage or {}).values()|list %}{% if cov %}, median depth {{ "%.0f"|format(cov[0].median) }}×, lowest {{ cov[0].min }}×{% endif %}.
{% if not r.read_check.extra_read_variants %} No position where reads disagree with the assembly (no sign of a mixed population).{% endif %}</div>{% endif %}
{% if c.strand == '-' %}<p class="why">The assembly was in the opposite orientation, so it was reverse-complemented. It was also rotated to start at the same position as your map. This changes nothing biologically.</p>{% endif %}
{% if plot_main %}<img alt="Plasmid map with differences and read depth" src="data:image/png;base64,{{ plot_main }}">{% endif %}
{% if c.variants %}
<div class="card"><table><tr><th>Position</th><th>Change</th><th>Feature</th><th>What it means</th><th>Reads</th></tr>
{% for v in c.variants %}<tr><td class="num">{{ "{:,}".format(v.position) }}</td><td>{{ v.description }}</td><td>{{ v.features|join(", ") or "–" }}</td>
<td><span class="pill {{ sev_class[v.severity] }}">{{ sev_word[v.severity] }}</span> {{ v.effect }}{% for fl in v.flags %}<span class="flag">⚠ {{ fl }}</span>{% endfor %}</td>
<td class="num">{{ v.depth if v.depth is not none else "–" }}</td></tr>{% endfor %}
</table></div>
{% else %}<p>No differences found.</p>{% endif %}
{% endif %}

{% if r.amplicon %}
{% set a = r.amplicon %}
<h2>3 · Targets</h2>
{% if a.mapping.mapped_pct is not none %}<div class="msg {{ 'ok' if a.mapping.mapped_pct >= 80 else 'warn' }}">{{ a.mapping.mapped_pct }}% of reads matched your target sequences.
{% if a.mapping.mapped_pct < 80 %} A low value can mean contamination, primer-dimers, or that the reference does not match the amplified region.{% endif %}</div>{% endif %}
<div class="card"><table><tr><th>Target</th><th>Length</th><th>Average depth</th><th>% positions ≥ {{ a.min_depth }}×</th><th>Variants</th><th>Status</th></tr>
{% for t in a.targets %}<tr><td>{{ t.name }}</td><td class="num">{{ "{:,}".format(t.length) }}</td><td class="num">{{ "%.0f"|format(t.coverage.mean) }}×</td>
<td class="num">{{ t.coverage.pct_ok }}%</td><td class="num">{{ t.variants|length }}</td>
<td><span class="pill {{ t.status }}">{{ {'ok':'good','warn':'partial','fail':'too low'}[t.status] }}</span></td></tr>{% endfor %}
</table></div>
{% if plot_main %}<img alt="Read depth along each target" src="data:image/png;base64,{{ plot_main }}">{% endif %}
<p class="why">Note: mutations under the PCR primers cannot be detected, because primer sequence comes from the primer, not from your sample.</p>
<h2>4 · Variants</h2>
{% for t in a.targets if t.variants %}
<h3>{{ t.name }}</h3>
<div class="card"><table><tr><th>Position</th><th>Change</th><th>Feature</th><th>What it means</th><th>Reads</th></tr>
{% for v in t.variants %}<tr><td class="num">{{ "{:,}".format(v.position) }}</td><td>{{ v.description }}</td><td>{{ v.features|join(", ") or "–" }}</td>
<td><span class="pill {{ sev_class[v.severity] }}">{{ sev_word[v.severity] }}</span> {{ v.effect }}{% for fl in v.flags %}<span class="flag">⚠ {{ fl }}</span>{% endfor %}</td>
<td class="num">{{ v.depth if v.depth is not none else "–" }}</td></tr>{% endfor %}
</table></div>
{% else %}<p>No variants found in any target.</p>{% endfor %}
{% endif %}

{% if plot_cov %}<h2>Read depth along the assembly</h2>
<p class="why">Even depth means the assembly is well supported. Sudden drops can point to errors or problem regions.</p>
<img alt="Read depth along the assembly" src="data:image/png;base64,{{ plot_cov }}">{% endif %}

<h2>Files in the results folder</h2>
<div class="card"><table>
{% for name, what in files %}<tr><td><code>{{ name }}</code></td><td class="why">{{ what }}</td></tr>{% endfor %}
</table></div>

<h2>Words used in this report</h2>
<dl class="gl">{% for g in glossary %}<dt>{{ g.title }}</dt><dd>{{ g.what }} <span class="why">{{ g.why }}</span></dd>{% endfor %}</dl>

<h2>Technical details</h2>
<details><summary>Tool versions and exact commands (for methods sections and reproducibility)</summary>
<h3>Tool versions</h3><table>{% for t, v in (r.tool_versions or {}).items() %}<tr><td>{{ t }}</td><td><code>{{ v }}</code></td></tr>{% endfor %}</table>
<h3>Commands</h3>{% for st in r.steps %}<p><b>{{ st.step }}</b>{% if st.seconds %} <span class="muted">({{ st.seconds }} s)</span>{% endif %}</p>
{% for cmd in st.commands %}<pre>{{ cmd }}</pre>{% endfor %}{% if st.error %}<pre>{{ st.error }}</pre>{% endif %}{% endfor %}
<h3>Settings</h3><pre>{{ params }}</pre>
</details>
</main></body></html>
"""

WF_TITLES = {
    "plasmid_denovo": "De novo plasmid assembly",
    "amplicon_align": "Amplicon / targeted gene alignment",
    "plasmid_align": "Alignment to plasmid reference",
    "compare_assembly": "Comparison of an existing assembly with its map",
}
LAYOUT_TITLES = {"paired": "Illumina paired-end", "single": "Illumina single-end", "long": "Nanopore",
                 "hybrid": "Illumina + Nanopore", "unknown": "Unclear"}
SEV_CLASS = {"high": "bad", "moderate": "warn", "low": "ok", "": "info"}
SEV_WORD = {"high": "likely harmful", "moderate": "check", "low": "minor", "": "note"}

FILES = [
    ("report.html", "This report."),
    ("assembly_matched_to_reference.fasta", "Your assembled plasmid, oriented and rotated to match your map. Open it in SnapGene/Benchling."),
    ("final_assembly.fasta", "All assembled contigs as produced."),
    ("differences.csv", "The table of differences, for Excel."),
    ("variants_explained.csv", "The table of variants with explanations, for Excel."),
    ("variants.vcf", "Variants in standard VCF format (for other tools)."),
    ("aligned.bam", "Reads aligned to your targets. Open in IGV to look at them yourself."),
    ("reads_on_assembly.bam", "Reads aligned to the assembly. Open in IGV together with final_assembly.fasta."),
    ("fastp.html", "Detailed read quality report from fastp."),
    ("provenance.json", "Machine-readable record of tools, versions and commands."),
    ("run.log", "Technical log."),
]


def _glossary_for(results: dict) -> list:
    keys = ["fastq", "depth", "mean_q"]
    insp = results.get("inspection") or {}
    layout = insp.get("layout")
    if layout in ("paired", "single", "long", "hybrid"):
        keys.append(layout)
    if results.get("workflow") == "plasmid_denovo":
        keys += ["de_novo", "circular", "comparison", "subsample", "trimming"]
    else:
        keys += ["alignment", "trimming", "allele_fraction"]
    if results.get("comparison") or results.get("amplicon"):
        keys += ["variant", "frameshift"]
    if layout in ("long", "hybrid"):
        keys += ["homopolymer", "plasmid_size"]
    if (results.get("reference") or {}).get("annotated"):
        keys.append("genbank")
    seen, out = set(), []
    for k in keys:
        if k in GLOSSARY and k not in seen:
            seen.add(k)
            out.append(explain(k))
    return out


def write_report(results: dict, plot_inputs: Optional[dict], path: Path) -> None:
    from . import __version__

    plot_len = plot_main = plot_cov = None
    insp = results.get("inspection") or {}
    try:
        if insp.get("files"):
            plot_len = plots.length_hist_figure(insp["files"])
        if plot_inputs:
            kind = plot_inputs.get("kind")
            if kind == "plasmid":
                plot_main = plots.plasmid_figure(plot_inputs["ref"], plot_inputs["variants"], plot_inputs.get("depth"))
            elif kind == "assembly" and plot_inputs.get("depth") is not None:
                plot_cov = plots.coverage_figure(plot_inputs["depth"], plot_inputs["name"])
            elif kind == "targets":
                plot_main = plots.targets_figure(plot_inputs["depth"], plot_inputs["refs"], plot_inputs["variants"],
                                                 plot_inputs["min_depth"])
    except Exception as exc:  # noqa: BLE001 - a broken plot must never kill the report
        results.setdefault("warnings", []).append(f"A figure could not be drawn: {exc}")

    verdict = None
    if results.get("comparison"):
        verdict = results["comparison"].get("verdict")
    elif results.get("amplicon"):
        verdict = results["amplicon"].get("verdict")
    elif results.get("contigs"):
        circ = [c for c in results["contigs"] if c.get("circular")]
        verdict = ({"level": "ok", "text": f"Assembly finished: {len(circ)} circular contig(s). Add a reference map to check for differences."}
                   if circ else {"level": "warn", "text": "Assembly finished, but no contig was confirmed circular. The plasmid may be incomplete."})

    outdir = Path(path).parent
    files = [(n, w) for n, w in FILES if (outdir / n).exists() or n == "report.html"]
    env = Environment(autoescape=True)
    tpl = env.from_string(TEMPLATE)
    out = tpl.render(
        r=results, version=__version__, wf_title=WF_TITLES.get(results.get("workflow"), results.get("workflow")),
        layout_title=LAYOUT_TITLES.get(insp.get("layout"), "–"), verdict=verdict, plot_len=plot_len,
        plot_main=plot_main, plot_cov=plot_cov, sev_class=SEV_CLASS, sev_word=SEV_WORD,
        files=files, glossary=_glossary_for(results), params=str(results.get("params", {})),
    )
    Path(path).write_text(out, encoding="utf-8")
