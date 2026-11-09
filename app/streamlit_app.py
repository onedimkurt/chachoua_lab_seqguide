"""SeqGuide graphical interface: a step-by-step wizard.

Run with:  streamlit run app/streamlit_app.py
"""

from __future__ import annotations

import base64
import os
import sys
import time
from pathlib import Path

import streamlit as st

try:  # embedding a full HTML page; API location differs between Streamlit versions
    from streamlit.components.v1 import html as embed_html
except ImportError:  # pragma: no cover
    def embed_html(content, height=None, scrolling=True):
        st.html(content)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from seqguide import decide, inspect, plots, refmap, workflows  # noqa: E402
from seqguide.explain import GLOSSARY, explain  # noqa: E402

BROWSER = sys.platform == "emscripten"   # running inside the web browser (stlite / Pyodide)
DATA_DIR = Path(os.environ.get("SEQGUIDE_DATA", Path.home() / "SeqGuide"))
HOST_DIR = os.environ.get("SEQGUIDE_HOST_DIR", str(DATA_DIR))  # how the folder looks on the user's computer
INPUT, RUNS = DATA_DIR / "input", DATA_DIR / "runs"
INPUT.mkdir(parents=True, exist_ok=True)
RUNS.mkdir(parents=True, exist_ok=True)
EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".gb", ".gbk", ".genbank", ".dna")
STEPS = ["Choose data", "Understand data", "Decide", "Run", "Results"]
LEVEL_FN = {"ok": st.success, "info": st.info, "warn": st.warning, "block": st.error}

st.set_page_config(page_title="SeqGuide", page_icon="🧬", layout="wide")
ss = st.session_state
ss.setdefault("step", 0)


_DL_HTML = """
<button id="b" style="font:14px -apple-system,Segoe UI,Roboto,sans-serif;padding:7px 14px;border:1px solid #ccc;
 border-radius:8px;background:#fff;cursor:pointer">⬇ __LABEL__</button>
<script>
const bytes = Uint8Array.from(atob("__B64__"), c => c.charCodeAt(0));
const url = URL.createObjectURL(new Blob([bytes], {type: "__MIME__"}));
document.getElementById("b").onclick = () => {
  // create the link in the main page so the download is not blocked by the frame
  let d = document; try { if (window.parent && window.parent.document) d = window.parent.document; } catch (e) {}
  const a = d.createElement("a"); a.href = url; a.download = "__NAME__";
  d.body.appendChild(a); a.click(); a.remove();
};
</script>"""


def download(col, label: str, data: bytes, file_name: str, mime: str = "application/octet-stream") -> None:
    """Download button. In the browser version stlite's st.download_button cannot serve files, so a small
    button builds the file in memory (works for large files, unlike plain data: links)."""
    if BROWSER:
        html = (_DL_HTML.replace("__B64__", base64.b64encode(data).decode()).replace("__LABEL__", label)
                .replace("__MIME__", mime).replace("__NAME__", file_name.replace('"', "")))
        with col:
            embed_html(html, height=48)
    else:
        col.download_button(f"⬇ {label}", data, file_name=file_name, mime=mime)


def compact(n: int) -> str:
    return f"{n / 1e6:.2f} M" if n >= 1e6 else (f"{n / 1e3:.1f} k" if n >= 1e4 else f"{n:,}")


def go(step: int) -> None:
    ss.step = step
    st.rerun()


def tip(key: str) -> str:
    e = explain(key)
    return f"{e['what']} {e['why']}"


def show_msg(m) -> None:
    text = m.text if hasattr(m, "text") else m["text"]
    level = m.level if hasattr(m, "level") else m["level"]
    key = m.learn_more if hasattr(m, "learn_more") else m.get("learn_more")
    LEVEL_FN.get(level, st.info)(text, icon={"ok": "✅", "info": "ℹ️", "warn": "⚠️", "block": "⛔"}[level])
    if key and key in GLOSSARY:
        with st.expander(f"Learn more: {GLOSSARY[key][0]}"):
            st.write(GLOSSARY[key][1])
            st.caption(GLOSSARY[key][2])


# ------------------------------------------------------------------ sidebar
with st.sidebar:
    st.title("🧬 SeqGuide")
    st.caption("Guided analysis of plasmid and amplicon sequencing data")
    for i, name in enumerate(STEPS):
        mark = "✅" if i < ss.step else ("👉" if i == ss.step else "⬜")
        st.write(f"{mark} {i + 1}. {name}")
    st.divider()
    if st.button("Start over"):
        for k in list(ss.keys()):
            del ss[k]
        st.rerun()
    with st.expander("Glossary"):
        for k, (title, what, why) in GLOSSARY.items():
            st.markdown(f"**{title}** — {what} *{why}*")

# Everything below lives in one container that is replaced on every rerun, so content from the
# previous step never lingers (faded) below the current one.
with st.empty().container():
    # ------------------------------------------------------------------ step 1: choose data
    if ss.step == 0:
        st.header("1 · Which files do you want to analyse?")
        if BROWSER:
            st.markdown(
                "Add your files below. **They stay on this computer**: SeqGuide runs inside your web browser and "
                "nothing is uploaded to a server.  \n"
                "SeqGuide reads **FASTQ** (raw reads, `.fastq` / `.fastq.gz`), and references as **GenBank** (`.gb`), "
                "**SnapGene** (`.dna`) or **FASTA**."
            )
            ups = st.file_uploader("Choose or drop files", accept_multiple_files=True)
            if ups:
                new_files = [u for u in ups if not (INPUT / u.name).exists()]
                for u in new_files:
                    (INPUT / u.name).write_bytes(u.getbuffer())
        else:
            st.markdown(
                f"Put your files in the **`input`** folder inside **`{HOST_DIR}`** on your computer, or upload them here.  \n"
                "SeqGuide reads **FASTQ** (raw reads, `.fastq` / `.fastq.gz`), and references as **GenBank** (`.gb`), "
                "**SnapGene** (`.dna`) or **FASTA**."
            )
            with st.expander("Upload files (fine for small files; for large FASTQ files copy them into the folder instead)"):
                ups = st.file_uploader("Drop files here", accept_multiple_files=True)
                if ups and st.button("Save uploaded files"):
                    for u in ups:
                        (INPUT / u.name).write_bytes(u.getbuffer())
                    st.success(f"Saved {len(ups)} file(s).")
                    st.rerun()

        if st.button("🎓 Load practice data (a demo plasmid with 3 hidden mutations)"):
            from seqguide import demo
            with st.spinner("Creating practice files… (about 30 seconds)"):
                demo.make_demo(INPUT / "demo")
            st.success("Practice files added. Tick the two pDemo_illumina files and pDemo_expected_map.gb below.")

        files = sorted(p for p in INPUT.rglob("*") if p.is_file() and p.name.lower().endswith(EXTS))
        if not files:
            st.info("No files added yet." if BROWSER else "No sequencing files found yet in the input folder.")
        else:
            st.subheader("Files found")
            st.caption("Tick the reads of ONE sample (e.g. both R1 and R2) plus its reference map, if you have one.")
            chosen = []
            for p in files:
                rel = p.relative_to(INPUT)
                size = p.stat().st_size / 1e6
                if st.checkbox(f"{rel}  ({size:,.1f} MB)", key=f"f_{rel}"):
                    chosen.append(p)
            if st.button("Inspect selected files ➜", type="primary", disabled=not chosen):
                with st.spinner("Reading your files… (large files: up to a minute)"):
                    ss.dataset = inspect.inspect_dataset(chosen)
                go(1)

    # ------------------------------------------------------------------ step 2: understand
    elif ss.step == 1:
        ds = ss.dataset
        st.header("2 · What is in your data?")
        if not ds.samples:
            fastas = [Path(r.path) for r in ds.references if r.format == "fasta"]
            maps = [Path(r.path) for r in ds.references]
            if fastas and len(maps) >= 2:
                st.info("No raw reads (FASTQ) were selected, but you can still **compare an existing assembly with its map** "
                        "(for example your own earlier assembly, or one from a sequencing service).")
                asm = st.selectbox("Which file is your assembly?", fastas, format_func=lambda p: p.name)
                ref = st.selectbox("Which file is the expected map?", [m for m in maps if m != asm], format_func=lambda p: p.name)
                st.caption("Without reads, SeqGuide cannot check read depth, so differences are listed but not confirmed by reads.")
                if st.button("Compare ➜", type="primary"):
                    name = inspect.sample_name(asm.name)
                    ss.pop("results", None)
                    ss.pop("job", None)
                    ss.cfg = workflows.RunConfig(
                        workflow="compare_assembly", sample=name, outdir=RUNS / f"{name}_{time.strftime('%Y%m%d-%H%M%S')}",
                        reference=str(ref), params={"assembly": str(asm)})
                    go(3)
            else:
                st.error("No raw read files (FASTQ) were selected. SeqGuide needs FASTQ reads, or an assembly plus a map.")
            if st.button("⬅ Back"):
                go(0)
            st.stop()
        for o in ds.other:
            st.warning(f"{o.name}: " + " ".join(o.problems))

        names = [s.sample for s in ds.samples]
        if len(names) > 1:
            st.info(f"Your selection contains {len(names)} samples. Analyse one at a time; pick it here.")
        idx = st.selectbox("Sample", range(len(names)), format_func=lambda i: names[i]) if len(names) > 1 else 0
        sample = ds.samples[idx]
        ss.sample_idx = idx

        layout_key = sample.layout if sample.layout in GLOSSARY else None
        title = {"paired": "Illumina paired-end reads", "single": "Illumina single-end reads", "long": "Nanopore long reads",
                 "hybrid": "Illumina + Nanopore (hybrid)", "unknown": "Unclear file combination"}[sample.layout]
        st.subheader(f"Sample “{sample.sample}”: {title}")
        if layout_key:
            e = explain(layout_key)
            st.caption(f"{e['what']} {e['why']}")
        if sample.pair_check:
            (st.success if sample.pair_check["ok"] else st.error)(sample.pair_check["message"])

        for f in sample.files:
            st.markdown(f"**{f.name}** — {f.platform_reason}")
            c = st.columns(5)
            reads = f.estimated_total_reads or f.n_reads
            c[0].metric("Reads", compact(reads) + ("*" if f.sampled else ""), help=tip("n_reads"))
            c[1].metric("Average length", f"{f.mean_len:,.0f} bp", help=tip("read_length"))
            c[2].metric("Quality", f"Q{f.mean_q:.0f}" if f.mean_q is not None else "–", help=tip("mean_q"))
            if f.q30_percent is not None and f.platform in ("illumina", "short"):
                c[3].metric("Bases ≥ Q30", f"{f.q30_percent:.0f}%", help=tip("q30"))
            else:
                c[3].metric("N50", f"{f.n50:,} bp", help=tip("n50"))
            if f.plasmid_size_estimate:
                c[4].metric("Plasmid size (from reads)", f"{f.plasmid_size_estimate:,} bp", help=tip("plasmid_size"))
            else:
                c[4].metric("GC content", f"{f.gc_percent:.0f}%", help=tip("gc"))
            if f.sampled:
                st.caption("* estimated: only the first part of this large file was read for this overview.")
            if f.dimer_peak:
                st.info("A second peak at twice the plasmid size: plasmid dimers. Common and usually harmless.")
            for p in f.problems:
                st.warning(p)
        st.image(base64.b64decode(plots.length_hist_figure([f.to_dict() for f in sample.files])),
                 caption="How long the reads are. For whole-plasmid Nanopore data the main peak is your plasmid size.")
        if ds.references:
            st.write("📄 Reference file(s) found: " + ", ".join(f"**{r.name}**" for r in ds.references))

        c1, c2 = st.columns([1, 5])
        if c1.button("⬅ Back"):
            go(0)
        if c2.button("Continue ➜", type="primary"):
            go(2)

    # ------------------------------------------------------------------ step 3: decide
    elif ss.step == 2:
        ds = ss.dataset
        sample = ds.samples[ss.get("sample_idx", 0)]
        st.header("3 · What do you want to find out?")
        st.caption("Answer a few questions; SeqGuide recommends an analysis and explains why.")

        stype = st.radio("**What did you sequence?**", list(decide.SAMPLE_TYPES), format_func=decide.SAMPLE_TYPES.get,
                         key="q_type")
        goal = st.radio("**What is your goal?**", list(decide.GOALS), format_func=decide.GOALS.get, key="q_goal")
        reftype = st.radio("**Do you have the expected sequence (reference)?**", list(decide.REFERENCE_TYPES),
                           format_func=decide.REFERENCE_TYPES.get, key="q_ref",
                           help="The sequence you *expect*: your plasmid map from SnapGene/Benchling/Addgene, or the gene from NCBI.")
        ref_path, ref_len = None, None
        if reftype != "none":
            candidates = [Path(r.path) for r in ds.references]
            others = sorted(p for p in INPUT.rglob("*") if p.is_file() and p.suffix.lower() in (".gb", ".gbk", ".genbank", ".dna", ".fasta", ".fa", ".fna"))
            options = list(dict.fromkeys(candidates + others))
            if not options:
                st.error("No reference file found. Add your .gb / .dna / .fasta file in step 1." if BROWSER else
                         "No reference file found. Put your .gb / .dna / .fasta file in the input folder (step 1).")
            else:
                ref_path = st.selectbox("Reference file", options, format_func=lambda p: p.name)
                try:
                    refs = refmap.load_references(ref_path)
                    ref_len = sum(len(r.seq) for r in refs)
                    n_feat = sum(len(r.features) for r in refs)
                    st.caption(f"{len(refs)} sequence(s), {ref_len:,} bp in total, {n_feat} annotated feature(s).")
                    if reftype == "annotated" and n_feat == 0:
                        st.warning("This file has no annotated features; effects on genes cannot be predicted.")
                except Exception as exc:  # noqa: BLE001
                    st.error(f"Could not read this reference: {exc}")
                    ref_path = None

        default_size = ref_len or next((f.plasmid_size_estimate for f in sample.files if f.plasmid_size_estimate), 0) or 0
        size_label = "Expected plasmid size (bp)" if stype == "plasmid" else "Total length of your targets (bp)"
        size = st.number_input(size_label, min_value=0, value=int(default_size), step=100,
                               help="Used to estimate depth. Filled in from your reference when possible. 0 = unknown.")
        params_extra = {}
        have_asm = None
        if stype == "plasmid" and reftype != "none":
            if st.checkbox("I already have an assembly of this plasmid (skip assembling, just compare and check reads)"):
                fastas = sorted(p for p in INPUT.rglob("*") if p.is_file() and p.suffix.lower() in (".fasta", ".fa", ".fna"))
                if fastas:
                    have_asm = st.selectbox("Assembly file", fastas, format_func=lambda p: str(p.relative_to(INPUT)))
                else:
                    st.warning("No FASTA file found. Add your assembly in step 1.")
        if stype == "amplicon":
            ploidy = st.radio("**Where does the DNA come from?**", [2, 1],
                              format_func={2: "Genomic DNA of a diploid organism (human, mouse, plant…) — can be heterozygous",
                                           1: "Haploid or cloned DNA (plasmid, bacteria, yeast haploid, a single clone)"}.get,
                              help="This decides whether a variant in ~50% of reads is reported as heterozygous (diploid) "
                                   "or as a possible mixed population (haploid).")
            params_extra["ploidy"] = ploidy

        rec = decide.recommend(sample, stype, reftype, goal, expected_size=size or None)
        st.divider()
        st.subheader(f"Recommendation: {rec.title}")
        for m in rec.reasons:
            show_msg(m)
        for m in rec.messages:
            show_msg(m)

        workflow = rec.workflow
        if rec.alternatives and not rec.blocked:
            with st.expander("Other options (and why they are not recommended here)"):
                for alt in rec.alternatives:
                    st.markdown(f"**{alt.title}** — {alt.summary}")
                    st.markdown("👍 " + " · ".join(alt.pros) + "  \n👎 " + " · ".join(alt.cons))
                if st.checkbox(f"Use “{rec.alternatives[0].title}” instead (advanced)"):
                    workflow = rec.alternatives[0].workflow
                    st.warning("You chose a non-recommended option. That's fine. The report will note it.")

        threads = 1
        if not BROWSER:
            with st.expander("Advanced settings (defaults are fine for most users)"):
                threads = st.slider("CPU cores to use", 1, os.cpu_count() or 4, max(1, (os.cpu_count() or 2) - 1))

        if have_asm and not rec.blocked:
            workflow = "compare_assembly"
            params_extra["assembly"] = str(have_asm)
            st.info("Using your existing assembly: SeqGuide will compare it with the map and check it against your reads.")
        needs_ref = workflow in ("amplicon_align", "plasmid_align", "compare_assembly")
        can_run = workflow is not None and not rec.blocked and (ref_path or not needs_ref) and (reftype == "none" or ref_path)
        c1, c2 = st.columns([1, 5])
        if c1.button("⬅ Back"):
            go(1)
        if c2.button("Start analysis ➜", type="primary", disabled=not can_run):
            short = [f.path for f in sample.files if f.platform in ("illumina", "short")]
            long_ = [f.path for f in sample.files if f.platform in ("nanopore", "pacbio", "long")]
            params = dict(rec.params, **params_extra)
            if workflow == "plasmid_align":
                params["ploidy"] = 1
            ss.pop("results", None)
            ss.pop("run_error", None)
            ss.cfg = workflows.RunConfig(
                workflow=workflow, sample=sample.sample, outdir=RUNS / f"{sample.sample}_{time.strftime('%Y%m%d-%H%M%S')}",
                short_reads=short, long_reads=long_, reference=str(ref_path) if ref_path else None, params=params,
                threads=threads, inspection=sample.to_dict(), recommendation=rec.to_dict(),
                answers={"sample_type": stype, "goal": goal, "reference": reftype, "expected_size": size},
            )
            go(3)

    # ------------------------------------------------------------------ step 4: run
    elif ss.step == 3:
        cfg = ss.cfg
        st.header("4 · Running the analysis")
        st.caption("Keep this window open. Plasmids usually take 2–10 minutes; amplicons less.")
        if "job" not in ss and "results" not in ss and not ss.get("run_error"):
            try:
                ss.job = workflows.Job(cfg)
                ss.job_log = []
            except Exception as exc:  # noqa: BLE001  (e.g. tools missing for this workflow)
                ss.run_error = str(exc)
        if ss.get("run_error"):
            st.error(ss.run_error)
            if st.button("⬅ Back to questions"):
                ss.pop("run_error", None)
                go(2)
            st.stop()
        job = ss.job
        n = max(1, len(job.steps))
        st.progress(min(job.index, n) / n)
        if cfg.params.get("engine") == "builtin":
            st.caption("Using SeqGuide's built-in engine (no external tools needed). Its assembler gave base-for-base "
                       "the same results as Unicycler on the lab's validation plasmids.")
        for i, stp in enumerate(job.steps):
            mark = "✅" if i < job.index else ("⏳" if i == job.index else "⬜")
            st.markdown(f"{mark} **{stp.name}**")
            if i == job.index:
                st.info(stp.explanation + ("  \n⏳ Working… this step can take a minute or two in the browser."
                                           if BROWSER else ""))
        if ss.job_log:
            st.code("\n".join(ss.job_log[-10:]))
        if not job.done:
            # Two refreshes per step: the first only draws this progress page, the second runs the step.
            # (In the browser version the page can only redraw between script runs.)
            if ss.get("job_phase") != "run":
                ss.job_phase = "run"
                st.rerun()
            ss.job_phase = "show"
            job.set_log(lambda m: ss.job_log.append(m))
            job.step()
            st.rerun()
        ss.results = job.finish()
        ss.pop("job", None)
        go(4)

    # ------------------------------------------------------------------ step 5: results
    elif ss.step == 4:
        cfg = ss.cfg
        res = ss.results
        outdir = Path(cfg.outdir)
        st.header(f"5 · Results for {cfg.sample}")
        if BROWSER:
            st.caption("Download the files you need: they are not stored anywhere else and disappear when you close this tab.")
        else:
            rel = outdir.relative_to(DATA_DIR) if outdir.is_relative_to(DATA_DIR) else outdir
            st.caption(f"All files are saved in **{HOST_DIR}/{rel}**")
        report = outdir / "report.html"
        from seqguide import export
        zip_bytes = export.zip_results(outdir, cfg.sample)
        st.markdown("#### Download")
        download(st.container(), "Download everything (ZIP folder)", zip_bytes, f"{cfg.sample}_SeqGuide_results.zip", "application/zip")
        st.caption("One ZIP with the report, tables, your plasmid as an annotated GenBank file (SnapGene/Benchling), "
                   "and a **genome_browser** folder ready for IGV (FASTA + .fai, reads as BAM + .bai, session file). "
                   "Unzip it and see HOW_TO_OPEN.txt inside.")
        with st.expander("Download single files"):
            c = st.columns(3)
            if report.exists():
                download(c[0], "Report (HTML)", report.read_bytes(), f"{cfg.sample}_report.html", mime="text/html")
            for name, col in (("differences.csv", c[1]), ("variants_explained.csv", c[1]),
                              ("assembly_matched_to_reference.fasta", c[2])):
                p = outdir / name
                if p.exists():
                    download(col, name, p.read_bytes(), f"{cfg.sample}_{name}", "text/plain")
            gbs = sorted((outdir / "genome_browser").glob("*.gb"))
            for gbp in gbs:
                download(c[2], gbp.name + " (SnapGene)", gbp.read_bytes(), gbp.name, "text/plain")
        if report.exists():
            embed_html(report.read_text(encoding="utf-8"), height=1400, scrolling=True)
        if st.button("Analyse another sample"):
            for k in ("cfg", "results", "run_error", "job", "job_log"):
                ss.pop(k, None)
            go(0)
