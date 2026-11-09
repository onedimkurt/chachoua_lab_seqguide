"""Command-line entry points (used for batch runs and validation; the GUI is the main interface).

  python -m seqguide check                      # are all tools installed? prints versions
  python -m seqguide batch samples.tsv [--out DIR] [--threads N]

samples.tsv (tab-separated, header required):
  sample  sample_type  reference  reads  [assembly]  [workflow]  [expected_size]  [engine]
  - sample_type: plasmid | amplicon
  - reads: one or more FASTQ paths separated by ';' (R1;R2, or a Nanopore file, or both)
  - workflow (optional): force plasmid_denovo | plasmid_align | amplicon_align | compare_assembly
  Relative paths are resolved against the folder containing samples.tsv.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

from . import __version__, decide, inspect, refmap, workflows
from .runner import Runner, VERSION_CMDS

TOOLS = ["fastp", "seqkit", "unicycler", "spades.py", "flye", "medaka_consensus", "minimap2", "samtools", "bcftools"]


def cmd_check(_args) -> int:
    import shutil
    print(f"SeqGuide {__version__}")
    versions = Runner.tool_versions(TOOLS)
    missing = [t for t in TOOLS if shutil.which(t) is None]
    for t in TOOLS:
        print(f"  {'OK ' if t not in missing else 'MISSING'}  {t:18s} {versions.get(t, '')}")
    return 1 if missing else 0


def _ref_kind(path: str) -> str:
    return "annotated" if Path(path).suffix.lower() in (".gb", ".gbk", ".genbank", ".dna") else "sequence"


def cmd_batch(args) -> int:
    sheet = Path(args.samplesheet).resolve()
    base = sheet.parent
    out = Path(args.out) if args.out else base.parent / "runs" / f"batch_{time.strftime('%Y%m%d-%H%M%S')}"
    out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(sheet), delimiter="\t"))
    summary = []
    for row in rows:
        name = row["sample"].strip()
        t0 = time.time()
        reads = [str((base / p.strip()).resolve()) for p in row.get("reads", "").split(";") if p.strip()]
        ref = str((base / row["reference"].strip()).resolve()) if row.get("reference", "").strip() else None
        stype = row.get("sample_type", "plasmid").strip() or "plasmid"
        print(f"\n=== {name} ===", flush=True)
        rec_d, insp = {}, {}
        short, long_ = [], []
        params: dict = {}
        if reads:
            ds = inspect.inspect_dataset(reads)
            sample = ds.samples[0]
            insp = sample.to_dict()
            short = [f.path for f in sample.files if f.platform in ("illumina", "short")]
            long_ = [f.path for f in sample.files if f.platform in ("nanopore", "pacbio", "long")]
            size = int(row["expected_size"]) if row.get("expected_size", "").strip() else (
                sum(len(r.seq) for r in refmap.load_references(ref)) if ref else None)
            rec = decide.recommend(sample, stype, _ref_kind(ref) if ref else "none", "verify", expected_size=size)
            rec_d = rec.to_dict()
            params = dict(rec.params)
            workflow = rec.workflow
            for m in rec.messages:
                if m.level in ("warn", "block"):
                    print(f"  [{m.level}] {m.text}")
            if rec.blocked:
                print("  blocked by the decision rules; skipped")
                summary.append({"sample": name, "status": "blocked"})
                continue
        else:
            workflow = "compare_assembly"
        if row.get("assembly", "").strip():
            params["assembly"] = str((base / row["assembly"].strip()).resolve())
            workflow = "compare_assembly"
        if row.get("workflow", "").strip():
            workflow = row["workflow"].strip()
        if row.get("engine", "").strip():
            params["engine"] = row["engine"].strip()      # "builtin" forces SeqGuide's own engine
        if workflow == "plasmid_align":
            params["ploidy"] = 1
            params.setdefault("preset", "sr" if short else "map-ont")
        cfg = workflows.RunConfig(workflow=workflow, sample=name, outdir=out / name, short_reads=short, long_reads=long_,
                                  reference=ref, params=params, threads=args.threads or workflows.RunConfig.threads,
                                  inspection=insp, recommendation=rec_d,
                                  answers={"sample_type": stype, "reference": _ref_kind(ref) if ref else "none"})
        res = workflows.run(cfg, log=lambda m: print("  " + m, flush=True))
        s = {"sample": name, "workflow": workflow, "engine": res["params"].get("engine"),
             "status": "error" if res.get("error") else "ok",
             "error": res.get("error") or "", "minutes": round((time.time() - t0) / 60, 1)}
        if res.get("comparison"):
            c = res["comparison"]
            main = next((x for x in res["contigs"] if x["name"] == res["main_contig"]), {})
            s.update(ref_len=c["ref_len"], assembly_len=c["qry_len"], circular=main.get("circular"),
                     n_contigs=len(res["contigs"]), identity=c["identity"], ref_covered=c["ref_covered"],
                     n_differences=len(c["variants"]),
                     differences="; ".join(f"{v['description']} [{', '.join(v['features'])}] depth={v['depth']}"
                                           for v in c["variants"]),
                     verdict=c["verdict"]["text"])
        if res.get("amplicon"):
            a = res["amplicon"]
            vs = [v for t in a["targets"] for v in t["variants"]]
            s.update(mapped_pct=a["mapping"].get("mapped_pct"), n_differences=len(vs),
                     mean_depth=";".join(str(t["coverage"]["mean"]) for t in a["targets"]),
                     differences="; ".join(f"{v['description']} [{', '.join(v['features'])}] depth={v['depth']}" for v in vs),
                     verdict=a["verdict"]["text"])
        cov = res.get("coverage") or {}
        if cov:
            first = list(cov.values())[0]
            s.update(median_depth=first.get("median"), min_depth=first.get("min"))
        summary.append(s)
        print(f"  -> {s['status']} {s.get('n_differences', '')} difference(s) in {s['minutes']} min", flush=True)
    keys = list(dict.fromkeys(k for s in summary for k in s))
    with open(out / "batch_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    (out / "batch_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nSummary written to {out / 'batch_summary.tsv'}")
    return 0 if all(s.get("status") == "ok" for s in summary) else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="seqguide")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="check that all tools are installed")
    b = sub.add_parser("batch", help="run many samples from a samples.tsv")
    b.add_argument("samplesheet")
    b.add_argument("--out")
    b.add_argument("--threads", type=int)
    args = p.parse_args(argv)
    return {"check": cmd_check, "batch": cmd_batch}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
