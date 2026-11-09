"""Workflows: turn a recommendation into concrete, explained steps and run them.

  plasmid_denovo : trim -> (subsample) -> assemble -> (polish) -> map reads back -> compare -> report
  amplicon_align : trim -> map to targets -> call variants -> interpret -> report
  plasmid_align  : same as amplicon_align, with a circular plasmid reference (haploid)
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from . import effects, plasmid_align, refmap
from .plasmid_align import Variant
from .runner import Runner, Step
from .seqio import read_fasta, write_fasta


@dataclass
class RunConfig:
    workflow: str
    sample: str
    outdir: Path
    short_reads: list = field(default_factory=list)   # [R1, R2] or [R]
    long_reads: list = field(default_factory=list)
    reference: Optional[str] = None
    params: dict = field(default_factory=dict)
    threads: int = max(1, (os.cpu_count() or 2) - 1)
    inspection: dict = field(default_factory=dict)
    recommendation: dict = field(default_factory=dict)
    answers: dict = field(default_factory=dict)


# ====================================================================== parsers (pure Python)

def parse_depth(path: Path) -> dict:
    """samtools depth -a output -> {contig: np.array(depth)}"""
    data: dict = {}
    if not Path(path).exists():
        return {}
    with open(path) as fh:
        for line in fh:
            c, p, d = line.rstrip("\n").split("\t")[:3]
            data.setdefault(c, []).append((int(p), int(d)))
    out = {}
    for c, rows in data.items():
        n = max(p for p, _ in rows)
        arr = np.zeros(n, dtype=np.int64)
        for p, d in rows:
            arr[p - 1] = d
        out[c] = arr
    return out


def parse_flagstat(path: Path) -> dict:
    txt = Path(path).read_text() if Path(path).exists() else ""
    total = mapped = None
    for line in txt.splitlines():
        m = re.match(r"(\d+) \+ \d+ in total", line)
        if m:
            total = int(m.group(1))
        m = re.match(r"(\d+) \+ \d+ mapped \(", line)
        if m and mapped is None:
            mapped = int(m.group(1))
    pct = round(100.0 * mapped / total, 1) if total and mapped is not None else None
    return {"total": total, "mapped": mapped, "mapped_pct": pct}


def parse_vcf(path: Path) -> list:
    """Minimal VCF reader -> list of dicts (one ALT per row; run bcftools norm -m -any first)."""
    path = Path(path)
    if not path.exists():
        return []
    op = gzip.open if path.suffix == ".gz" else open
    rows = []
    with op(path, "rt") as fh:
        sample_cols = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_cols = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, _, ref, alt, qual, _, info = f[:8]
            row = {"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt.split(",")[0],
                   "qual": float(qual) if qual not in (".", "") else None, "info": info,
                   "dp": None, "af": None, "gt": None}
            if len(f) >= 10 and sample_cols:
                keys = f[8].split(":")
                vals = dict(zip(keys, f[9].split(":")))
                row["gt"] = vals.get("GT")
                if "AD" in vals and vals["AD"] not in (".", ""):
                    ad = [int(x) for x in vals["AD"].split(",") if x.isdigit()]
                    if len(ad) >= 2 and sum(ad):
                        row["af"] = round(ad[1] / sum(ad), 3)
                        row["dp"] = sum(ad)
                if row["dp"] is None and vals.get("DP", ".").isdigit():
                    row["dp"] = int(vals["DP"])
            if row["dp"] is None:
                m = re.search(r"(?:^|;)DP=(\d+)", info)
                if m:
                    row["dp"] = int(m.group(1))
            rows.append(row)
    return rows


def vcf_to_variant(row: dict) -> Variant:
    """Convert VCF alleles (with anchor base) to the plain representation used in SeqGuide."""
    ref, alt, pos0 = row["ref"], row["alt"], row["pos"] - 1
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    ref, alt, pos0 = ref[i:], alt[i:], pos0 + i
    j = 0
    while j < min(len(ref), len(alt)) and ref[-1 - j] == alt[-1 - j]:
        j += 1
    if j:
        ref, alt = ref[: len(ref) - j], alt[: len(alt) - j]
    v = Variant(pos0, ref, alt, pos0)
    v.depth = row.get("dp")
    return v


def parse_assembly_contigs(fasta: Path, info_file: Optional[Path] = None) -> list:
    """Read assembled contigs and whatever the assembler said about circularity / depth."""
    circ_info, depth_info = {}, {}
    if info_file and Path(info_file).exists():  # Flye assembly_info.txt
        for line in Path(info_file).read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) >= 4:
                circ_info[f[0]] = f[3].strip() == "Y"
                try:
                    depth_info[f[0]] = float(f[2])
                except ValueError:
                    pass
    contigs = []
    for r in read_fasta(fasta):
        desc = r.description
        circ = circ_info.get(r.id)
        if circ is None:
            m = re.search(r"circular=(true|false)", desc)
            circ = (m.group(1) == "true") if m else None
        dep = depth_info.get(r.id)
        if dep is None:
            m = re.search(r"depth=([\d.]+)x", desc) or re.search(r"_cov_([\d.]+)", r.id)  # Unicycler / SPAdes
            dep = float(m.group(1)) if m else None
        contigs.append({"name": r.id, "length": len(r.seq), "circular": circ, "depth": dep, "seq": r.seq})
    return contigs


def depth_summary(arr: np.ndarray, min_depth: int) -> dict:
    if arr is None or not len(arr):
        return {"mean": 0, "min": 0, "pct_ok": 0, "low_regions": []}
    low = arr < min_depth
    regions, start = [], None
    for i, flag in enumerate(low):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            regions.append((start + 1, i))
            start = None
    if start is not None:
        regions.append((start + 1, len(arr)))
    return {"mean": round(float(arr.mean()), 1), "min": int(arr.min()), "median": float(np.median(arr)),
            "pct_ok": round(100.0 * float((~low).mean()), 1), "low_regions": regions[:20]}


# ====================================================================== step builders

def _fastp_steps(cfg: RunConfig, state: dict) -> list:
    r = cfg.short_reads
    if not r:
        return []
    expl = ("Removes adapter sequences and low-quality read ends, so sequencing artefacts are not mistaken for mutations.")
    if len(r) == 2:
        cmd = ["fastp", "-i", r[0], "-I", r[1], "-o", "trim_R1.fastq.gz", "-O", "trim_R2.fastq.gz",
               "--detect_adapter_for_pe", "-w", str(min(cfg.threads, 16)), "-j", "fastp.json", "-h", "fastp.html"]
        state["short"] = ["trim_R1.fastq.gz", "trim_R2.fastq.gz"]
    else:
        cmd = ["fastp", "-i", r[0], "-o", "trim.fastq.gz", "-w", str(min(cfg.threads, 16)), "-j", "fastp.json", "-h", "fastp.html"]
        state["short"] = ["trim.fastq.gz"]
    return [Step("Quality trimming", expl, [cmd], tools=["fastp"])]


def _subsample_steps(cfg: RunConfig, state: dict) -> list:
    steps = []
    target = cfg.params.get("subsample_depth", 150)
    for kind in ("short", "long"):
        d = cfg.params.get(f"depth_{kind}")
        files = state.get(kind) or []
        if not d or not files or d <= 3 * target:
            continue
        frac = round(target / d, 4)
        new = []
        cmds = []
        for f in files:
            out = f"sub_{Path(f).name}"
            if not out.endswith(".gz"):
                out += ".gz"
            cmds.append(["seqkit", "sample", "-p", str(frac), "-s", "11", "-o", out, f])
            new.append(out)
        for c in cmds:
            steps.append(Step(f"Subsampling {kind} reads" + (f" ({Path(c[-1]).name})" if len(cmds) > 1 else ""),
                              f"Your data has about {d:,.0f}× depth. Assemblers work best around {target}×, so a random "
                              f"{frac:.1%} of reads is kept (same reads for R1 and R2).", [c], tools=["seqkit"]))
        state[kind] = new
    return steps


def plasmid_denovo_steps(cfg: RunConfig, state: dict) -> list:
    steps = _fastp_steps(cfg, state)
    if cfg.long_reads:
        state["long"] = ["long_filtered.fastq.gz"]
        steps.append(Step("Long-read filtering",
                          "Removes reads shorter than 1,000 bp (broken fragments), keeping mostly full-length plasmid copies.",
                          [["seqkit", "seq", "-m", "1000", "-o", "long_filtered.fastq.gz", cfg.long_reads[0]]], tools=["seqkit"]))
    state["support_short"] = list(state.get("short") or [])   # all trimmed reads, for read support
    state["support_long"] = list(state.get("long") or [])
    steps += _subsample_steps(cfg, state)
    asm = cfg.params.get("assembler", "unicycler")
    t = str(cfg.threads)
    if asm in ("unicycler", "unicycler-hybrid"):
        cmd = ["unicycler", "-o", "assembly", "-t", t]
        s = state.get("short", [])
        if len(s) == 2:
            cmd += ["-1", s[0], "-2", s[1]]
        elif s:
            cmd += ["-s", s[0]]
        if asm == "unicycler-hybrid" and state.get("long"):
            cmd += ["-l", state["long"][0]]
        steps.append(Step("De novo assembly (Unicycler)",
                          "Builds the plasmid sequence from the reads alone and checks whether its ends join into a circle. "
                          "Unicycler is designed for small circular sequences such as plasmids.",
                          [cmd], tools=["unicycler", "spades.py"]))
        state["assembly_fasta"] = "assembly/assembly.fasta"
        state["assembly_info"] = None
    else:
        cmd = ["flye", "--nano-hq", state["long"][0], "--out-dir", "flye", "--threads", t]
        size = cfg.params.get("expected_size")
        if size:
            cmd += ["--genome-size", str(size), "--asm-coverage", "50"]
        steps.append(Step("De novo assembly (Flye)",
                          "Builds the plasmid sequence from long reads alone. Long reads often span the whole plasmid, "
                          "so the structure comes out complete.", [cmd], tools=["flye"]))
        steps.append(Step("Polishing (Medaka)",
                          "Corrects small errors in the Nanopore assembly using a model of how the sequencer makes mistakes.",
                          [["medaka_consensus", "-i", state["long"][0], "-d", "flye/assembly.fasta", "-o", "medaka", "-t", t]],
                          tools=["medaka_consensus"], optional=True))
        state["assembly_fasta"] = "medaka/consensus.fasta"
        state["assembly_fallback"] = "flye/assembly.fasta"
        state["assembly_info"] = "flye/assembly_info.txt"

    steps.append(Step("Collect assembly", "Reads the assembled contigs and checks size and circularity.",
                      func=lambda: _collect_assembly(cfg, state)))
    reads = state.get("support_short") or state.get("support_long")
    preset = "sr" if state.get("support_short") else "map-ont"
    steps.append(Step("Map reads back to assembly",
                      "Places all reads (not only the subsample used for assembly) on the finished assembly to confirm "
                      "that every position is supported by many reads.",
                      [["minimap2", "-ax", preset, "-t", t, "final_assembly.fasta", *reads],
                       ["samtools", "sort", "-o", "reads_on_assembly.bam", "-"]], tools=["minimap2", "samtools"]))
    steps.append(Step("Index alignment", "Technical step needed to read the alignment quickly.",
                      [["samtools", "index", "reads_on_assembly.bam"]], tools=["samtools"]))
    steps.append(Step("Measure depth", "Counts how many reads cover each position of the assembly.",
                      [["samtools", "depth", "-a", "reads_on_assembly.bam"]], stdout_to="assembly_depth.txt", tools=["samtools"]))
    if cfg.reference:
        steps.append(Step("Compare to reference",
                          "Lines up the assembled plasmid with your expected map, lists every difference, "
                          "and explains which feature it hits.", func=lambda: _compare_plasmid(cfg, state)))
    else:
        steps.append(Step("Summarise assembly", "Summarises the assembly (no reference given, so no comparison).",
                          func=lambda: _assembly_only(cfg, state)))
    return steps


def align_steps(cfg: RunConfig, state: dict) -> list:
    steps = _fastp_steps(cfg, state)
    if cfg.long_reads and not cfg.short_reads:
        state["long"] = list(cfg.long_reads)
    reads = state.get("short") or state.get("long")
    preset = cfg.params.get("preset") or ("sr" if state.get("short") else "map-ont")
    t = str(cfg.threads)
    ploidy = cfg.params.get("ploidy", 1 if cfg.workflow == "plasmid_align" else 2)
    steps.append(Step("Prepare reference", "Converts your reference into the plain format the aligner needs.",
                      func=lambda: _prepare_reference(cfg, state)))
    steps.append(Step("Index reference", "Technical step so tools can look up reference positions quickly.",
                      [["samtools", "faidx", "ref.fa"]], tools=["samtools"]))
    steps.append(Step("Align reads to reference",
                      "Places every read at its matching position on your target sequence(s).",
                      [["minimap2", "-ax", preset, "-t", t, "ref.fa", *reads],
                       ["samtools", "sort", "-o", "aligned.bam", "-"]], tools=["minimap2", "samtools"]))
    steps.append(Step("Index alignment", "Technical step needed to read the alignment quickly.",
                      [["samtools", "index", "aligned.bam"]], tools=["samtools"]))
    steps.append(Step("Mapping statistics", "Counts how many reads matched your targets.",
                      [["samtools", "flagstat", "aligned.bam"]], stdout_to="flagstat.txt", tools=["samtools"]))
    steps.append(Step("Measure depth", "Counts how many reads cover each position of each target.",
                      [["samtools", "depth", "-a", "aligned.bam"]], stdout_to="depth.txt", tools=["samtools"]))
    call = ["bcftools", "call", "-mv", "-Ou"]
    if ploidy == 1:
        call += ["--ploidy", "1"]
    mpileup = ["bcftools", "mpileup", "-f", "ref.fa", "-d", "100000", "-a", "FORMAT/AD,FORMAT/DP", "-Ou", "aligned.bam"]
    if preset == "map-ont":
        mpileup[2:2] = ["-X", "ont"]  # error model tuned for Nanopore reads
    steps.append(Step("Find variants",
                      "Looks at every position and decides, using all reads together, whether your sample differs from the reference. "
                      + ("Assumes one copy per cell (plasmid / haploid)." if ploidy == 1 else
                         "Assumes two copies per cell (diploid genomic DNA), so heterozygous changes (~50% of reads) are reported."),
                      [mpileup,
                       call,
                       ["bcftools", "norm", "-f", "ref.fa", "-m", "-any", "-Ov", "-"]],
                      stdout_to="variants.vcf", tools=["bcftools"]))
    steps.append(Step("Interpret variants", "Explains each variant: where it is, which feature it hits, and how well it is supported.",
                      func=lambda: _interpret_alignment(cfg, state)))
    return steps


def compare_assembly_steps(cfg: RunConfig, state: dict) -> list:
    """For users who already have an assembly (own run, or a sequencing service's result)."""
    state["assembly_fasta"] = cfg.params["assembly"]
    state["assembly_info"] = None
    steps = [Step("Load assembly", "Reads your existing assembly and checks its contigs.", func=lambda: _collect_assembly(cfg, state))]
    if cfg.short_reads:
        steps.append(Step("Compare to reference and check reads",
                          "Lines up the assembled plasmid with your expected map and lists every difference. Your reads are then "
                          "placed on the map to confirm each difference and to measure depth along the plasmid.",
                          func=lambda: _compare_plasmid(cfg, state)))
    else:
        steps.append(Step("Compare to reference", "Lines up the assembled plasmid with your expected map, lists every difference, "
                                                  "and explains which feature it hits.", func=lambda: _compare_plasmid(cfg, state)))
    return steps


def build_steps(cfg: RunConfig, state: dict) -> list:
    if cfg.workflow == "compare_assembly":
        return compare_assembly_steps(cfg, state)
    if cfg.workflow == "plasmid_denovo":
        return plasmid_denovo_steps(cfg, state)
    if cfg.workflow in ("amplicon_align", "plasmid_align"):
        return align_steps(cfg, state)
    raise ValueError(f"Unknown workflow {cfg.workflow}")


# ====================================================================== python steps

def _wd(cfg: RunConfig) -> Path:
    return Path(cfg.outdir)


def _collect_assembly(cfg: RunConfig, state: dict) -> None:
    wd = _wd(cfg)
    fa = Path(state["assembly_fasta"])
    fa = fa if fa.is_absolute() else wd / fa
    if not fa.exists() and state.get("assembly_fallback"):
        fa = wd / state["assembly_fallback"]
        state.setdefault("warnings", []).append("Polishing did not run, so the unpolished assembly is used. "
                                                "Expect more small errors in runs like AAAAAA.")
    if not fa.exists() or fa.stat().st_size == 0:
        raise RuntimeError("The assembler produced no sequence. Usually this means too few reads or very low quality.")
    info = wd / state["assembly_info"] if state.get("assembly_info") else None
    contigs = parse_assembly_contigs(fa, info)
    contigs.sort(key=lambda c: -c["length"])
    write_fasta([(c["name"], c["seq"]) for c in contigs], wd / "final_assembly.fasta")
    state["contigs"] = contigs


def _attach_depth(variants: list, depth: Optional[np.ndarray], comp: plasmid_align.Comparison, contig_len: int) -> None:
    if depth is None:
        return
    for v in variants:
        q = (comp.rotation + v.qry_pos) % contig_len
        if comp.strand == "-":
            q = contig_len - 1 - q
        if 0 <= q < len(depth):
            v.depth = int(depth[q])


def _compare_plasmid(cfg: RunConfig, state: dict) -> None:
    wd = _wd(cfg)
    refs = refmap.load_references(cfg.reference, assume_circular=True)
    ref = refs[0]
    if len(refs) > 1:
        state.setdefault("warnings", []).append(f"The reference file has {len(refs)} sequences; only the first "
                                                f"({ref.name}) is used for comparison.")
    depth = parse_depth(wd / "assembly_depth.txt")
    contigs = [c for c in state["contigs"] if c["length"] >= 500]
    best = None
    for c in contigs[:10]:
        comp = plasmid_align.compare(c["seq"], ref.seq, circular=bool(c["circular"]) or c["circular"] is None,
                                     qry_name=c["name"], ref_name=ref.name)
        if best is None or comp.ref_covered > best[1].ref_covered:
            best = (c, comp)
    if best is None:
        raise RuntimeError("No assembled contig longer than 500 bp to compare.")
    c, comp = best
    platform = "nanopore" if (cfg.long_reads and not cfg.short_reads) else "illumina"
    effects.annotate(ref, comp.variants, platform=platform)
    if comp.trimmed_overlap:
        c["circular"] = True
        c["note"] = f"{comp.trimmed_overlap} bp end overlap removed (circular)"
    _attach_depth(comp.variants, depth.get(c["name"]), comp, comp.qry_len)
    min_d = 20 if platform == "nanopore" else 30
    rc = None
    if not depth and cfg.short_reads:
        # No read alignment from external tools (e.g. "compare an existing assembly"): use the built-in read check
        from .readcheck import check_reads
        rc = check_reads(cfg.short_reads[0], ref.seq)
        for v in comp.variants:
            if v.ref_pos < len(rc.depth):
                v.depth = int(rc.depth[v.ref_pos])
            if v.kind == "SNP" and v.depth:
                n_alt = rc.alt_counts.get(v.ref_pos, {}).get(v.alt, 0)
                good = max(1, v.depth - int(rc.lowq[v.ref_pos]))  # reads with a high-quality base here
                v.flags.append(f"Confirmed by reads: {n_alt} of {good} high-quality reads ({n_alt / good:.0%}) carry it.")
        known = {v.ref_pos for v in comp.variants}
        extra = [x for x in rc.variants() if x["pos"] - 1 not in known]
        state["read_check"] = {"reads": rc.n_reads, "pct_mapped": rc.pct_mapped, "extra_read_variants": extra}
        if extra:
            state.setdefault("warnings", []).append(
                f"{len(extra)} position(s) where 5% or more of reads differ from the reference but the assembly does not, "
                "e.g. " + ", ".join(f"{x['pos']} ({x['af']:.0%} {x['alt']})" for x in extra[:5]) +
                ". This can mean a mixed population (two clones) or a sequencing artefact.")
    for v in comp.variants:
        if v.depth is not None and v.depth < min_d:
            v.flags.append(f"Only {v.depth} reads cover this position; treat it as uncertain.")
    # rotated assembly is directly comparable with the reference map
    write_fasta([(f"{cfg.sample}_assembly_rotated_to_{ref.name}", comp.rotated_seq)], wd / "assembly_matched_to_reference.fasta")
    _write_variants_csv(comp.variants, wd / "differences.csv")
    for oc in state["contigs"]:
        oc["matched_reference"] = oc["name"] == c["name"]
        d = depth.get(oc["name"])
        oc["mean_depth"] = round(float(d.mean()), 1) if d is not None and len(d) else oc.get("depth")
    state["main_contig"] = c["name"]
    state["comparison"] = comp.to_dict()
    state["comparison"]["verdict"] = effects.verdict(comp.variants, comp.ref_covered, comp.identity, comp.notes)
    state["reference"] = _ref_summary(ref)
    dep = rc.depth if rc is not None else depth.get(c["name"])
    state["coverage"] = {c["name"]: depth_summary(dep, min_d)} if dep is not None else None
    if rc is not None:
        state["coverage"][c["name"]]["source"] = "reads aligned to the reference (built-in read check)"
    state["plot_inputs"] = {"kind": "plasmid", "ref": ref, "variants": comp.variants,
                            "depth": rc.depth if rc is not None else
                            _depth_in_ref_coords(depth.get(c["name"]), comp, comp.qry_len, len(ref.seq))}


def _depth_in_ref_coords(depth, comp, contig_len: int, ref_len: int):
    """Project assembly depth onto reference coordinates via the anchor map (approximate, for plotting)."""
    if depth is None or not comp.ref_to_qry:
        return None
    rpos = np.array([0] + [a[0] for a in comp.ref_to_qry] + [ref_len - 1])
    qpos = np.array([0] + [a[1] for a in comp.ref_to_qry] + [contig_len - 1])
    x = np.arange(ref_len)
    q = np.clip(np.interp(x, rpos, qpos).astype(int), 0, contig_len - 1)
    q = (comp.rotation + q) % contig_len
    if comp.strand == "-":
        q = contig_len - 1 - q
    q = np.clip(q, 0, len(depth) - 1)
    return depth[q]


def _assembly_only(cfg: RunConfig, state: dict) -> None:
    wd = _wd(cfg)
    depth = parse_depth(wd / "assembly_depth.txt")
    for c in state["contigs"]:
        d = depth.get(c["name"])
        c["mean_depth"] = round(float(d.mean()), 1) if d is not None and len(d) else c.get("depth")
    main = state["contigs"][0]
    state["main_contig"] = main["name"]
    state["coverage"] = {main["name"]: depth_summary(depth.get(main["name"]), 30)}
    state["plot_inputs"] = {"kind": "assembly", "depth": depth.get(main["name"]), "name": main["name"]}


def _prepare_reference(cfg: RunConfig, state: dict) -> None:
    refs = refmap.load_references(cfg.reference)
    write_fasta([(r.name, r.seq) for r in refs], _wd(cfg) / "ref.fa")
    state["refs"] = refs


def _interpret_alignment(cfg: RunConfig, state: dict) -> None:
    wd = _wd(cfg)
    refs = {r.name: r for r in state["refs"]}
    rows = parse_vcf(wd / "variants.vcf")
    depth = parse_depth(wd / "depth.txt")
    min_d = cfg.params.get("min_depth", 100 if cfg.workflow == "amplicon_align" else 30)
    platform = "nanopore" if (cfg.long_reads and not cfg.short_reads) else "illumina"
    ploidy = cfg.params.get("ploidy", 1 if cfg.workflow == "plasmid_align" else 2)
    by_target: dict = {name: [] for name in refs}
    for row in rows:
        ref = refs.get(row["chrom"])
        if not ref:
            continue
        v = vcf_to_variant(row)
        v.qry_pos = v.ref_pos
        by_target[row["chrom"]].append((v, row))
    targets = []
    all_vars = []
    for name, ref in refs.items():
        items = by_target[name]
        vs = [v for v, _ in items]
        for v in vs:
            v.homopolymer = plasmid_align.homopolymer_len(ref.seq, v)
        effects.annotate(ref, vs, platform=platform)
        for v, row in items:
            af = row.get("af")
            gt = row.get("gt") or ""
            if af is not None:
                v.flags.append(f"Seen in {af:.0%} of reads.")
                if ploidy == 2 and gt in ("0/1", "0|1", "1|0"):
                    v.flags.append("Heterozygous: one of the two copies carries it.")
                elif ploidy == 1 and 0.2 < af < 0.8:
                    v.flags.append("Not in all reads: possibly a mixed population (two clones). Consider re-streaking.")
            if v.depth is not None and v.depth < min_d:
                v.flags.append(f"Low depth ({v.depth} reads; ≥ {min_d} recommended): uncertain.")
            if row.get("qual") is not None and row["qual"] < 20:
                v.flags.append(f"Low confidence score (QUAL {row['qual']:.0f}).")
        cov = depth_summary(depth.get(name), min_d)
        status = "ok" if cov["pct_ok"] >= 95 else ("warn" if cov["pct_ok"] >= 50 else "fail")
        targets.append({"name": name, "length": len(ref.seq), "annotated": ref.annotated, "coverage": cov,
                        "status": status, "variants": [v.to_dict() for v in vs]})
        all_vars.extend(vs)
    _write_variants_csv(all_vars, wd / "variants_explained.csv", with_target=[t["name"] for t in targets for _ in t["variants"]])
    state["amplicon"] = {"targets": targets, "mapping": parse_flagstat(wd / "flagstat.txt"), "min_depth": min_d,
                         "ploidy": ploidy}
    state["plot_inputs"] = {"kind": "targets", "depth": depth, "refs": state["refs"], "min_depth": min_d,
                            "variants": {t["name"]: [v for v, _ in by_target[t["name"]]] for t in targets}}
    n_high = sum(1 for v in all_vars if v.severity == "high")
    failed = [t["name"] for t in targets if t["status"] == "fail"]
    if failed:
        verdict = {"level": "warn", "text": f"{len(failed)} target(s) have too little coverage to analyse: {', '.join(failed[:5])}."}
    elif n_high:
        verdict = {"level": "block", "text": f"{n_high} variant(s) likely disrupt a protein or feature."}
    elif all_vars:
        verdict = {"level": "warn", "text": f"{len(all_vars)} variant(s) found. Review the effects below."}
    else:
        verdict = {"level": "ok", "text": "No differences from the reference in well-covered regions."}
    state["amplicon"]["verdict"] = verdict


def _ref_summary(ref: refmap.Reference) -> dict:
    return {"name": ref.name, "length": len(ref.seq), "circular": ref.circular, "annotated": ref.annotated,
            "features": [{"label": f.label, "type": f.type, "start": f.start + 1, "end": f.end, "strand": f.strand}
                         for f in ref.features]}


def _write_variants_csv(variants: list, path: Path, with_target: Optional[list] = None) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        head = (["target"] if with_target else []) + ["position", "type", "change", "reference", "sample", "feature",
                                                      "effect", "severity", "read_depth", "notes"]
        w.writerow(head)
        for i, v in enumerate(variants):
            row = ([with_target[i]] if with_target else []) + [
                v.ref_pos + 1, v.kind, v.describe(), v.ref[:50], v.alt[:50], "; ".join(v.features), v.effect,
                v.severity, v.depth if v.depth is not None else "", " ".join(v.flags)]
            w.writerow(row)




# ====================================================================== built-in engine (no external tools)
# Used automatically when Unicycler/minimap2/... are not installed, e.g. in the browser version.

BUILTIN_WORKFLOWS = ("plasmid_denovo", "compare_assembly", "amplicon_align", "plasmid_align")


def _read_records(paths: list):
    from .seqio import read_fastq
    for p in paths:
        for r in read_fastq(p):
            yield r.seq, r.qual


def builtin_steps(cfg: RunConfig, state: dict) -> list:
    if cfg.workflow == "compare_assembly":
        return compare_assembly_steps(cfg, state)
    if cfg.workflow == "plasmid_denovo":
        long_mode = bool(cfg.long_reads and not cfg.short_reads)
        steps = [Step("De novo assembly (built-in)",
                      ("Builds the plasmid from long reads: a complete read is taken as a draft and corrected with all "
                       "other full-length reads." if long_mode else
                       "Cleans the reads (adapters, poly-G tails, low-quality ends), then builds the plasmid sequence from "
                       "the reads alone with SeqGuide's own assembler, resolves repeats and closes the circle. Validated "
                       "base-for-base against Unicycler."),
                      func=lambda: _builtin_assemble(cfg, state))]
        if cfg.reference:
            steps.append(Step("Compare to reference and check reads",
                              "Lines up the assembled plasmid with your expected map, lists every difference and checks "
                              "it against the reads.", func=lambda: _compare_plasmid(cfg, state)))
        else:
            steps.append(Step("Check the assembly against the reads",
                              "Places the reads on the finished assembly to confirm that every position is supported.",
                              func=lambda: _builtin_assembly_only(cfg, state)))
        return steps
    return [Step("Prepare reference", "Reads your target sequence(s).", func=lambda: _prepare_reference(cfg, state)),
            Step("Align reads and find variants (built-in)",
                 "Places every read on your target sequence(s), measures depth at every position and finds positions "
                 "where your sample differs, with the share of reads that show each change.",
                 func=lambda: _builtin_align(cfg, state))]


def _builtin_assemble(cfg: RunConfig, state: dict) -> None:
    from . import pyassembler
    size = cfg.params.get("expected_size")
    if cfg.long_reads and not cfg.short_reads:
        res = pyassembler.assemble_long(_read_records(cfg.long_reads), expected_size=size, log=lambda m: None)
    else:
        res = pyassembler.assemble(_read_records(cfg.short_reads), expected_size=size, log=lambda m: None)
    if not res or not res.sequence:
        raise RuntimeError("The assembler could not build a sequence. " + " ".join(res.notes if res else []))
    state["contigs"] = [{"name": c["name"], "length": len(c["seq"]), "circular": c["circular"], "depth": c["coverage"],
                         "seq": c["seq"]} for c in res.contigs]
    state["assembly_stats"] = res.stats
    state.setdefault("warnings", []).extend(res.notes)
    for mv in res.minor_variants:
        state["warnings"].append(f"About {mv['fraction']:.0%} of reads carry a different version of one region: "
                                 "possibly a mixed population (two clones). Consider re-streaking a single colony.")
    write_fasta([(c["name"], c["seq"]) for c in state["contigs"]], _wd(cfg) / "final_assembly.fasta")


def _builtin_assembly_only(cfg: RunConfig, state: dict) -> None:
    main = state["contigs"][0]
    state["main_contig"] = main["name"]
    if cfg.short_reads:
        from .readcheck import map_and_pileup
        rc = map_and_pileup(_read_records(cfg.short_reads), main["seq"], circular=bool(main["circular"]))
        main["mean_depth"] = round(float(rc.depth.mean()), 1)
        state["coverage"] = {main["name"]: depth_summary(rc.depth, 30)}
        state["plot_inputs"] = {"kind": "assembly", "depth": rc.depth, "name": main["name"]}


def _builtin_align(cfg: RunConfig, state: dict) -> None:
    from .readcheck import map_and_pileup
    min_d = cfg.params.get("min_depth", 100 if cfg.workflow == "amplicon_align" else 30)
    ploidy = cfg.params.get("ploidy", 1 if cfg.workflow == "plasmid_align" else 2)
    reads = cfg.short_reads or cfg.long_reads
    targets, all_vars, plot_vars, depths = [], [], {}, {}
    total = mapped_any = 0
    for ref in state["refs"]:
        circ = cfg.workflow == "plasmid_align" or ref.circular
        rc = map_and_pileup(_read_records(reads), ref.seq, circular=circ)
        total = max(total, rc.n_reads)
        mapped_any += rc.n_mapped
        depths[ref.name] = rc.depth
        vs = []
        for x in rc.variants(min_af=0.15, min_alt=3):
            pos = x["pos"] - 1
            v = Variant(pos, ref.seq[pos], x["alt"], pos)
            v.depth = x["depth"]
            v.homopolymer = 0
            af = x["af"]
            v.flags.append(f"Seen in {af:.0%} of reads.")
            if ploidy == 2 and af < 0.8:
                v.flags.append("Heterozygous: one of the two copies carries it." if af >= 0.3 else
                               "Seen in a minority of reads: possible mosaic, contamination or sequencing artefact.")
            elif ploidy == 1 and af < 0.8:
                v.flags.append("Not in all reads: possibly a mixed population (two clones). Consider re-streaking.")
            if x["depth"] < min_d:
                v.flags.append(f"Low depth ({x['depth']} reads; >= {min_d} recommended): uncertain.")
            vs.append(v)
        effects.annotate(ref, vs)
        for pos, n in sorted(getattr(rc, "clips", {}).items()):
            d = int(rc.depth[pos]) if pos < len(rc.depth) else 0
            if n >= 10 and d and n / (n + d) >= 0.2 and 50 < pos < len(ref.seq) - 50:
                state.setdefault("warnings", []).append(
                    f"{ref.name}: {n} reads stop matching at position {pos + 1}. This usually means an insertion or "
                    "deletion there, which the built-in engine cannot spell out. Check it with the full (Docker) version.")
        cov = depth_summary(rc.depth, min_d)
        status = "ok" if cov["pct_ok"] >= 95 else ("warn" if cov["pct_ok"] >= 50 else "fail")
        targets.append({"name": ref.name, "length": len(ref.seq), "annotated": ref.annotated, "coverage": cov,
                        "status": status, "variants": [v.to_dict() for v in vs]})
        all_vars.extend(vs)
        plot_vars[ref.name] = vs
    pct = round(100.0 * min(mapped_any, total) / total, 1) if total else None
    _write_variants_csv(all_vars, _wd(cfg) / "variants_explained.csv",
                        with_target=[t["name"] for t in targets for _ in t["variants"]])
    state["amplicon"] = {"targets": targets, "mapping": {"total": total, "mapped": mapped_any, "mapped_pct": pct},
                         "min_depth": min_d, "ploidy": ploidy}
    state["plot_inputs"] = {"kind": "targets", "depth": depths, "refs": state["refs"], "min_depth": min_d,
                            "variants": plot_vars}
    n_high = sum(1 for v in all_vars if v.severity == "high")
    failed = [t["name"] for t in targets if t["status"] == "fail"]
    if failed:
        verdict = {"level": "warn", "text": f"{len(failed)} target(s) have too little coverage to analyse: {', '.join(failed[:5])}."}
    elif n_high:
        verdict = {"level": "block", "text": f"{n_high} variant(s) likely disrupt a protein or feature."}
    elif all_vars:
        verdict = {"level": "warn", "text": f"{len(all_vars)} variant(s) found. Review the effects below."}
    else:
        verdict = {"level": "ok", "text": "No differences from the reference in well-covered regions."}
    state["amplicon"]["verdict"] = verdict


# ====================================================================== entry point

class Job:
    """One analysis, run step by step (the GUI runs one step per screen refresh so progress stays visible)."""

    def __init__(self, cfg: RunConfig, log: Optional[Callable[[str], None]] = None):
        cfg.outdir = Path(cfg.outdir)
        cfg.outdir.mkdir(parents=True, exist_ok=True)
        self.cfg = cfg
        self.state: dict = {"warnings": []}
        engine = cfg.params.get("engine", "auto")
        steps = build_steps(cfg, self.state) if engine != "builtin" else []
        self.runner = Runner(cfg.outdir, log=log)
        missing = self.runner.missing_tools(steps) if steps else ["(built-in requested)"]
        if missing:
            if cfg.workflow not in BUILTIN_WORKFLOWS or engine == "tools":
                raise RuntimeError("These tools are not installed: " + ", ".join(missing) +
                                   ". Start SeqGuide with the Docker launcher, which includes all of them.")
            self.state = {"warnings": []}
            steps = builtin_steps(cfg, self.state)
            cfg.params["engine"] = "builtin"
        else:
            cfg.params["engine"] = "tools"
        self.steps = steps
        self.index = 0
        self.error: Optional[str] = None
        self.started = time.strftime("%Y-%m-%d %H:%M %Z")

    @property
    def done(self) -> bool:
        return self.error is not None or self.index >= len(self.steps)

    @property
    def current(self) -> Optional[Step]:
        return self.steps[self.index] if self.index < len(self.steps) else None

    def set_log(self, log: Optional[Callable[[str], None]]) -> None:
        self.runner._log_cb = log

    def step(self) -> None:
        if self.done:
            return
        try:
            self.runner.run_one(self.steps[self.index], self.steps)
        except Exception as exc:  # noqa: BLE001
            self.error = str(exc)
        self.index += 1

    def finish(self) -> dict:
        from . import report  # local import keeps matplotlib out of light-weight imports
        cfg, state, runner = self.cfg, self.state, self.runner
        runner.save_provenance(self.steps)
        results = {
            "sample": cfg.sample, "workflow": cfg.workflow, "started": self.started,
            "finished": time.strftime("%Y-%m-%d %H:%M %Z"),
            "inspection": cfg.inspection, "recommendation": cfg.recommendation, "answers": cfg.answers,
            "params": {k: v for k, v in cfg.params.items()}, "steps": runner.records,
            "warnings": state.get("warnings", []) + runner.warnings, "error": self.error,
            "contigs": [{k: v for k, v in c.items() if k != "seq"} for c in state.get("contigs", [])],
            "main_contig": state.get("main_contig"), "comparison": state.get("comparison"),
            "reference": state.get("reference"), "coverage": state.get("coverage"), "amplicon": state.get("amplicon"),
            "read_check": state.get("read_check"), "assembly_stats": state.get("assembly_stats"),
        }
        prov = cfg.outdir / "provenance.json"
        if prov.exists():
            results["tool_versions"] = json.loads(prov.read_text()).get("tool_versions", {})
        (cfg.outdir / "results.json").write_text(json.dumps(results, indent=2, default=str))
        report.write_report(results, state.get("plot_inputs"), cfg.outdir / "report.html")
        return results


def run(cfg: RunConfig, log: Optional[Callable[[str], None]] = None,
        on_step: Optional[Callable] = None) -> dict:
    job = Job(cfg, log=log)
    while not job.done:
        if on_step:
            on_step(job.index, len(job.steps), job.current)
        job.step()
    return job.finish()
