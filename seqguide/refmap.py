"""Load reference sequences: FASTA, GenBank (.gb/.gbk) and SnapGene (.dna).

GenBank is parsed by a small built-in parser, so annotated references work
even without Biopython. SnapGene files need Biopython (installed in the
Docker image).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from .seqio import read_fasta, sniff_format

LABEL_KEYS = ("label", "gene", "product", "note", "standard_name", "locus_tag")
SKIP_TYPES = {"source"}


@dataclass
class Feature:
    type: str
    label: str
    parts: list               # list of (start0, end, strand), in biological (5'->3') order
    qualifiers: dict = field(default_factory=dict)

    @property
    def strand(self) -> int:
        return self.parts[0][2] if self.parts else 1

    @property
    def start(self) -> int:
        return min(p[0] for p in self.parts)

    @property
    def end(self) -> int:
        return max(p[1] for p in self.parts)

    def overlaps(self, s: int, e: int) -> bool:
        """Half-open interval overlap; a zero-length interval (insertion point) counts if inside."""
        for ps, pe, _ in self.parts:
            if s == e:
                if ps < s < pe:
                    return True
            elif s < pe and ps < e:
                return True
        return False

    def length(self) -> int:
        return sum(pe - ps for ps, pe, _ in self.parts)


@dataclass
class Reference:
    name: str
    seq: str
    circular: bool
    features: list = field(default_factory=list)
    source_format: str = "fasta"

    @property
    def annotated(self) -> bool:
        return bool(self.features)


# ------------------------------------------------------------ location parsing

def parse_location(loc: str) -> list:
    """Parse a GenBank location string into [(start0, end, strand), ...] in biological order."""
    loc = re.sub(r"\s+", "", loc)
    return _parse_loc(loc, 1)


def _split_top(s: str) -> list:
    parts, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts


def _parse_loc(s: str, strand: int) -> list:
    if s.startswith("complement(") and s.endswith(")"):
        inner = _parse_loc(s[len("complement("):-1], 1)
        return [(a, b, -st) for a, b, st in reversed(inner)]
    for kw in ("join(", "order("):
        if s.startswith(kw) and s.endswith(")"):
            out = []
            for p in _split_top(s[len(kw):-1]):
                out.extend(_parse_loc(p, 1))
            return out
    if ":" in s:  # reference to another record; not supported
        return []
    s = s.replace("<", "").replace(">", "")
    if ".." in s:
        a, b = s.split("..", 1)
        return [(int(a) - 1, int(b), strand)]
    if "^" in s:
        a = s.split("^", 1)[0]
        return [(int(a), int(a), strand)]
    if s.isdigit():
        return [(int(s) - 1, int(s), strand)]
    return []


def _label(ftype: str, quals: dict) -> str:
    for k in LABEL_KEYS:
        if k in quals and quals[k]:
            return quals[k][:60]
    return ftype


# ------------------------------------------------------------ GenBank parser

def parse_genbank(path: str | Path) -> list:
    refs = []
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    for chunk in re.split(r"^//\s*$", text, flags=re.M):
        if "LOCUS" not in chunk:
            continue
        refs.append(_parse_gb_record(chunk))
    return refs


def _parse_gb_record(chunk: str) -> Reference:
    lines = chunk.splitlines()
    name, circular = "reference", False
    features_raw = []
    seq_parts = []
    section = None
    for line in lines:
        if line.startswith("LOCUS"):
            toks = line.split()
            name = toks[1] if len(toks) > 1 else name
            circular = "circular" in line.lower()
            section = None
            continue
        if line.startswith("FEATURES"):
            section = "features"
            continue
        if line.startswith("ORIGIN"):
            section = "origin"
            continue
        if section == "features":
            if line[:1].strip():  # a new top-level keyword ends the feature table
                section = None
                continue
            key = line[5:21].strip()
            body = line[21:].rstrip()
            if key:
                features_raw.append([key, body, []])
            elif features_raw:
                if body.startswith("/"):
                    features_raw[-1][2].append(body)
                elif features_raw[-1][2]:
                    features_raw[-1][2][-1] += (" " if not features_raw[-1][2][-1].endswith("-") else "") + body
                else:
                    features_raw[-1][1] += body
        elif section == "origin":
            seq_parts.append(re.sub(r"[^A-Za-z]", "", line))
    seq = "".join(seq_parts).upper()
    feats = []
    for ftype, loc, qlines in features_raw:
        if ftype in SKIP_TYPES:
            continue
        quals = {}
        for q in qlines:
            m = re.match(r'/([^=]+)(?:=(.*))?', q)
            if not m:
                continue
            k, v = m.group(1), (m.group(2) or "")
            v = v.strip().strip('"')
            if k == "translation":
                v = v.replace(" ", "")
            quals.setdefault(k, v)
        try:
            parts = parse_location(loc)
        except ValueError:
            parts = []
        if parts:
            feats.append(Feature(ftype, _label(ftype, quals), parts, quals))
    return Reference(name=name, seq=seq, circular=circular, features=feats, source_format="genbank")


# ------------------------------------------------------------ Biopython path

def _from_biopython(path: Path, fmt: str) -> list:
    try:
        from Bio import SeqIO  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Reading SnapGene (.dna) files needs Biopython. Export the map as GenBank (.gb) from SnapGene, "
                           "or run SeqGuide via Docker where Biopython is installed.") from exc
    refs = []
    for rec in SeqIO.parse(str(path), fmt):
        feats = []
        for f in rec.features:
            if f.type in SKIP_TYPES:
                continue
            quals = {k: (v[0] if isinstance(v, list) and v else str(v)) for k, v in f.qualifiers.items()}
            parts = [(int(p.start), int(p.end), int(p.strand or 1)) for p in f.location.parts]
            feats.append(Feature(f.type, _label(f.type, quals), parts, quals))
        topo = rec.annotations.get("topology", "")
        refs.append(Reference(rec.id or rec.name, str(rec.seq).upper(), topo == "circular", feats, fmt))
    return refs


# ------------------------------------------------------------ public API

def load_references(path: str | Path, assume_circular: bool | None = None) -> list:
    path = Path(path)
    fmt = sniff_format(path)
    if fmt == "genbank":
        refs = parse_genbank(path)
    elif fmt == "snapgene":
        refs = _from_biopython(path, "snapgene")
    elif fmt == "fasta":
        refs = [Reference(r.id, r.seq.upper(), False, [], "fasta") for r in read_fasta(path)]
        if assume_circular is None:
            assume_circular = any(w in " ".join(r.description for r in read_fasta(path)).lower() for w in ("circular", "plasmid"))
    else:
        raise ValueError(f"{path.name}: not a reference format (expected FASTA, GenBank or SnapGene).")
    if assume_circular is not None:
        for r in refs:
            if r.source_format == "fasta":
                r.circular = assume_circular
    if not refs or not refs[0].seq:
        raise ValueError(f"{path.name}: no sequence found.")
    return refs
