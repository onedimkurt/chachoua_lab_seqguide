"""Minimal, dependency-free sequence file readers."""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def open_text(path: str | Path) -> io.TextIOBase:
    """Open plain or gzip-compressed text, detected by magic bytes (not extension)."""
    path = Path(path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def is_gzipped(path: str | Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


@dataclass
class Record:
    id: str
    seq: str
    qual: Optional[str] = None  # None for FASTA
    description: str = ""


def read_fasta(path: str | Path) -> Iterator[Record]:
    name, desc, chunks = None, "", []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield Record(name, "".join(chunks), None, desc)
                header = line[1:].strip()
                parts = header.split(None, 1)
                name = parts[0] if parts else ""
                desc = parts[1] if len(parts) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield Record(name, "".join(chunks), None, desc)


def read_fastq(path: str | Path) -> Iterator[Record]:
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            header = header.rstrip("\r\n")
            if not header:
                continue
            if not header.startswith("@"):
                raise ValueError(f"Not a valid FASTQ record header: {header[:50]!r}")
            seq = fh.readline().rstrip("\r\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\r\n")
            if not plus.startswith("+"):
                raise ValueError("Malformed FASTQ: expected '+' line")
            parts = header[1:].split(None, 1)
            yield Record(parts[0] if parts else "", seq, qual, parts[1] if len(parts) > 1 else "")


def sniff_format(path: str | Path) -> str:
    """Return 'fastq', 'fasta', 'genbank', 'snapgene', 'bam', 'sam', 'vcf' or 'unknown'."""
    path = Path(path)
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head[:2] == b"\x1f\x8b":
        # gzip; BAM is BGZF (also gzip) so peek inside
        with gzip.open(path, "rb") as g:
            inner = g.read(4)
        if inner == b"BAM\x01":
            return "bam"
        with open_text(path) as fh:
            first = _first_nonempty(fh)
        return _classify_text(first)
    if head and head[0] == 0x09:  # SnapGene .dna files start with byte 0x09
        if path.suffix.lower() == ".dna":
            return "snapgene"
    with open_text(path) as fh:
        first = _first_nonempty(fh)
    return _classify_text(first)


def _first_nonempty(fh) -> str:
    for line in fh:
        if line.strip():
            return line
    return ""


def _classify_text(first: str) -> str:
    if first.startswith("@HD") or first.startswith("@SQ"):
        return "sam"
    if first.startswith("@"):
        return "fastq"
    if first.startswith(">"):
        return "fasta"
    if first.startswith("LOCUS"):
        return "genbank"
    if first.startswith("##fileformat=VCF"):
        return "vcf"
    return "unknown"


def write_fasta(records: list[tuple[str, str]], path: str | Path, width: int = 70) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
