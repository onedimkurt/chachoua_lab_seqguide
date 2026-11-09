"""Step 1 of the wizard: "What data do I have?"

Everything here is pure Python + numpy so it runs on any laptop without
bioinformatics tools installed.
"""

from __future__ import annotations

import os
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from .seqio import is_gzipped, read_fasta, read_fastq, sniff_format

DEFAULT_MAX_READS = 500_000

# Common paired-end naming patterns, most specific first.
_PAIR_PATTERNS = [
    re.compile(r"^(?P<sample>.+?)(?P<sep>[._-])R(?P<mate>[12])(?P<rest>[._-]\d{3})?(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
    re.compile(r"^(?P<sample>.+?)(?P<sep>[._-])(?P<mate>[12])(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
]
_FASTQ_EXT = re.compile(r"\.f(ast)?q(\.gz)?$", re.I)


@dataclass
class FileReport:
    path: str
    name: str
    format: str
    compressed: bool
    size_mb: float
    n_reads: int = 0
    total_bases: int = 0
    min_len: int = 0
    max_len: int = 0
    mean_len: float = 0.0
    median_len: float = 0.0
    n50: int = 0
    gc_percent: float = 0.0
    mean_q: Optional[float] = None
    q30_percent: Optional[float] = None
    quality_encoding: Optional[str] = None
    platform: str = "unknown"          # "illumina", "nanopore", "pacbio", "long", "short", "unknown"
    platform_reason: str = ""
    sampled: bool = False              # True if only the first max_reads were read
    estimated_total_reads: Optional[int] = None
    length_hist: dict = field(default_factory=dict)  # {"edges": [...], "counts": [...]}
    plasmid_size_estimate: Optional[int] = None
    dimer_peak: bool = False
    first_ids: list = field(default_factory=list)
    problems: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class SampleGroup:
    sample: str
    files: list           # list[FileReport]
    layout: str           # "paired", "single", "long", "hybrid", "reference", "unknown"
    pair_check: Optional[dict] = None  # result of check_pair

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


@dataclass
class DatasetSummary:
    samples: list          # list[SampleGroup]
    references: list       # list[FileReport] for FASTA/GenBank/SnapGene references
    other: list            # unrecognised files

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------- per file

def n50(lengths: np.ndarray) -> int:
    if lengths.size == 0:
        return 0
    s = np.sort(lengths)[::-1]
    c = np.cumsum(s)
    return int(s[np.searchsorted(c, c[-1] / 2.0)])


def guess_platform(lengths: np.ndarray, first_header: str, mean_q: Optional[float]) -> tuple[str, str]:
    if lengths.size == 0:
        return "unknown", "No reads found."
    med = float(np.median(lengths))
    mx = int(lengths.max())
    h = first_header
    if "runid=" in h or re.search(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", h):
        return "nanopore", "Read names have the Oxford Nanopore format (run ID / UUID)."
    if re.match(r"^m\d+[eU]?_\d+_\d+/\d+", h):
        return "pacbio", "Read names have the PacBio format."
    if h.count(":") >= 4 and mx <= 600:
        return "illumina", f"Read names have the Illumina format and reads are short (up to {mx} bp)."
    if med >= 800 or mx > 2000:
        return "long", f"Reads are long (median {med:.0f} bp, longest {mx} bp), typical of Nanopore/PacBio."
    if mx <= 600:
        return "short", f"All reads are short (up to {mx} bp), typical of Illumina."
    return "unknown", "Read lengths do not clearly match a platform."


def estimate_circular_size(lengths: np.ndarray, bin_size: int = 100) -> tuple[Optional[int], bool]:
    """For whole-plasmid long-read data the read-length histogram peaks at the plasmid size.

    Returns (size_estimate, dimer_peak_present).
    """
    lengths = lengths[lengths >= 1000]
    if lengths.size < 20:
        return None, False
    edges = np.arange(0, int(lengths.max()) + bin_size * 2, bin_size)
    counts, edges = np.histogram(lengths, bins=edges)
    peak = int(np.argmax(counts))
    # A peak must be meaningful: at least 3x the median non-empty bin
    nonzero = counts[counts > 0]
    if counts[peak] < 3 * np.median(nonzero):
        return None, False
    lo, hi = edges[peak], edges[peak + 1]
    in_peak = lengths[(lengths >= lo - bin_size) & (lengths < hi + bin_size)]
    size = int(np.median(in_peak))
    # Dimer: bin around 2x size holds >= 5% of the main peak
    d_lo, d_hi = 1.9 * size, 2.1 * size
    dimer = int(((lengths >= d_lo) & (lengths <= d_hi)).sum())
    return size, dimer >= max(5, 0.05 * in_peak.size)


def _length_hist(lengths: np.ndarray, nbins: int = 40) -> dict:
    if lengths.size == 0:
        return {"edges": [], "counts": []}
    lo, hi = int(lengths.min()), int(lengths.max())
    if lo == hi:
        return {"edges": [lo, hi + 1], "counts": [int(lengths.size)]}
    counts, edges = np.histogram(lengths, bins=min(nbins, hi - lo + 1))
    return {"edges": [float(e) for e in edges], "counts": [int(c) for c in counts]}


def inspect_file(path: str | Path, max_reads: int = DEFAULT_MAX_READS) -> FileReport:
    path = Path(path)
    fmt = sniff_format(path)
    rep = FileReport(
        path=str(path),
        name=path.name,
        format=fmt,
        compressed=is_gzipped(path),
        size_mb=round(os.path.getsize(path) / 1e6, 2),
    )
    if fmt not in ("fastq", "fasta"):
        if fmt in ("genbank", "snapgene"):
            rep.platform_reason = "Annotated sequence file (usually a plasmid map or reference)."
        elif fmt in ("bam", "sam"):
            rep.problems.append("This is an alignment file (BAM/SAM). SeqGuide currently starts from raw reads (FASTQ).")
        elif fmt == "vcf":
            rep.problems.append("This is a variant file (VCF), i.e. already-analysed results.")
        else:
            rep.problems.append("File type not recognised. Expected FASTQ, FASTA, GenBank (.gb) or SnapGene (.dna).")
        return rep

    reader = read_fastq if fmt == "fastq" else read_fasta
    lengths: list[int] = []
    gc = 0
    qual_sum = 0
    q30 = 0
    qual_min, qual_max = 255, 0
    qbuf: list[bytes] = []
    qbuf_len = 0
    first_header = ""
    bytes_read_estimate = 0

    def flush():
        nonlocal qual_sum, q30, qual_min, qual_max, qbuf, qbuf_len
        if not qbuf:
            return
        arr = np.frombuffer(b"".join(qbuf), dtype=np.uint8)
        qual_sum += int(arr.sum(dtype=np.int64))
        q30 += int((arr >= 63).sum())
        qual_min = min(qual_min, int(arr.min()))
        qual_max = max(qual_max, int(arr.max()))
        qbuf, qbuf_len = [], 0

    try:
        for i, rec in enumerate(reader(path)):
            if i >= max_reads:
                rep.sampled = True
                break
            if i == 0:
                first_header = rec.id + " " + rec.description
            if i < 1000:
                rep.first_ids.append(rec.id)
            L = len(rec.seq)
            lengths.append(L)
            s = rec.seq.upper()
            gc += s.count("G") + s.count("C")
            bytes_read_estimate += 2 * L + len(rec.id) + 6
            if rec.qual is not None:
                if len(rec.qual) != L:
                    rep.problems.append(f"Read {rec.id}: sequence and quality lengths differ (file may be truncated).")
                    break
                qbuf.append(rec.qual.encode("ascii", "replace"))
                qbuf_len += L
                if qbuf_len > 5_000_000:
                    flush()
        flush()
    except (ValueError, EOFError, OSError) as exc:
        rep.problems.append(f"File could not be read completely: {exc}. It may be truncated or corrupted.")

    arr = np.asarray(lengths, dtype=np.int64)
    rep.n_reads = int(arr.size)
    rep.total_bases = int(arr.sum())
    if arr.size:
        rep.min_len, rep.max_len = int(arr.min()), int(arr.max())
        rep.mean_len = round(float(arr.mean()), 1)
        rep.median_len = float(np.median(arr))
        rep.n50 = n50(arr)
        rep.gc_percent = round(100.0 * gc / max(1, rep.total_bases), 1)
        rep.length_hist = _length_hist(arr)
    if fmt == "fastq" and rep.total_bases:
        if qual_min < 59:
            rep.quality_encoding = "phred33"
            offset = 33
        elif qual_min >= 64:
            rep.quality_encoding = "phred64 (old Illumina)"
            offset = 64
            rep.problems.append("Quality scores use the old Phred+64 encoding; they will be converted.")
        else:
            rep.quality_encoding = "phred33"
            offset = 33
        rep.mean_q = round(qual_sum / rep.total_bases - offset, 1)
        rep.q30_percent = round(100.0 * q30 / rep.total_bases, 1) if offset == 33 else None

    if rep.sampled:
        # Extrapolate total read count from file size
        raw_mb = rep.size_mb * (3.5 if rep.compressed else 1.0)
        per_read = bytes_read_estimate / max(1, rep.n_reads)
        rep.estimated_total_reads = int(raw_mb * 1e6 / per_read)

    if fmt == "fastq":
        rep.platform, rep.platform_reason = guess_platform(arr, first_header, rep.mean_q)
    else:
        rep.platform, rep.platform_reason = "unknown", "FASTA files have no quality scores; this looks like assembled sequence or a reference."
        if rep.n_reads <= 50 and rep.max_len >= 1000:
            rep.platform_reason = "A few long sequences without quality scores: most likely a reference or an assembly."

    if rep.platform in ("nanopore", "pacbio", "long"):
        rep.plasmid_size_estimate, rep.dimer_peak = estimate_circular_size(arr)

    if rep.n_reads == 0:
        rep.problems.append("The file contains no reads.")
    return rep


# ---------------------------------------------------------------- pairing

def pair_key(filename: str) -> Optional[tuple[str, str]]:
    """Return (sample, mate) if the filename looks like one half of a pair."""
    for pat in _PAIR_PATTERNS:
        m = pat.match(filename)
        if m:
            return m.group("sample"), m.group("mate")
    return None


def _norm_id(read_id: str) -> str:
    return re.sub(r"/[12]$", "", read_id)


def check_pair(r1: FileReport, r2: FileReport) -> dict:
    """Verify two files are mates: same read names, in the same order, same count."""
    ids1 = [_norm_id(x) for x in r1.first_ids]
    ids2 = [_norm_id(x) for x in r2.first_ids]
    n = min(len(ids1), len(ids2))
    matching = sum(1 for a, b in zip(ids1[:n], ids2[:n]) if a == b)
    same_order = n > 0 and matching == n
    counts_equal = r1.n_reads == r2.n_reads
    ok = same_order and (counts_equal or (r1.sampled and r2.sampled))
    if ok:
        msg = "The two files are proper mates: read names match one-to-one."
    elif not same_order:
        msg = (f"Only {matching} of the first {n} read names match between the two files. "
               "They may not belong together, or one file was filtered separately.")
    else:
        msg = (f"Read names match but the files have different read counts ({r1.n_reads} vs {r2.n_reads}). "
               "One file may be truncated.")
    return {"ok": ok, "matching_ids": matching, "checked": n, "counts_equal": counts_equal, "message": msg}


def sample_name(filename: str) -> str:
    k = pair_key(filename)
    if k:
        name = k[0]
    else:
        name = _FASTQ_EXT.sub("", filename)
        name = re.sub(r"\.(fa|fasta|fna|gb|gbk|genbank|dna)(\.gz)?$", "", name, flags=re.I)
    # Illumina adds sample-sheet number and lane: BA137_S1_L001 -> BA137
    return re.sub(r"_S\d+(_L\d{3})?$", "", name)


def inspect_dataset(paths: list[str | Path], max_reads: int = DEFAULT_MAX_READS) -> DatasetSummary:
    reports = [inspect_file(p, max_reads=max_reads) for p in paths]
    reads = [r for r in reports if r.format == "fastq"]
    refs = [r for r in reports if r.format in ("fasta", "genbank", "snapgene")]
    other = [r for r in reports if r not in reads and r not in refs]

    groups: dict[str, list[FileReport]] = {}
    for r in reads:
        groups.setdefault(sample_name(r.name), []).append(r)

    samples = []
    for sname, files in sorted(groups.items()):
        files.sort(key=lambda f: f.name)
        shorts = [f for f in files if f.platform in ("illumina", "short")]
        longs = [f for f in files if f.platform in ("nanopore", "pacbio", "long")]
        pair = None
        if len(shorts) == 2 and all(pair_key(f.name) for f in shorts):
            pair = check_pair(shorts[0], shorts[1])
            layout = "hybrid" if longs else "paired"
        elif len(shorts) == 1 and not longs:
            layout = "single"
        elif longs and not shorts:
            layout = "long"
        elif longs and len(shorts) == 1:
            layout = "hybrid"
        else:
            layout = "unknown"
        samples.append(SampleGroup(sample=sname, files=files, layout=layout, pair_check=pair))
    return DatasetSummary(samples=samples, references=refs, other=other)


def depth(total_bases: int, target_size: Optional[int]) -> Optional[float]:
    if not target_size:
        return None
    return round(total_bases / target_size, 1)
