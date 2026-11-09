"""Pure-Python read support check for small circular references (plasmids).

Maps short reads to a plasmid with exact k-mer seeds + ungapped extension and builds a pileup:
per-position depth and alternative-base counts. Used when minimap2/samtools are not available,
and as an independent check of an assembly or a reference.

Limitations (by design, to stay simple): ungapped only, so reads spanning an indel are not
counted. Reads from repeats (e.g. the two identical LTRs of lentiviral plasmids) are placed on one
copy chosen at random, as bwa/minimap2 do. Adapter or low-quality 3' tails are clipped.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .seqio import read_fastq, revcomp


@dataclass
class ReadCheck:
    n_reads: int
    n_mapped: int
    depth: np.ndarray
    alt_counts: dict = field(default_factory=dict)   # pos0 -> {base: count}
    lowq: np.ndarray = None                           # bases below the quality cut-off, per position

    @property
    def pct_mapped(self) -> float:
        return round(100.0 * self.n_mapped / max(1, self.n_reads), 1)

    def variants(self, min_af: float = 0.05, min_alt: int = 3) -> list:
        out = []
        for pos, counts in sorted(self.alt_counts.items()):
            d = int(self.depth[pos] - (self.lowq[pos] if self.lowq is not None else 0))
            for base, n in counts.items():
                if n >= min_alt and d and n / d >= min_af:
                    out.append({"pos": pos + 1, "alt": base, "alt_reads": n, "depth": d, "af": round(n / d, 3)})
        return out


def _iter_fastq(path):
    for rec in read_fastq(path):
        yield rec.seq, rec.qual


def map_and_pileup(records, ref: str, k: int = 19, max_mm_frac: float = 0.05, min_bq: int = 20,
                   max_reads: int | None = None, circular: bool = True) -> ReadCheck:
    """records: iterable of (sequence, quality-or-None). Quality None = every base trusted.
    circular=False for linear targets (amplicons): reads must fit inside the target."""
    ref = ref.upper()
    L = len(ref)
    ext = ref + ref[:400] if circular else ref  # circular: reads may run across the origin
    index: dict = {}
    for i in range(L):
        index.setdefault(ext[i:i + k], []).append(i)
    rng = np.random.default_rng(1)
    depth_diff = np.zeros(L + 1, dtype=np.int64)
    wrap = np.zeros(L + 1, dtype=np.int64)
    alt: dict = {}
    clips: dict = {}
    lowq = np.zeros(L, dtype=np.int64)
    n = mapped = 0
    for seq, qual in records:
        if max_reads and n >= max_reads:
            break
        n += 1
        seq = seq.upper()
        rl = len(seq)
        if rl < k + 5:
            continue
        rq = qual if qual else None
        allowed = max(2, int(max_mm_frac * rl))
        cands = []
        seen = set()
        for strand, s, q in ((1, seq, rq), (-1, revcomp(seq), rq[::-1] if rq else None)):
            for off in (0, rl // 2 - k // 2, rl - k):
                for p in index.get(s[off:off + k], ())[:8]:
                    if not circular and p - off < 0:
                        continue
                    start = (p - off) % L
                    if (strand, start) in seen:
                        continue
                    seen.add((strand, start))
                    target = ext[start:start + rl]
                    alen = rl
                    if len(target) < rl:
                        if circular or len(target) < 40:
                            continue
                        alen = len(target)          # read runs past the end of a linear target
                    mm = [i for i in range(alen) if s[i] != target[i]]
                    if len(mm) > allowed:
                        j = next((j for j in range(len(mm)) if (len(mm) - j) / (rl - mm[j]) > 0.2), None)
                        if j is None or j > allowed or mm[j] < 40:
                            continue
                        mm, alen = mm[:j], mm[j]
                        clipped_at = (start + alen) % L
                        clips[clipped_at] = clips.get(clipped_at, 0) + 1
                    cands.append((alen, -len(mm), rng.random(), start, s, q, mm))
        if not cands:
            continue
        best = max(cands)
        alen, start, s, q, mm = best[0], best[3], best[4], best[5], best[6]
        mapped += 1
        end = start + alen
        if end <= L:
            depth_diff[start] += 1
            depth_diff[end] -= 1
        else:
            depth_diff[start] += 1
            depth_diff[L] -= 1
            wrap[0] += 1
            wrap[end - L] -= 1
        if q:
            qa = np.frombuffer(q[:alen].encode("ascii", "replace"), dtype=np.uint8)
            bad = np.nonzero(qa < 33 + min_bq)[0]
            if bad.size:
                np.add.at(lowq, (start + bad) % L, 1)
        for i in mm:
            if (q and ord(q[i]) - 33 < min_bq) or s[i] == "N":
                continue
            pos = (start + i) % L
            d = alt.setdefault(pos, {})
            d[s[i]] = d.get(s[i], 0) + 1
    depth = np.cumsum(depth_diff)[:L] + np.cumsum(wrap)[:L]
    rc = ReadCheck(n, mapped, depth, alt, lowq)
    rc.clips = clips
    return rc


def check_reads(fastq: str | Path, ref: str, k: int = 19, max_mm_frac: float = 0.05, min_bq: int = 20,
                max_reads: int | None = None) -> ReadCheck:
    return map_and_pileup(_iter_fastq(fastq), ref, k, max_mm_frac, min_bq, max_reads)


def pileup_reads(seqs: list, ref: str):
    """Pile up plain sequences (already cleaned) on a circular reference -> (depth array, alt counts)."""
    rc = map_and_pileup(((s, None) for s in seqs), ref)
    return rc.depth, rc.alt_counts
