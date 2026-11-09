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
                   max_reads: int | None = None, circular: bool = True,
                   alignments: list | None = None) -> ReadCheck:
    """records: iterable of (sequence, quality-or-None[, name]). Quality None = every base trusted.
    circular=False for linear targets (amplicons): reads must fit inside the target.
    alignments: if a list is given, one tuple per placed read is appended:
        (pos0, reverse, read_seq_as_aligned, qual_or_None, left_clip, aligned_len, mapq, name)"""
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
    for item in records:
        seq, qual = item[0], item[1]
        if max_reads and n >= max_reads:
            break
        n += 1
        name = item[2] if len(item) > 2 else f"read{n}"
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
                    lead = 0
                    # 5' end that does not match (read starts inside an insertion / junction): clip it
                    for j in range(len(mm) - 1, 0, -1):
                        if (j + 1) / (mm[j] + 1) > 0.2:
                            lead = mm[j] + 1
                            break
                    if lead:
                        mm = [m for m in mm if m >= lead]
                    # 3' end that does not match (adapter, junction, insertion): clip it
                    clip_pos = None
                    j = next((j for j in range(len(mm) - 1) if (len(mm) - j) / (alen - mm[j]) > 0.2), None)
                    if j is not None:
                        clip_pos = (start + mm[j]) % L
                        mm, alen = mm[:j], mm[j]
                    if alen - lead < 40 or len(mm) > allowed:
                        continue
                    if lead:
                        mm = [m - lead for m in mm]
                        s2, q2 = s[lead:alen], (q[lead:alen] if q else None)
                        cands.append((alen - lead, -len(mm), rng.random(), (start + lead) % L, s2, q2, mm, clip_pos,
                                      strand, s, q, lead))
                        continue
                    cands.append((alen, -len(mm), rng.random(), start, s, q, mm, clip_pos, strand, s, q, 0))
        if not cands:
            continue
        cands.sort(reverse=True)
        best = cands[0]
        alen, start, s, q, mm = best[0], best[3], best[4], best[5], best[6]
        if alignments is not None:
            tie = len(cands) > 1 and cands[1][0] == best[0] and cands[1][1] == best[1] and cands[1][3] != best[3]
            alignments.append((start, best[8] == -1, best[9], best[10], best[11], alen, 0 if tie else 60, name))
        if best[7] is not None:
            clips[best[7]] = clips.get(best[7], 0) + 1
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
