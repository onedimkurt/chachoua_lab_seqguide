"""Compare an assembled plasmid with its reference map.

Pure-Python, plasmid-scale (up to ~100 kb) aligner:
  1. orient (forward / reverse complement) and rotate the circular assembly so
     it starts where the reference starts;
  2. anchor the two sequences with unique shared k-mers;
  3. fill the small gaps between anchors with an affine-gap global alignment;
  4. turn the result into a list of differences (variants).

Large gaps are reported as structural differences (insertion, deletion,
inversion, replaced region) instead of being aligned base by base.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Optional

from .seqio import revcomp

K = 21
MAX_DP_CELLS = 400_000   # gaps bigger than this are reported as structural events


@dataclass
class Variant:
    ref_pos: int            # 0-based position in the reference (insertions: inserted *before* this position)
    ref: str                # reference allele ("" for insertion)
    alt: str                # sample allele ("" for deletion)
    qry_pos: int            # 0-based position in the (rotated) assembly
    kind: str = ""          # SNP, MNP, INS, DEL, COMPLEX, INVERSION
    structural: bool = False
    # filled in later by effects.annotate / read support
    features: list = field(default_factory=list)
    effect: str = ""
    severity: str = ""
    homopolymer: int = 0
    depth: Optional[int] = None
    flags: list = field(default_factory=list)

    def __post_init__(self):
        if not self.kind:
            self.kind = classify(self.ref, self.alt)

    @property
    def size(self) -> int:
        return max(len(self.ref), len(self.alt))

    def describe(self) -> str:
        p = self.ref_pos + 1
        if self.kind == "SNP":
            return f"{self.ref}{p}{self.alt}"
        if self.kind == "INS":
            ins = self.alt if len(self.alt) <= 12 else f"{len(self.alt)} bp"
            return f"insertion of {ins} after position {self.ref_pos}"
        if self.kind == "DEL":
            if len(self.ref) == 1:
                return f"deletion of {self.ref} at position {p}"
            return f"deletion of {len(self.ref)} bp at {p}–{p + len(self.ref) - 1}"
        if self.kind == "INVERSION":
            return f"inverted region {p}–{p + len(self.ref) - 1} ({len(self.ref)} bp)"
        if self.kind == "MNP":
            return f"{self.ref}>{self.alt} at {p}–{p + len(self.ref) - 1}"
        return f"region {p}–{p + max(len(self.ref), 1) - 1} ({len(self.ref)} bp) replaced by {len(self.alt)} bp"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["description"] = self.describe()
        d["position"] = self.ref_pos + 1
        return d


def classify(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNP"
    if not ref:
        return "INS"
    if not alt:
        return "DEL"
    if len(ref) == len(alt):
        return "MNP"
    return "COMPLEX"


@dataclass
class Comparison:
    ref_name: str
    ref_len: int
    qry_name: str
    qry_len: int
    strand: str              # "+" or "-" (assembly was reverse-complemented)
    rotation: int            # assembly position that became position 1
    identity: float          # % identical over aligned positions
    ref_covered: float       # % of reference matched by the assembly
    variants: list           # list[Variant]
    rotated_seq: str = ""
    notes: list = field(default_factory=list)
    ref_to_qry: list = field(default_factory=list)  # sparse anchor map [(ref_pos, qry_pos)] for coordinate conversion
    trimmed_overlap: int = 0

    def to_dict(self) -> dict:
        d = asdict(self)
        d["variants"] = [v.to_dict() for v in self.variants]
        d.pop("rotated_seq")
        d.pop("ref_to_qry")
        return d


# ------------------------------------------------------------------ k-mers

def _kmer_positions(seq: str, k: int, circular: bool) -> dict:
    s = seq + (seq[: k - 1] if circular else "")
    pos: dict = {}
    for i in range(len(s) - k + 1):
        km = s[i : i + k]
        if "N" in km:
            continue
        pos.setdefault(km, []).append(i)
    return pos


def _unique(pos: dict) -> dict:
    return {k: v[0] for k, v in pos.items() if len(v) == 1}


# ------------------------------------------------------------------ orientation

def orient_and_rotate(qry: str, ref: str, circular: bool = True, k: int = K) -> tuple:
    """Return list of candidate (strand, rotation, rotated_seq), best guess first."""
    ref_u = _unique(_kmer_positions(ref, k, False))
    cands = []
    for strand, s in (("+", qry), ("-", revcomp(qry))):
        q_u = _unique(_kmer_positions(s, k, circular))
        hits = sorted((rp, q_u[km]) for km, rp in ref_u.items() if km in q_u)
        cands.append((len(hits), strand, s, hits))
    cands.sort(key=lambda c: -c[0])
    n_hits, strand, s, hits = cands[0]
    if n_hits == 0:
        return [(strand, 0, s)], 0
    out = []
    if circular:
        L = len(s)
        offsets = Counter((q - r) % L for r, q in hits)
        first_r, first_q = hits[0]
        options = [(first_q - first_r) % L, first_q % L, offsets.most_common(1)[0][0]]
        seen = set()
        for o in options:
            if o in seen:
                continue
            seen.add(o)
            out.append((strand, o, s[o:] + s[:o]))
    else:
        out.append((strand, 0, s))
    return out, n_hits


# ------------------------------------------------------------------ gap alignment

def _gotoh(a: str, b: str, match=2, mismatch=-3, gap_open=-6, gap_ext=-1) -> list:
    """Global affine-gap alignment. Returns list of (a_char or '', b_char or '') columns."""
    n, m = len(a), len(b)
    if n == 0:
        return [("", c) for c in b]
    if m == 0:
        return [(c, "") for c in a]
    NEG = float("-inf")
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]   # consumes a (deletion from reference)
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]   # consumes b (insertion)
    PM = [[0] * (m + 1) for _ in range(n + 1)]
    PX = [[0] * (m + 1) for _ in range(n + 1)]
    PY = [[0] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = gap_open + (i - 1) * gap_ext
        PX[i][0] = 1 if i > 1 else 0
    for j in range(1, m + 1):
        Y[0][j] = gap_open + (j - 1) * gap_ext
        PY[0][j] = 2 if j > 1 else 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        Mi, Mp, Xi, Xp, Yi, Yp = M[i], M[i - 1], X[i], X[i - 1], Y[i], Y[i - 1]
        PMi, PXi, PYi = PM[i], PX[i], PY[i]
        for j in range(1, m + 1):
            s = match if ai == b[j - 1] else mismatch
            m0, m1, m2 = Mp[j - 1], Xp[j - 1], Yp[j - 1]
            if m0 >= m1 and m0 >= m2:
                Mi[j], PMi[j] = m0 + s, 0
            elif m1 >= m2:
                Mi[j], PMi[j] = m1 + s, 1
            else:
                Mi[j], PMi[j] = m2 + s, 2
            x0, x1, x2 = Mp[j] + gap_open, Xp[j] + gap_ext, Yp[j] + gap_open
            if x1 >= x0 and x1 >= x2:
                Xi[j], PXi[j] = x1, 1
            elif x0 >= x2:
                Xi[j], PXi[j] = x0, 0
            else:
                Xi[j], PXi[j] = x2, 2
            y0, y1, y2 = Mi[j - 1] + gap_open, Xi[j - 1] + gap_open, Yi[j - 1] + gap_ext
            if y2 >= y0 and y2 >= y1:
                Yi[j], PYi[j] = y2, 2
            elif y0 >= y1:
                Yi[j], PYi[j] = y0, 0
            else:
                Yi[j], PYi[j] = y1, 1
    scores = (M[n][m], X[n][m], Y[n][m])
    state = max(range(3), key=lambda t: scores[t])
    cols = []
    i, j = n, m
    while i > 0 or j > 0:
        if state == 0:
            p = PM[i][j]
            cols.append((a[i - 1], b[j - 1]))
            i, j = i - 1, j - 1
        elif state == 1:
            p = PX[i][j]
            cols.append((a[i - 1], ""))
            i -= 1
        else:
            p = PY[i][j]
            cols.append(("", b[j - 1]))
            j -= 1
        state = p
    cols.reverse()
    return cols


def _cols_to_variants(cols: list, r0: int, q0: int) -> list:
    """Turn alignment columns into merged variants. r0/q0 are the global start coordinates."""
    out = []
    r, q = r0, q0
    cur = None  # [ref_pos, ref, alt, qry_pos]
    for a, b in cols:
        if a and b and a == b:
            if cur:
                out.append(Variant(cur[0], cur[1], cur[2], cur[3]))
                cur = None
            r += 1
            q += 1
            continue
        if cur is None:
            cur = [r, "", "", q]
        cur[1] += a
        cur[2] += b
        if a:
            r += 1
        if b:
            q += 1
    if cur:
        out.append(Variant(cur[0], cur[1], cur[2], cur[3]))
    # split "mismatch then gap" blocks only if clearly separable; keep COMPLEX otherwise
    return out


def _looks_inverted(a: str, b: str, k: int = 15) -> bool:
    if min(len(a), len(b)) < 50 or not 0.8 < len(a) / len(b) < 1.25:
        return False
    rb = revcomp(b)
    ka = {a[i:i + k] for i in range(len(a) - k + 1)}
    kb = [rb[i:i + k] for i in range(len(rb) - k + 1)]
    shared = sum(1 for x in kb if x in ka)
    return shared / max(1, len(kb)) > 0.5


def _align_gap(ref: str, qry: str, r0: int, r1: int, q0: int, q1: int) -> list:
    a, b = ref[r0:r1], qry[q0:q1]
    if not a and not b:
        return []
    if not a or not b or len(a) * len(b) <= MAX_DP_CELLS:
        return _cols_to_variants(_gotoh(a, b), r0, q0)
    v = Variant(r0, a, b, q0, kind="INVERSION" if _looks_inverted(a, b) else "COMPLEX", structural=True)
    return [v]


# ------------------------------------------------------------------ chaining

def _chain(hits: list) -> list:
    """Longest chain of anchors increasing in both ref and qry (O(n log n) LIS)."""
    import bisect
    hits = sorted(hits)
    tails, tails_idx, prev = [], [], [-1] * len(hits)
    for idx, (_, q) in enumerate(hits):
        p = bisect.bisect_left(tails, q)
        if p == len(tails):
            tails.append(q)
            tails_idx.append(idx)
        else:
            tails[p] = q
            tails_idx[p] = idx
        prev[idx] = tails_idx[p - 1] if p > 0 else -1
    chain = []
    i = tails_idx[-1] if tails_idx else -1
    while i >= 0:
        chain.append(hits[i])
        i = prev[i]
    return chain[::-1]


def align_linear(ref: str, qry: str, k: int = K) -> tuple:
    """Align two (already oriented & rotated) sequences end to end.

    Returns (variants, matched_bases, anchor_map).
    """
    ref_u = _unique(_kmer_positions(ref, k, False))
    q_u = _unique(_kmer_positions(qry, k, False))
    hits = [(rp, q_u[km]) for km, rp in ref_u.items() if km in q_u]
    chain = _chain(hits)
    variants, anchors = [], []
    r = q = 0
    matched = 0
    for ra, qa in chain:
        if qa - ra == q - r and ra <= r:          # continues current diagonal
            ext = ra + k - r
            if ext > 0:
                matched += ext
                r, q = r + ext, q + ext
            continue
        if ra >= r and qa >= q:
            variants.extend(_align_gap(ref, qry, r, ra, q, qa))
            anchors.append((ra, qa))
            matched += k
            r, q = ra + k, qa + k
        # else: anchor overlaps an earlier one on a different diagonal; skip it
    variants.extend(_align_gap(ref, qry, r, len(ref), q, len(qry)))
    # count matches inside DP-aligned gaps: total ref length minus bases involved in variants
    return variants, matched, anchors


def left_normalize(v: Variant, ref: str) -> Variant:
    """Shift simple indels to their left-most equivalent position (standard representation)."""
    if v.kind not in ("INS", "DEL") or v.structural:
        return v
    seq = v.ref if v.kind == "DEL" else v.alt
    pos, qpos = v.ref_pos, v.qry_pos
    while pos > 0 and ref[pos - 1] == seq[-1]:
        seq = seq[-1] + seq[:-1]
        pos -= 1
        qpos -= 1
    if v.kind == "DEL":
        v.ref = seq
    else:
        v.alt = seq
    v.ref_pos, v.qry_pos = pos, qpos
    return v


def homopolymer_len(ref: str, v: Variant) -> int:
    """Length of the longest single-letter run touching the variant in the reference."""
    if v.kind not in ("INS", "DEL"):
        return 0
    seq = v.ref if v.kind == "DEL" else v.alt
    if len(set(seq)) != 1:
        return 0
    base = seq[0]
    i = v.ref_pos
    j = v.ref_pos + len(v.ref)
    while i > 0 and ref[i - 1] == base:
        i -= 1
    while j < len(ref) and ref[j] == base:
        j += 1
    return j - i


def trim_circular_overlap(seq: str, min_overlap: int = 20, max_overlap: int = 500) -> tuple:
    """Remove the duplicated sequence that de Bruijn assemblers (SPAdes, Velvet...) leave at both ends of a
    circular contig (usually k-1 or k bp). Returns (trimmed_seq, overlap_length)."""
    best = 0
    for L in range(min(max_overlap, len(seq) // 2), min_overlap - 1, -1):
        if seq[:L] == seq[-L:]:
            best = L
            break
    return (seq[:-best], best) if best else (seq, 0)


def compare(qry: str, ref: str, circular: bool = True, qry_name: str = "assembly", ref_name: str = "reference") -> Comparison:
    qry, ref = qry.upper(), ref.upper()
    overlap = 0
    if circular:
        qry, overlap = trim_circular_overlap(qry)
    cands, n_hits = orient_and_rotate(qry, ref, circular)
    best = None
    for strand, rot, seq in cands:
        variants, matched, anchors = align_linear(ref, seq)
        edit = sum(max(len(v.ref), len(v.alt)) for v in variants)
        if best is None or edit < best[0]:
            best = (edit, strand, rot, seq, variants, anchors)
    edit, strand, rot, seq, variants, anchors = best
    variants = [left_normalize(v, ref) for v in variants]
    for v in variants:
        v.homopolymer = homopolymer_len(ref, v)

    unmatched_ref = sum(len(v.ref) for v in variants)
    mism = sum(len(v.ref) for v in variants if v.kind in ("SNP", "MNP"))
    aligned = max(1, len(ref) - unmatched_ref + mism)
    identity = 100.0 * (aligned - mism) / aligned
    covered = 100.0 * (len(ref) - unmatched_ref + mism) / max(1, len(ref))
    notes = []
    if overlap:
        notes.append(f"The assembly's two ends overlapped by {overlap} bp, which is how assemblers like SPAdes mark a circle. "
                     "The duplicate copy was removed before comparison (this is not a real insertion).")
    if n_hits == 0:
        notes.append("The assembly shares no sequence with the reference — this is probably a different plasmid.")
    ratio = len(qry) / max(1, len(ref))
    if 1.9 < ratio < 2.1:
        notes.append("The assembly is about twice the reference length; it may be a dimer (two copies joined).")
    return Comparison(ref_name, len(ref), qry_name, len(qry), strand, rot, round(identity, 3), round(covered, 2),
                      variants, seq, notes, anchors, overlap)
