"""A small, dependency-free de novo assembler for plasmids from short reads (Illumina).

Written in pure Python (+ the standard library) so it also runs inside a web browser
(Pyodide / stlite), where Unicycler and SPAdes cannot run.

Approach (the same ideas Unicycler/SPAdes use, scaled down to plasmids):
  1. clean reads: cut adapters, poly-G tails and low-quality 3' ends;
  2. count k-mers on both strands and keep "solid" k-mers (seen often enough to be real);
  3. build a de Bruijn graph and compact it into unitigs;
  4. clean the graph: remove dead-end tips, pop error bubbles (recording minor variants), drop
     low-coverage pieces (host DNA, errors);
  5. estimate how many copies of each unitig the plasmid contains from its coverage and find the
     circular walk that uses every unitig that many times (this resolves repeats such as the two
     LTR copies of a lentiviral plasmid);
  6. try several k values and keep the first that gives a clean circle;
  7. polish: place all reads on the result and fix any position where most reads disagree.

Results are always checked against the reads afterwards (depth along the plasmid).
"""

from __future__ import annotations

import statistics
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Optional

from .seqio import revcomp

ADAPTERS = ("AGATCGGAAGAGC", "CTGTCTCTTATACACATCT")  # Illumina TruSeq / Nextera


# ------------------------------------------------------------------ read cleaning

def clean_read(seq: str, qual: Optional[str], min_q: int = 20, min_len: int = 40) -> Optional[str]:
    s = seq.upper()
    cut = len(s)
    for a in ADAPTERS:
        i = s.find(a)
        if i != -1:
            cut = min(cut, i)
    # poly-G tails (2-colour chemistry: NextSeq / NovaSeq)
    j = cut
    while j > 0 and s[j - 1] == "G":
        j -= 1
    if cut - j >= 10:
        cut = j
    if qual:
        # trim 3' end while a 4-base window averages below min_q
        q = [ord(c) - 33 for c in qual[:cut]]
        while cut >= 4 and sum(q[cut - 4:cut]) / 4 < min_q:
            cut -= 1
    s = s[:cut]
    if len(s) < min_len or s.count("N") > 2:
        return None
    return s


# ------------------------------------------------------------------ graph

@dataclass
class Unitig:
    uid: int
    seq: str
    cov: float          # mean k-mer count
    twin: int = -1      # uid of the reverse-complement unitig
    succ: list = field(default_factory=list)
    pred: list = field(default_factory=list)

    @property
    def nkmers(self) -> int:
        return len(self.seq)  # placeholder, set properly in build (len - k + 1)


@dataclass
class AssemblyResult:
    sequence: Optional[str]             # main circular (or longest) contig
    circular: bool
    k: Optional[int]
    contigs: list                       # [{"name", "seq", "circular", "coverage", "copies"}]
    notes: list
    minor_variants: list                # bubbles popped: possible mixed population
    stats: dict


class DeBruijn:
    def __init__(self, counts: dict, k: int):
        self.k = k
        self.counts = counts
        self.unitigs: list = []
        self._build()

    # -- helpers
    def _succ(self, km: str) -> list:
        base = km[1:]
        return [base + b for b in "ACGT" if base + b in self.counts]

    def _pred(self, km: str) -> list:
        base = km[:-1]
        return [b + base for b in "ACGT" if b + base in self.counts]

    def _build(self) -> None:
        k = self.k
        counts = self.counts
        visited = set()
        seqs = []

        def extend(start: str) -> list:
            path = [start]
            visited.add(start)
            cur = start
            while True:
                s = self._succ(cur)
                if len(s) != 1:
                    break
                nxt = s[0]
                if len(self._pred(nxt)) != 1 or nxt in visited:
                    break
                path.append(nxt)
                visited.add(nxt)
                cur = nxt
            return path

        for km in counts:
            if km in visited:
                continue
            p = self._pred(km)
            if len(p) == 1 and len(self._succ(p[0])) == 1:
                continue  # interior k-mer; reached from its unitig start
            seqs.append(extend(km))
        # what is left are perfect cycles (e.g. a plasmid with no branching at all)
        for km in counts:
            if km not in visited:
                seqs.append(extend(km))

        self.unitigs = []
        end_index, start_index = {}, {}
        for path in seqs:
            seq = path[0] + "".join(x[-1] for x in path[1:])
            cov = sum(counts[x] for x in path) / len(path)
            u = Unitig(len(self.unitigs), seq, cov)
            self.unitigs.append(u)
            start_index[path[0]] = u.uid
            end_index[path[-1]] = u.uid
        for u in self.unitigs:
            last = u.seq[-k:]
            for nxt in self._succ(last):
                if nxt in start_index:
                    u.succ.append(start_index[nxt])
            first = u.seq[:k]
            for prv in self._pred(first):
                if prv in end_index:
                    u.pred.append(end_index[prv])
        owner = {}
        for u in self.unitigs:
            for i in range(len(u.seq) - k + 1):
                owner[u.seq[i:i + k]] = u.uid
        for u in self.unitigs:  # works for linear unitigs and for perfect cycles (any rotation)
            u.twin = owner.get(revcomp(u.seq[:k]), -1)

    def nk(self, u: Unitig) -> int:
        return len(u.seq) - self.k + 1


# ------------------------------------------------------------------ cleaning

def _remove(counts: dict, seq: str, k: int) -> None:
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        counts.pop(km, None)
        counts.pop(revcomp(km), None)


def _clean(counts: dict, k: int, minor: list, rounds: int = 6) -> DeBruijn:
    g = DeBruijn(counts, k)
    for _ in range(rounds):
        changed = False
        covs = [u.cov for u in g.unitigs for _ in range(max(1, g.nk(u) // 50))]
        base = statistics.median(covs) if covs else 0
        for u in g.unitigs:
            n = g.nk(u)
            dead_end = (not u.succ) or (not u.pred)
            # tips: short dead ends
            if dead_end and n < 2 * k and (u.succ or u.pred):
                _remove(counts, u.seq, k)
                changed = True
            # low coverage pieces (errors, host DNA, chimeric artefacts)
            elif u.cov < 0.15 * base or (n < 3 * k and u.cov < 0.35 * base):
                _remove(counts, u.seq, k)
                changed = True
        # bubbles: two unitigs with the same single predecessor and successor
        seen = set()
        for u in g.unitigs:
            if len(u.succ) != 2:
                continue
            a, b = (g.unitigs[x] for x in u.succ)
            if a.uid in seen or b.uid in seen or a.twin == b.uid:
                continue  # a sequence and its own reverse complement (inverted repeat) are not a bubble
            if a.succ and b.succ and a.succ == b.succ and abs(len(a.seq) - len(b.seq)) <= 5 and len(a.pred) == len(b.pred) == 1:
                weak, strong = (a, b) if a.cov < b.cov else (b, a)
                frac = weak.cov / (weak.cov + strong.cov)
                if frac < 0.45:
                    minor.append({"fraction": round(frac, 3), "minor": weak.seq, "major": strong.seq})
                _remove(counts, weak.seq[1:-1] if len(weak.seq) > 2 * k else weak.seq, k)
                seen.update({a.uid, b.uid, a.twin, b.twin})
                changed = True
        if not changed:
            break
        g = DeBruijn(counts, k)
    return g


# ------------------------------------------------------------------ path finding

def _base_coverage(g: DeBruijn) -> float:
    """Coverage of single-copy sequence: length-weighted median of unitig coverage."""
    vals = []
    for u in g.unitigs:
        vals.extend([u.cov] * max(1, g.nk(u)))
    return statistics.median(vals) if vals else 0.0


def _components(g: DeBruijn) -> list:
    seen, comps = set(), []
    for u in g.unitigs:
        if u.uid in seen:
            continue
        stack, comp = [u.uid], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            ux = g.unitigs[x]
            stack.extend(ux.succ + ux.pred + ([ux.twin] if ux.twin >= 0 else []))
        seen |= comp
        comps.append(comp)
    return comps


def _circuits(g: DeBruijn, comp: set, base: float, limit: int = 2000) -> list:
    """Circular walks through the graph, best coverage fit first.

    Every unitig with clear coverage (>= 0.5x single copy) must be used at least once; repeats
    may be used more often, up to one more than their coverage suggests. Among all walks the ones
    whose use-counts best match the coverage come first. Weak unitigs (< 0.5x) are optional.
    """
    pairs = {}
    for x in comp:
        t = g.unitigs[x].twin
        pairs[x] = min(x, t) if t >= 0 else x
    est, cap, required = {}, {}, set()
    for x in comp:
        key = pairs[x]
        m = g.unitigs[x].cov / base if base else 1.0
        est[key] = m
        if m >= 0.5:
            required.add(key)
            cap[key] = max(1, round(m)) + 1
        else:
            cap[key] = 0
    if not required or sum(cap.values()) > 60:
        return []
    start_key = max(required, key=lambda kk: (g.nk(g.unitigs[kk]), -est[kk]))
    found = []
    steps = [0]

    def dfs(path: list, used: dict):
        steps[0] += 1
        if steps[0] > 200000 or len(found) >= limit:
            return
        cur = path[-1]
        if path[0] in g.unitigs[cur].succ and all(used.get(r, 0) >= 1 for r in required):
            found.append((list(path), dict(used)))
        for nxt in g.unitigs[cur].succ:
            key = pairs.get(nxt)
            if key is None or used.get(key, 0) >= cap[key]:
                continue
            if nxt == path[0] and used.get(key, 0) >= 1:
                continue  # returning to the start closes the circle; handled above
            used[key] = used.get(key, 0) + 1
            path.append(nxt)
            dfs(path, used)
            path.pop()
            used[key] -= 1

    dfs([start_key], {start_key: 1})

    def misfit(item):
        _, used = item
        return sum((used.get(kk, 0) - est[kk]) ** 2 * g.nk(g.unitigs[kk]) for kk in est)

    found.sort(key=misfit)
    uniq, seen, best_fit = [], set(), None
    for item in found:
        fit = misfit(item)
        if best_fit is None:
            best_fit = fit
        if fit > best_fit * 1.5 + 1:   # keep only walks that fit (almost) as well as the best
            break
        seq = _walk_seq(g, item[0])
        if seq not in seen and revcomp(seq) not in seen:
            seen.add(seq)
            uniq.append(item[0])
    return uniq


def _walk_seq(g: DeBruijn, walk: list) -> str:
    k = g.k
    seq = g.unitigs[walk[0]].seq
    for x in walk[1:]:
        seq += g.unitigs[x].seq[k - 1:]
    # the walk returns to its start: drop the duplicated (k-1) overlap
    return seq[: len(seq) - (k - 1)]


# ------------------------------------------------------------------ polishing

def polish(seq: str, reads: list, rounds: int = 2) -> tuple:
    """Ungapped read pile-up on the circular sequence; fix positions where most reads disagree."""
    from .readcheck import pileup_reads
    fixes = 0
    for _ in range(rounds):
        depth, alt = pileup_reads(reads, seq)
        changed = False
        s = list(seq)
        for pos, counts in alt.items():
            d = depth[pos]
            if not d:
                continue
            base, n = max(counts.items(), key=lambda kv: kv[1])
            if n > 0.6 * d and base in "ACGT":
                s[pos] = base
                fixes += 1
                changed = True
        seq = "".join(s)
        if not changed:
            break
    depth, _ = pileup_reads(reads, seq)
    return seq, fixes, depth


# ------------------------------------------------------------------ main entry

def assemble(reads: Iterable, expected_size: Optional[int] = None, target_depth: int = 100,
             k_values: Optional[list] = None, log=print) -> AssemblyResult:
    """reads: iterable of (seq, qual) tuples (qual may be None)."""
    cleaned = []
    total = 0
    for seq, qual in reads:
        s = clean_read(seq, qual)
        if s:
            cleaned.append(s)
            total += len(s)
    if not cleaned:
        return AssemblyResult(None, False, None, [], ["No usable reads after cleaning."], [], {})
    lens = sorted(len(s) for s in cleaned)
    med_len = lens[len(lens) // 2]
    # subsample deterministically to ~target_depth over the expected size (or ~3 Mb if unknown)
    budget = (expected_size or 30000) * target_depth
    if total > budget:
        step = total / budget
        cleaned = [cleaned[int(i * step)] for i in range(int(len(cleaned) / step))]
    log(f"Using {len(cleaned):,} cleaned reads (median length {med_len} bp).")

    if not k_values:
        top = min(med_len - 8, 99)
        k_values = sorted({kk for kk in (top, top - 14, 61, 51, 41, 31) if 21 <= kk <= top}, reverse=True)
        k_values = [kk if kk % 2 else kk - 1 for kk in k_values]  # odd k: no palindromic k-mers
        k_values = sorted(set(k_values), reverse=True)

    best = None
    for k in k_values:
        cnt = Counter()
        for s in cleaned:
            if len(s) < k:
                continue
            cnt.update(s[i:i + k] for i in range(len(s) - k + 1))
            r = revcomp(s)
            cnt.update(r[i:i + k] for i in range(len(r) - k + 1))
        cnt.pop("", None)
        if not cnt:
            continue
        hist = sorted(v for v in cnt.values() if v >= 3)
        med = hist[len(hist) // 2] if hist else 0
        thr = max(3, int(0.1 * med))
        counts = {km: c for km, c in cnt.items() if c >= thr and "N" not in km}
        minor: list = []
        g = _clean(counts, k, minor)
        base = _base_coverage(g)
        comps = sorted(_components(g), key=lambda c: -sum(len(g.unitigs[x].seq) for x in c))
        contigs, notes = [], []
        main_seq, main_circ = None, False
        done_twins = set()
        for ci, comp in enumerate(comps):
            if any(x in done_twins for x in comp):
                continue
            done_twins |= {g.unitigs[x].twin for x in comp} | comp
            comp_cov = sum(g.unitigs[x].cov * len(g.unitigs[x].seq) for x in comp) / max(1, sum(len(g.unitigs[x].seq) for x in comp))
            circuits = _circuits(g, comp, base if ci == 0 else comp_cov)
            if circuits:
                seqs = [_walk_seq(g, w) for w in circuits]
                seq = seqs[0]
                circ = True
                if len(circuits) > 1:
                    notes.append(f"k={k}: {len(circuits)} different circular arrangements fit the reads "
                                 "(repeats longer than the reads). They are reported as alternatives.")
            else:
                longest = max(comp, key=lambda x: len(g.unitigs[x].seq))
                seq, circ, seqs = g.unitigs[longest].seq, False, None
            contigs.append({"name": f"contig_{len(contigs) + 1}", "seq": seq, "circular": circ,
                            "coverage": round(comp_cov, 1), "alternatives": seqs[1:] if seqs else []})
            if main_seq is None:
                main_seq, main_circ = seq, circ
        stats = {"k": k, "unitigs": len(g.unitigs), "base_kmer_coverage": round(base, 1),
                 "solid_threshold": thr, "reads_used": len(cleaned)}
        size_ok = (expected_size is None or main_seq is None or
                   abs(len(main_seq) - expected_size) <= 0.25 * expected_size)
        log(f"k={k}: {len(g.unitigs)} graph pieces, main contig "
            f"{len(main_seq) if main_seq else 0:,} bp, circular={main_circ}")
        result = AssemblyResult(main_seq, main_circ, k, contigs, notes, minor, stats)
        if main_circ and size_ok and not (contigs and contigs[0]["alternatives"]):
            best = result
            break
        if best is None or (main_circ and not best.circular):
            best = result
    if best and best.sequence:
        seq, fixes, depth = polish(best.sequence, cleaned)
        if fixes:
            best.notes.append(f"Polishing corrected {fixes} position(s) using the reads.")
        best.sequence = seq
        best.contigs[0]["seq"] = seq
        best.stats["polish_fixes"] = fixes
        best.stats["min_depth"] = int(min(depth)) if len(depth) else 0
        best.stats["median_depth"] = int(statistics.median(depth)) if len(depth) else 0
        if best.stats["min_depth"] == 0:
            best.notes.append("Some positions are not covered by any read: the circle may contain an error there.")
    return best


# ================================================================== long reads (Nanopore)

def assemble_long(reads: Iterable, expected_size: Optional[int] = None, max_reads: int = 60, rounds: int = 3,
                  log=print) -> AssemblyResult:
    """Whole-plasmid Nanopore data: most reads are complete (linearised) copies of the plasmid.

    1. find the plasmid size from the read-length peak;
    2. take a high-quality full-length read as the draft;
    3. align other full-length reads to it and replace every position where most reads disagree;
       repeat until nothing changes.
    """
    from . import plasmid_align
    from .inspect import estimate_circular_size
    import numpy as np

    rs = [(s.upper(), q) for s, q in reads if len(s) >= 1000]
    if not rs:
        return AssemblyResult(None, False, None, [], ["No reads of at least 1,000 bp."], [], {})
    size, _ = estimate_circular_size(np.array([len(s) for s, _ in rs]))
    size = size or expected_size
    if not size:
        return AssemblyResult(None, False, None, [], ["Could not find a plasmid-size peak in the read lengths."], [], {})
    full = [(s, q) for s, q in rs if abs(len(s) - size) <= 0.04 * size]
    if len(full) < 5:
        return AssemblyResult(None, False, None, [], [f"Only {len(full)} full-length reads (~{size:,} bp); need at least 5."], [], {})

    def mean_q(q):
        return sum(ord(c) - 33 for c in q) / len(q) if q else 0

    full.sort(key=lambda r: -mean_q(r[1]))
    draft = full[0][0]
    pool = full[1:max_reads + 1]
    log(f"Plasmid size from reads: ~{size:,} bp; {len(full)} full-length reads; polishing with {len(pool)}.")
    total_fixes = 0
    for rnd in range(rounds):
        votes: dict = {}
        used = 0
        for s, _ in pool:
            c = plasmid_align.compare(s, draft, circular=True)
            if c.ref_covered < 90:
                continue
            used += 1
            for v in c.variants:
                if v.structural:
                    continue
                key = (v.ref_pos, v.ref, v.alt)
                votes[key] = votes.get(key, 0) + 1
        accepted = [kv for kv, n in votes.items() if n > 0.5 * used]
        # apply right to left so coordinates stay valid; skip overlapping edits
        accepted.sort(key=lambda x: -x[0])
        seq, last = draft, None
        applied = 0
        for pos, ref, alt in accepted:
            if last is not None and pos + len(ref) > last:
                continue
            seq = seq[:pos] + alt + seq[pos + len(ref):]
            last = pos
            applied += 1
        log(f"Polishing round {rnd + 1}: {used} reads aligned, {applied} correction(s).")
        total_fixes += applied
        draft = seq
        if not applied:
            break
    contigs = [{"name": "contig_1", "seq": draft, "circular": True, "coverage": len(full), "alternatives": []}]
    notes = ["Nanopore consensus: positions in runs of the same letter (e.g. AAAAAA) are the least certain."]
    return AssemblyResult(draft, True, None, contigs, notes, [],
                          {"full_length_reads": len(full), "reads_used": len(pool), "polish_fixes": total_fixes,
                           "size_from_reads": size})
