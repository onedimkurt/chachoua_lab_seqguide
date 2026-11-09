"""Synthetic test/demo data: a plasmid map, a mutated 'assembly', and simulated reads."""

from __future__ import annotations

import gzip
import random
from pathlib import Path

from seqguide.effects import CODON_TABLE
from seqguide.seqio import revcomp

SENSE = [c for c, a in CODON_TABLE.items() if a != "*"]


def random_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def random_orf(n_codons: int, rng: random.Random) -> str:
    return "ATG" + "".join(rng.choice(SENSE) for _ in range(n_codons - 2)) + "TAA"


def build_plasmid(seed: int = 1) -> tuple:
    """Return (sequence, features) where features = [(type, label, start0, end, strand)]."""
    rng = random.Random(seed)
    parts, feats = [], []
    pos = 0

    def add(seq, ftype=None, label=None, strand=1):
        nonlocal pos
        if ftype:
            feats.append((ftype, label, pos, pos + len(seq), strand))
        parts.append(seq if strand == 1 else revcomp(seq))
        pos += len(seq)

    add(random_seq(200, rng))
    add(random_seq(60, rng), "promoter", "lac promoter")
    add(random_seq(20, rng))
    add(random_orf(240, rng), "CDS", "GFP")           # 720 bp
    add(random_seq(300, rng))
    add(random_seq(589, rng), "rep_origin", "ori")
    add(random_seq(250, rng))
    add(random_orf(287, rng), "CDS", "AmpR", strand=-1)  # 861 bp, minus strand
    add(random_seq(150, rng))
    add(random_seq(40, rng) + "AAAAAAA" + random_seq(40, rng))  # homopolymer region
    add(random_seq(1500, rng))
    return "".join(parts), feats


def write_genbank(path: Path, name: str, seq: str, feats: list) -> None:
    lines = [f"LOCUS       {name:<16}{len(seq):>11} bp    DNA     circular SYN 01-JAN-2026",
             "DEFINITION  synthetic test plasmid.", "FEATURES             Location/Qualifiers",
             f"     source          1..{len(seq)}", '                     /organism="synthetic"']
    for ftype, label, s, e, strand in feats:
        loc = f"{s + 1}..{e}"
        if strand == -1:
            loc = f"complement({loc})"
        lines.append(f"     {ftype:<16}{loc}")
        lines.append(f'                     /label="{label}"')
    lines.append("ORIGIN")
    low = seq.lower()
    for i in range(0, len(low), 60):
        chunk = " ".join(low[j:j + 10] for j in range(i, min(i + 60, len(low)), 10))
        lines.append(f"{i + 1:>9} {chunk}")
    lines.append("//")
    path.write_text("\n".join(lines) + "\n")


def mutate(seq: str, feats: list) -> tuple:
    """Apply known edits. Returns (mutated_seq, expected) with expected edits in reference coordinates."""
    f = {lab: (s, e, st) for _, lab, s, e, st in feats}
    gs, ge, _ = f["GFP"]
    as_, ae, _ = f["AmpR"]
    s = list(seq)
    expected = {}
    # 1) missense SNP in GFP codon 10 (first base)
    p = gs + 27
    old = s[p]
    s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[old]
    expected["snp"] = p
    # 2) 1 bp deletion inside AmpR (minus strand) -> frameshift
    pdel = as_ + 300
    expected["del"] = pdel
    # 3) 300 bp insertion in the intergenic spacer after ori
    ori_end = f["ori"][1]
    pins = ori_end + 100
    expected["ins"] = pins
    ins = random_seq(300, random.Random(99))
    assert pins < pdel
    del s[pdel]                 # apply right-most edit first so coordinates stay valid
    s[pins:pins] = list(ins)
    return "".join(s), expected


def sim_illumina_pairs(seq: str, n_pairs: int, out_prefix: Path, read_len: int = 150, insert: int = 350,
                       seed: int = 3, gz: bool = True) -> tuple:
    rng = random.Random(seed)
    circ = seq + seq[:insert]
    op = gzip.open if gz else open
    ext = ".fastq.gz" if gz else ".fastq"
    p1, p2 = Path(f"{out_prefix}_R1{ext}"), Path(f"{out_prefix}_R2{ext}")
    with op(p1, "wt") as f1, op(p2, "wt") as f2:
        for i in range(n_pairs):
            start = rng.randrange(len(seq))
            frag = circ[start:start + insert]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            r1, r2 = frag[:read_len], revcomp(frag)[:read_len]
            q = "".join(rng.choice("FFFF:F,") for _ in range(read_len))
            name = f"M00123:42:000000000-ABCDE:1:1101:{1000 + i}:{2000 + i}"
            f1.write(f"@{name} 1:N:0:1\n{r1}\n+\n{q}\n")
            f2.write(f"@{name} 2:N:0:1\n{r2}\n+\n{q}\n")
    return p1, p2


def sim_nanopore(seq: str, n: int, path: Path, seed: int = 5, error: float = 0.03) -> Path:
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for i in range(n):
            cut = rng.randrange(len(seq))
            read = seq[cut:] + seq[:cut]
            if rng.random() < 0.08:
                read = read + read  # dimer
            elif rng.random() < 0.15:
                read = read[: rng.randrange(500, len(read))]  # fragment
            if rng.random() < 0.5:
                read = revcomp(read)
            out = []
            for b in read:
                r = rng.random()
                if r < error / 3:
                    continue
                if r < 2 * error / 3:
                    out.append(rng.choice("ACGT"))
                elif r < error:
                    out.append(rng.choice("ACGT") + b)
                else:
                    out.append(b)
            read = "".join(out)
            q = "".join(chr(33 + rng.randint(8, 25)) for _ in read)
            fh.write(f"@0a1b2c3d-{i:04d}-4e5f-8a9b-0c1d2e3f4a5b runid=abc123 read={i}\n{read}\n+\n{q}\n")
    return path


def make_demo(folder: Path) -> list:
    """Write a realistic practice data set into `folder` (used by the 'Load demo data' button).

    The reads come from a plasmid that carries three known changes relative to the map, so a
    correct analysis should report: a missense change in GFP, a frameshift in AmpR and a
    300 bp insertion between ori and AmpR.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    seq, feats = build_plasmid()
    write_genbank(folder / "pDemo_expected_map.gb", "pDemo", seq, feats)
    mut, _ = mutate(seq, feats)
    r1, r2 = sim_illumina_pairs(mut, 20000, folder / "pDemo_illumina")
    np_ = sim_nanopore(mut, 600, folder / "pDemo_nanopore.fastq")
    return [folder / "pDemo_expected_map.gb", r1, r2, np_]
