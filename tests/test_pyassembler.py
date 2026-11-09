import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from seqguide import demo, plasmid_align, pyassembler  # noqa: E402
from seqguide.seqio import read_fastq, revcomp  # noqa: E402


def _reads(seq, n_pairs, tmp_path, name="r", seed=3):
    r1, r2 = demo.sim_illumina_pairs(seq, n_pairs, tmp_path / name, seed=seed, gz=False)
    return [(r.seq, r.qual) for f in (r1, r2) for r in read_fastq(f)]


def _same_circle(a, b):
    return len(a) == len(b) and (b in a + a or revcomp(b) in a + a)


def test_reconstructs_mutated_plasmid_exactly(tmp_path):
    seq, feats = demo.build_plasmid()
    mut, _ = demo.mutate(seq, feats)          # SNP, 1-bp deletion, 300-bp insertion
    res = pyassembler.assemble(_reads(mut, 4000, tmp_path), log=lambda m: None)
    assert res.circular and _same_circle(res.sequence, mut)
    kinds = sorted(v.kind for v in plasmid_align.compare(res.sequence, seq).variants)
    assert kinds == ["DEL", "INS", "SNP"]


def test_without_expected_size(tmp_path):
    seq, _ = demo.build_plasmid(seed=5)
    res = pyassembler.assemble(_reads(seq, 3000, tmp_path), expected_size=None, log=lambda m: None)
    assert res.circular and _same_circle(res.sequence, seq)


def test_low_depth_still_circular(tmp_path):
    seq, _ = demo.build_plasmid(seed=7)
    res = pyassembler.assemble(_reads(seq, 400, tmp_path), expected_size=len(seq), log=lambda m: None)  # ~25x
    assert res.circular and _same_circle(res.sequence, seq)


def test_mixed_population_reports_minor_variant(tmp_path):
    seq, _ = demo.build_plasmid(seed=9)
    p = 2000
    minor = seq[:p] + ("A" if seq[p] != "A" else "C") + seq[p + 1:]
    reads = _reads(seq, 2800, tmp_path, "maj", seed=1) + _reads(minor, 1200, tmp_path, "min", seed=2)
    res = pyassembler.assemble(reads, expected_size=len(seq), log=lambda m: None)
    assert res.circular and _same_circle(res.sequence, seq)      # majority sequence
    assert res.minor_variants and 0.15 < res.minor_variants[0]["fraction"] < 0.45


def test_direct_repeat_longer_than_reads(tmp_path):
    rng = random.Random(11)
    rep = demo.random_seq(600, rng)                               # like the two LTRs
    seq = demo.random_seq(2500, rng) + rep + demo.random_seq(2000, rng) + rep + demo.random_seq(1200, rng)
    res = pyassembler.assemble(_reads(seq, 4000, tmp_path), expected_size=len(seq), log=lambda m: None)
    assert res.circular and _same_circle(res.sequence, seq)


def test_inverted_repeat_is_flagged_not_hidden(tmp_path):
    rng = random.Random(13)
    rep = demo.random_seq(600, rng)
    seq = demo.random_seq(2500, rng) + rep + demo.random_seq(2000, rng) + revcomp(rep) + demo.random_seq(1200, rng)
    res = pyassembler.assemble(_reads(seq, 4000, tmp_path), expected_size=len(seq), log=lambda m: None)
    assert res.circular and len(res.sequence) == len(seq)
    alts = res.contigs[0]["alternatives"]
    assert _same_circle(res.sequence, seq) or any(_same_circle(a, seq) for a in alts)
    assert alts and any("arrangements" in n for n in res.notes)


def test_nanopore_draft_and_polish(tmp_path):
    seq, feats = demo.build_plasmid(seed=21)
    mut, _ = demo.mutate(seq, feats)
    fq = demo.sim_nanopore(mut, 60, tmp_path / "np.fastq", seed=4, error=0.06)   # noisy reads, ~50 full length
    res = pyassembler.assemble_long([(r.seq, r.qual) for r in read_fastq(fq)], log=lambda m: None)
    assert res.circular and _same_circle(res.sequence, mut)
