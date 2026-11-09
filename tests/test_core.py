import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from seqguide import decide, effects, inspect, plasmid_align, refmap  # noqa: E402
from seqguide.seqio import revcomp  # noqa: E402
from seqguide import demo as simulate  # noqa: E402


@pytest.fixture(scope="module")
def plasmid(tmp_path_factory):
    d = tmp_path_factory.mktemp("plasmid")
    seq, feats = simulate.build_plasmid()
    gb = d / "pTest.gb"
    simulate.write_genbank(gb, "pTest", seq, feats)
    return seq, feats, gb, d


# ---------------------------------------------------------------- references

def test_genbank_parse(plasmid):
    seq, feats, gb, _ = plasmid
    ref = refmap.load_references(gb)[0]
    assert ref.seq == seq and ref.circular
    labels = {f.label: f for f in ref.features}
    assert set(labels) == {"lac promoter", "GFP", "ori", "AmpR"}
    assert labels["AmpR"].strand == -1
    s, e = next((f[2], f[3]) for f in feats if f[1] == "AmpR")
    assert (labels["AmpR"].start, labels["AmpR"].end) == (s, e)
    # AmpR translates to a full ORF on the minus strand
    prot = effects.translate(effects._cds_seq(ref.seq, labels["AmpR"], effects._cds_positions(labels["AmpR"])))
    assert prot.startswith("M") and prot.endswith("*") and prot.count("*") == 1


def test_location_parser():
    assert refmap.parse_location("complement(join(1..10,21..30))") == [(20, 30, -1), (0, 10, -1)]
    assert refmap.parse_location("join(5000..5386,1..200)") == [(4999, 5386, 1), (0, 200, 1)]
    assert refmap.parse_location("<1..>50") == [(0, 50, 1)]


# ---------------------------------------------------------------- comparison

def test_identical_rotated_revcomp(plasmid):
    seq, _, _, _ = plasmid
    asm = revcomp(seq[1234:] + seq[:1234])
    cmp = plasmid_align.compare(asm, seq)
    assert cmp.variants == []
    assert cmp.strand == "-"
    assert cmp.identity == 100.0
    assert cmp.rotated_seq == seq


def test_detects_known_edits(plasmid):
    seq, feats, gb, _ = plasmid
    mut, exp = simulate.mutate(seq, feats)
    asm = revcomp(mut[2500:] + mut[:2500])
    ref = refmap.load_references(gb)[0]
    cmp = plasmid_align.compare(asm, ref.seq)
    effects.annotate(ref, cmp.variants)
    kinds = {v.kind: v for v in cmp.variants}
    assert set(kinds) == {"SNP", "DEL", "INS"}, [v.describe() for v in cmp.variants]
    assert kinds["SNP"].ref_pos == exp["snp"]
    assert kinds["INS"].size == 300
    assert abs(kinds["DEL"].ref_pos - exp["del"]) <= 6  # left-normalisation may shift within a repeat
    assert kinds["SNP"].severity in ("moderate", "high") and "GFP" in kinds["SNP"].effect
    assert "Frameshift" in kinds["DEL"].effect and "AmpR" in kinds["DEL"].effect
    assert kinds["INS"].severity == "low"  # intergenic
    v = effects.verdict(cmp.variants, cmp.ref_covered, cmp.identity, cmp.notes)
    assert v["level"] == "block"


def test_start_region_deleted(plasmid):
    seq, _, _, _ = plasmid
    asm = seq[150:]  # first 150 bp missing
    cmp = plasmid_align.compare(asm[700:] + asm[:700], seq)
    assert len(cmp.variants) == 1 and cmp.variants[0].kind == "DEL" and len(cmp.variants[0].ref) == 150


def test_stop_gain_and_homopolymer(plasmid):
    seq, feats, gb, _ = plasmid
    ref = refmap.load_references(gb)[0]
    gfp = next(f for f in ref.features if f.label == "GFP")
    # turn codon 5 into TAA
    p = gfp.start + 12
    mut = seq[:p] + "TAA" + seq[p + 3:]
    # homopolymer 1-bp deletion
    hp = seq.index("AAAAAAA")
    mut = mut[:hp] + mut[hp + 1:]
    cmp = plasmid_align.compare(mut, seq)
    effects.annotate(ref, cmp.variants, platform="nanopore")
    stop = [v for v in cmp.variants if "stop" in v.effect]
    assert stop and stop[0].severity == "high"
    hpv = [v for v in cmp.variants if v.kind == "DEL"]
    assert hpv and hpv[0].homopolymer == 7 and hpv[0].flags


def test_inversion_reported(plasmid):
    seq, _, _, _ = plasmid
    mut = seq[:3000] + revcomp(seq[3000:4000]) + seq[4000:]
    cmp = plasmid_align.compare(mut, seq)
    assert any(v.kind == "INVERSION" for v in cmp.variants), [v.describe() for v in cmp.variants]


def test_unrelated_sequence():
    rng = random.Random(7)
    a = simulate.random_seq(3000, rng)
    b = simulate.random_seq(3000, rng)
    cmp = plasmid_align.compare(a, b)
    assert cmp.notes and cmp.ref_covered < 20


# ---------------------------------------------------------------- inspection

def test_inspect_illumina_pairs(plasmid, tmp_path):
    seq, *_ = plasmid
    r1, r2 = simulate.sim_illumina_pairs(seq, 2000, tmp_path / "S1", gz=True)
    ds = inspect.inspect_dataset([r1, r2])
    assert len(ds.samples) == 1
    s = ds.samples[0]
    assert s.sample == "S1" and s.layout == "paired" and s.pair_check["ok"]
    f = s.files[0]
    assert f.platform == "illumina" and f.n_reads == 2000 and f.max_len == 150 and f.compressed
    assert f.mean_q > 25
    rec = decide.recommend(s, "plasmid", "annotated", "verify", expected_size=len(seq))
    assert rec.workflow == "plasmid_denovo" and not rec.blocked
    assert rec.params["assembler"] == "unicycler"
    assert 100 < rec.params["depth_short"] < 140  # 2000*2*150/5600 ≈ 107


def test_mismatched_pairs_blocked(plasmid, tmp_path):
    seq, *_ = plasmid
    a1, _ = simulate.sim_illumina_pairs(seq, 300, tmp_path / "A", seed=1, gz=False)
    _, b2 = simulate.sim_illumina_pairs(seq, 280, tmp_path / "B", seed=2, gz=False)
    b2.rename(tmp_path / "A_R2.fastq")
    ds = inspect.inspect_dataset([a1, tmp_path / "A_R2.fastq"])
    s = ds.samples[0]
    assert not s.pair_check["ok"]
    rec = decide.recommend(s, "plasmid", "annotated", "verify", expected_size=len(seq))
    assert rec.blocked


def test_inspect_nanopore_size_and_dimer(plasmid, tmp_path):
    seq, *_ = plasmid
    p = simulate.sim_nanopore(seq, 400, tmp_path / "np.fastq")
    rep = inspect.inspect_file(p)
    assert rep.platform == "nanopore"
    assert abs(rep.plasmid_size_estimate - len(seq)) / len(seq) < 0.05
    assert rep.dimer_peak
    ds = inspect.inspect_dataset([p])
    rec = decide.recommend(ds.samples[0], "plasmid", "annotated", "verify", expected_size=len(seq) + 2000)
    assert any("different" in m.text for m in rec.messages)  # size mismatch warning
    assert rec.params["assembler"] == "flye+medaka"


def test_amplicon_needs_reference(plasmid, tmp_path):
    seq, *_ = plasmid
    r1, r2 = simulate.sim_illumina_pairs(seq[:500], 500, tmp_path / "amp", gz=False)
    s = inspect.inspect_dataset([r1, r2]).samples[0]
    assert decide.recommend(s, "amplicon", "none", "find_mutations").blocked
    rec = decide.recommend(s, "amplicon", "sequence", "find_mutations", expected_size=500)
    assert rec.workflow == "amplicon_align" and rec.params["preset"] == "sr"


def test_reference_files_classified(plasmid, tmp_path):
    _, _, gb, _ = plasmid
    ds = inspect.inspect_dataset([gb])
    assert ds.references and ds.references[0].format == "genbank"


def test_spades_end_overlap_is_not_an_insertion(plasmid):
    seq, *_ = plasmid
    contig = seq[1500:] + seq[:1500]
    contig = contig + contig[:77]          # SPAdes-style k-mer overlap at the ends
    cmp = plasmid_align.compare(revcomp(contig), seq)
    assert cmp.variants == [] and cmp.trimmed_overlap == 77 and cmp.notes


def test_illumina_sample_names():
    assert inspect.sample_name("BA137_S1_L001_R1_001.fastq.gz") == "BA137"
    assert inspect.sample_name("pX_S12_R2.fastq") == "pX"


def test_readcheck_finds_snp_and_depth(plasmid, tmp_path):
    from seqguide.readcheck import check_reads
    seq, feats, *_ = plasmid
    mut = seq[:1000] + ("A" if seq[1000] != "A" else "C") + seq[1001:]
    r1, _ = simulate.sim_illumina_pairs(mut, 3000, tmp_path / "rc", gz=False)
    rc = check_reads(r1, seq)
    assert rc.pct_mapped > 95
    v = rc.variants()
    assert len(v) == 1 and v[0]["pos"] == 1001 and v[0]["af"] > 0.95
    assert rc.depth.min() > 20  # including across the origin of the circle


def test_plasmid_figure_draws_depth_panel(plasmid):
    import numpy as np
    from seqguide import plots
    seq, _, gb, _ = plasmid
    ref = refmap.load_references(gb)[0]
    import matplotlib.pyplot as plt
    called = {}
    orig = plots._depth_axis
    plots._depth_axis = lambda *a, **k: called.setdefault("yes", True) or orig(*a, **k)
    try:
        plots.plasmid_figure(ref, [], np.full(len(seq), 100))
    finally:
        plots._depth_axis = orig
    assert called.get("yes")
