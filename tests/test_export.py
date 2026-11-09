import gzip
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from seqguide import demo, effects, export, plasmid_align, refmap  # noqa: E402
from seqguide.readcheck import map_and_pileup  # noqa: E402
from seqguide.seqio import read_fastq, revcomp  # noqa: E402


def _read_bam(path):
    """Minimal BAM reader (BGZF is multi-member gzip) used to check what we wrote."""
    data = gzip.decompress(Path(path).read_bytes())
    assert data[:4] == b"BAM\x01"
    l_text = struct.unpack_from("<i", data, 4)[0]
    off = 8 + l_text
    n_ref = struct.unpack_from("<i", data, off)[0]
    off += 4
    refs = []
    for _ in range(n_ref):
        ln = struct.unpack_from("<i", data, off)[0]
        name = data[off + 4: off + 4 + ln - 1].decode()
        L = struct.unpack_from("<i", data, off + 4 + ln)[0]
        refs.append((name, L))
        off += 8 + ln
    recs = []
    while off < len(data):
        bs = struct.unpack_from("<i", data, off)[0]
        ref_id, pos, l_rn, mapq, bin_, n_cig, flag, l_seq = struct.unpack_from("<iiBBHHHi", data, off + 4)
        p = off + 4 + 32
        name = data[p:p + l_rn - 1].decode()
        p += l_rn
        cig = [struct.unpack_from("<I", data, p + 4 * i)[0] for i in range(n_cig)]
        p += 4 * n_cig
        seq = "".join("=ACMGRSVTWYHKDBN"[(data[p + i // 2] >> (4 if i % 2 == 0 else 0)) & 15] for i in range(l_seq))
        recs.append(dict(ref_id=ref_id, pos=pos, name=name, mapq=mapq, bin=bin_, flag=flag, seq=seq,
                         cigar=[(c >> 4, "MIDNSHP=X"[c & 15]) for c in cig]))
        off += 4 + bs
    return refs, recs


def test_reg2bin_matches_spec_examples():
    assert export.reg2bin(0, 1) == 4681
    assert export.reg2bin(16383, 16385) == 585          # crosses a 16 kb boundary -> level-4 bin
    assert export.reg2bin(0, 1 << 29) == 0


def test_bam_roundtrip_and_index(tmp_path):
    seq, feats = demo.build_plasmid()
    r1, r2 = demo.sim_illumina_pairs(seq, 1500, tmp_path / "S", gz=False)
    placed = []
    map_and_pileup(((r.seq, r.qual, r.id) for f in (r1, r2) for r in read_fastq(f)), seq, alignments=placed)
    aln = export.alignments_from_readcheck(placed, len(seq))
    n = export.write_bam(tmp_path / "x.bam", [("pT", len(seq))], aln)
    refs, recs = _read_bam(tmp_path / "x.bam")
    assert refs == [("pT", len(seq))] and len(recs) == n > 2500
    assert [r["pos"] for r in recs] == sorted(r["pos"] for r in recs)
    for r in recs[:500]:
        m = sum(k for k, op in r["cigar"] if op == "M")
        lead = r["cigar"][0][0] if r["cigar"][0][1] == "S" else 0
        aligned = r["seq"][lead:lead + m]
        assert aligned == seq[r["pos"]:r["pos"] + m]               # M part matches the reference exactly
        assert sum(k for k, op in r["cigar"] if op in "MIS") == len(r["seq"])
        assert r["pos"] + m <= len(seq)
        assert r["bin"] == export.reg2bin(r["pos"], r["pos"] + m)
    bai = (tmp_path / "x.bam.bai").read_bytes()
    assert bai[:4] == b"BAI\x01" and struct.unpack_from("<i", bai, 4)[0] == 1
    assert (tmp_path / "x.bam").read_bytes().endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))


def test_bundle_carries_features_to_the_sample(tmp_path):
    seq, feats = demo.build_plasmid()
    demo.write_genbank(tmp_path / "ref.gb", "pT", seq, feats)
    ref = refmap.load_references(tmp_path / "ref.gb")[0]
    mut, _ = demo.mutate(seq, feats)
    cmp = plasmid_align.compare(revcomp(mut), ref.seq)
    effects.annotate(ref, cmp.variants)
    moved = export.transfer_features(ref.features, cmp.variants, len(cmp.rotated_seq))
    by = {f.label: f for f in moved}
    # GFP keeps its sequence except the SNP; AmpR is one base shorter (frameshift deletion)
    g_ref = next(f for f in ref.features if f.label == "GFP")
    g_new = by["GFP"]
    a = ref.seq[g_ref.start:g_ref.end]
    b = cmp.rotated_seq[g_new.start:g_new.end]
    assert len(a) == len(b) and sum(x != y for x, y in zip(a, b)) == 1
    assert by["AmpR"].length() == next(f for f in ref.features if f.label == "AmpR").length() - 1
    assert by["ori"].start == next(f for f in ref.features if f.label == "ori").start   # before the insertion
    out = export.make_bundle(tmp_path / "gb", "demo sample",
                             [{"name": "pT", "seq": cmp.rotated_seq, "features": moved, "variants": cmp.variants,
                               "circular": True}])
    back = refmap.load_references(next(out.glob("*.gb")))[0]
    assert back.seq == cmp.rotated_seq and back.circular
    labels = {f.label for f in back.features}
    assert {"GFP", "AmpR", "ori", "lac promoter"} <= labels
    assert sum(1 for f in back.features if f.type == "misc_difference") == 3
    fai = next(out.glob("*.fasta.fai")).read_text().split("\t")
    assert int(fai[1]) == len(cmp.rotated_seq)
    assert (out / "igv_session.xml").exists() and (out / "differences.bed").exists()
