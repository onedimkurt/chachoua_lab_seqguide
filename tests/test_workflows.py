"""End-to-end orchestration tests using fake tools (no bioinformatics software needed)."""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from seqguide import demo as simulate  # noqa: E402
from seqguide import decide, inspect, workflows  # noqa: E402
from seqguide.seqio import revcomp  # noqa: E402

FAKE = Path(__file__).resolve().parent / "fake_tools"


def _env(monkeypatch, **kw):
    monkeypatch.setenv("PATH", f"{FAKE}{os.pathsep}{os.environ['PATH']}")
    for k, v in kw.items():
        monkeypatch.setenv(k, str(v))


def test_plasmid_denovo_end_to_end(tmp_path, monkeypatch):
    seq, feats = simulate.build_plasmid()
    gb = tmp_path / "pTest.gb"
    simulate.write_genbank(gb, "pTest", seq, feats)
    mut, _ = simulate.mutate(seq, feats)
    asm = tmp_path / "asm.txt"
    asm.write_text(revcomp(mut[1000:] + mut[:1000]))
    r1, r2 = simulate.sim_illumina_pairs(seq, 3000, tmp_path / "S1")
    _env(monkeypatch, FAKE_ASSEMBLY=asm)

    sample = inspect.inspect_dataset([r1, r2]).samples[0]
    rec = decide.recommend(sample, "plasmid", "annotated", "verify", expected_size=len(seq))
    cfg = workflows.RunConfig("plasmid_denovo", "S1", tmp_path / "run", short_reads=[str(r1), str(r2)],
                              reference=str(gb), params=rec.params, inspection=sample.to_dict(),
                              recommendation=rec.to_dict(), threads=2)
    res = workflows.run(cfg)
    assert res["error"] is None, res["error"]
    cmp = res["comparison"]
    assert {v["kind"] for v in cmp["variants"]} == {"SNP", "DEL", "INS"}
    assert cmp["verdict"]["level"] == "block"
    assert all(v["depth"] == 150 or v["depth"] == 12 for v in cmp["variants"])
    out = tmp_path / "run"
    for f in ("report.html", "differences.csv", "assembly_matched_to_reference.fasta", "provenance.json", "results.json"):
        assert (out / f).exists(), f
    html = (out / "report.html").read_text()
    assert "Frameshift" in html and "data:image/png;base64" in html
    assert json.loads((out / "provenance.json").read_text())["tool_versions"]["unicycler"].endswith("fake")
    # rotated assembly equals the mutated plasmid in reference orientation
    rot = "".join((out / "assembly_matched_to_reference.fasta").read_text().splitlines()[1:])
    assert rot == mut


def test_amplicon_end_to_end(tmp_path, monkeypatch):
    seq, feats = simulate.build_plasmid()
    gb = tmp_path / "targets.gb"
    simulate.write_genbank(gb, "GFPamp", seq, feats)
    gfp = next(f for f in feats if f[1] == "GFP")
    p = gfp[2] + 27  # same codon as simulate.mutate SNP
    alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[seq[p]]
    vcf = tmp_path / "fake.vcf"
    vcf.write_text("##fileformat=VCF\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                   f"GFPamp\t{p + 1}\t.\t{seq[p]}\t{alt}\t220\t.\tDP=300\tGT:AD:DP\t0/1:160,140:300\n"
                   f"GFPamp\t{p + 100}\t.\t{seq[p + 99]}{seq[p + 100]}\t{seq[p + 99]}\t50\t.\tDP=40\tGT:AD:DP\t1/1:2,38:40\n")
    r1, r2 = simulate.sim_illumina_pairs(seq, 1000, tmp_path / "A1")
    _env(monkeypatch, FAKE_VCF=vcf)
    sample = inspect.inspect_dataset([r1, r2]).samples[0]
    rec = decide.recommend(sample, "amplicon", "annotated", "find_mutations", expected_size=len(seq))
    cfg = workflows.RunConfig(rec.workflow, "A1", tmp_path / "run", short_reads=[str(r1), str(r2)], reference=str(gb),
                              params=rec.params, inspection=sample.to_dict(), recommendation=rec.to_dict(), threads=2)
    res = workflows.run(cfg)
    assert res["error"] is None, res["error"]
    t = res["amplicon"]["targets"][0]
    assert len(t["variants"]) == 2
    snp, dl = sorted(t["variants"], key=lambda v: v["position"])
    assert "GFP" in snp["effect"] and any("Heterozygous" in f for f in snp["flags"])
    assert "Frameshift" in dl["effect"] and any("Low depth" in f for f in dl["flags"])
    assert res["amplicon"]["mapping"]["mapped_pct"] == 95.0
    assert (tmp_path / "run" / "report.html").exists()


def test_missing_tools_message(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", "/nonexistent")
    cfg = workflows.RunConfig("plasmid_denovo", "x", tmp_path, short_reads=["a", "b"],
                              params={"assembler": "unicycler", "engine": "tools"})
    try:
        workflows.run(cfg)
    except RuntimeError as e:
        assert "Docker" in str(e)
    else:
        raise AssertionError("expected RuntimeError")


def test_read_support_uses_all_reads_not_subsample(tmp_path):
    cfg = workflows.RunConfig("plasmid_denovo", "x", tmp_path, short_reads=["r1.fq.gz", "r2.fq.gz"], reference="ref.gb",
                              params={"assembler": "unicycler", "depth_short": 2000, "subsample_depth": 150})
    steps = workflows.build_steps(cfg, {})
    asm = next(s for s in steps if s.name.startswith("De novo"))
    mapb = next(s for s in steps if s.name.startswith("Map reads back"))
    assert any("sub_" in a for a in asm.cmds[0])                 # assembler gets the subsample
    assert not any("sub_" in a for a in mapb.cmds[0])            # read support gets all trimmed reads
    assert "trim_R1.fastq.gz" in mapb.cmds[0]


def test_falls_back_to_builtin_engine_without_tools(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", "/nonexistent")
    seq, feats = simulate.build_plasmid()
    gb = tmp_path / "ref.gb"
    simulate.write_genbank(gb, "pT", seq, feats)
    mut, _ = simulate.mutate(seq, feats)
    r1, r2 = simulate.sim_illumina_pairs(mut, 3000, tmp_path / "S")
    cfg = workflows.RunConfig("plasmid_denovo", "S", tmp_path / "run", short_reads=[str(r1), str(r2)], reference=str(gb),
                              params={"expected_size": len(seq)})
    res = workflows.run(cfg)
    assert res["error"] is None and res["params"]["engine"] == "builtin"
    assert sorted(v["kind"] for v in res["comparison"]["variants"]) == ["DEL", "INS", "SNP"]
    assert (tmp_path / "run" / "report.html").exists()
