#!/usr/bin/env python3
"""Stand-in for bioinformatics tools so the orchestration can be tested without them."""
import os, shutil, sys
from pathlib import Path

tool = Path(sys.argv[0]).name
args = sys.argv[1:]

def opt(flag, default=None):
    return args[args.index(flag) + 1] if flag in args else default

def fasta_names(path):
    out = []
    name, n = None, 0
    for line in open(path):
        if line.startswith(">"):
            if name: out.append((name, n))
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name: out.append((name, n))
    return out

if "--version" in args or args[:1] == ["version"]:
    print(f"{tool} 0.0-fake"); sys.exit(0)

if tool == "fastp":
    shutil.copy(opt("-i"), opt("-o"))
    if opt("-I"): shutil.copy(opt("-I"), opt("-O"))
    Path(opt("-j")).write_text("{}"); Path(opt("-h")).write_text("<html></html>")
elif tool == "seqkit":
    shutil.copy(args[-1], opt("-o"))
elif tool == "unicycler":
    out = Path(opt("-o")); out.mkdir(exist_ok=True)
    seq = Path(os.environ["FAKE_ASSEMBLY"]).read_text().strip()
    (out / "assembly.fasta").write_text(f">1 length={len(seq)} depth=1.00x circular=true\n{seq}\n>2 length=640 depth=0.05x\n{'ACGT'*160}\n")
elif tool == "minimap2":
    sys.stdout.write("@HD\tVN:1.6\n")
elif tool == "samtools":
    sub = args[0]
    if sub == "sort":
        sys.stdin.read(); Path(opt("-o")).write_text("BAM")
    elif sub in ("index", "faidx"):
        pass
    elif sub == "depth":
        fa = "final_assembly.fasta" if Path("final_assembly.fasta").exists() else "ref.fa"
        for name, n in fasta_names(fa):
            for i in range(1, n + 1):
                d = 150 if not (200 <= i < 260) else 12
                print(f"{name}\t{i}\t{d}")
    elif sub == "flagstat":
        print("10000 + 0 in total (QC-passed reads + QC-failed reads)\n9500 + 0 mapped (95.00% : N/A)")
elif tool == "bcftools":
    if args[0] == "norm":
        sys.stdin.read(); sys.stdout.write(Path(os.environ["FAKE_VCF"]).read_text())
    elif args[0] == "call":
        sys.stdout.write(sys.stdin.read())
    else:  # mpileup
        sys.stdout.write("##fileformat=VCF\n")
sys.exit(0)
