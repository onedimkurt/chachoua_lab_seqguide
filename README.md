# SeqGuide (Chachoua lab)

**Use it in your browser: https://onedimkurt.github.io/chachoua_lab_seqguide/** (Chrome, Edge or Safari; nothing to install,
your files never leave your computer).

A guided sequencing-analysis pipeline for the lab. It runs on a normal laptop (Mac or Windows) and is meant for people who are not bioinformaticians.

SeqGuide walks the user through five steps in the browser. At every step it explains **what** it found or is about to do, and **why**.

1. **Choose data:** pick FASTQ reads and, if available, the expected map (GenBank `.gb`, SnapGene `.dna` or FASTA).
2. **Understand data:** file type, Illumina vs Nanopore, paired-end check (do R1/R2 really belong together?), read count, length, quality, plasmid size from read lengths, dimers.
3. **Decide:** a few plain questions (what did you sequence, what is your goal, do you have a reference). The app then recommends an analysis with reasons and warnings (e.g. depth too low, size mismatch). It also shows the alternatives and their trade-offs.
4. **Run:** shows a progress bar and explains each step while it runs.
5. **Results:** a self-contained HTML report with a bottom-line verdict, a plasmid map showing the differences, read depth, a table of differences with protein-level effects, a glossary, and the exact tool versions and commands for your methods section.

There are two ways to run it:

* **In the browser** (link above): zero installation. Uses SeqGuide's built-in engine: its own plasmid assembler
  (validated base-for-base against Unicycler on 8 APEX2 plasmids), read check and variant finder.
* **Docker version** (below): the same app with the standard tools (fastp, Unicycler/SPAdes, Flye/Medaka, minimap2,
  samtools, bcftools). Use it for large data sets, amplicon insertions/deletions, or as an independent check.

Scope: **plasmids** (de novo assembly + comparison with the expected map) and **amplicons / targeted genes** (alignment + variant calling). Whole genomes and RNA-seq are out of scope for now.

---

## Installing (once per laptop)

1. Install **Docker Desktop** (https://www.docker.com/products/docker-desktop/) and start it.
   *Windows:* accept the WSL 2 option when the installer asks. *Apple Silicon Macs:* in Docker Desktop → Settings → General, turn on "Use Rosetta for x86/amd64 emulation".
   Docker Desktop is free for education and small organisations; check your institution's terms.
2. Copy this `seqguide` folder anywhere on the laptop.
3. Double-click **`start-mac.command`** (Mac) or **`start-windows.bat`** (Windows).
   The **first start takes 20–40 minutes** because it downloads and installs all tools. Later starts take seconds.
   The browser opens at http://localhost:8501.

Recommended laptop: 8 GB RAM minimum (16 GB is better), about 15 GB free disk for the tools.

## Daily use

* Put your files in **`SeqGuide/input`** in your home folder. Don't point SeqGuide at a folder on the Desktop or in Documents/Downloads: macOS blocks Docker from those folders unless you allow it in System Settings → Privacy & Security → Files and Folders. It is created on first start: `~/SeqGuide/input` on Mac, `C:\Users\<you>\SeqGuide\input` on Windows.
* Start SeqGuide and follow the steps.
* Results go to **`SeqGuide/runs/<sample>_<date-time>/`**. Open `report.html` in any browser to view the report, share it or print it.
* Stop SeqGuide with `stop-mac.command` / `stop-windows.bat`.
* New users: click **"Load practice data"** in step 1. It creates a demo plasmid with three hidden changes (a missense mutation in GFP, a frameshift in AmpR and a 300 bp insertion). A correct run finds all three.

---

## How it works

```
 FASTQ files ──► 1 Inspect ──► 2 Decide ──┬─► Plasmid: trim → subsample → assemble → (polish) → map reads back → compare to map → report
 (+ map)         (pure Python)  (rules)    └─► Amplicon: trim → align to targets → call variants → interpret → report
```

| Layer | What | Where |
|---|---|---|
| GUI | Streamlit wizard in the browser | `app/streamlit_app.py` |
| Inspection | format, platform, pairing, stats, plasmid size from read lengths | `seqguide/inspect.py` |
| Decision rules | recommendation + reasons + warnings; thresholds at the top of the file | `seqguide/decide.py` |
| Explanations | every "what / why" text shown to users | `seqguide/explain.py` |
| Workflows | builds and runs tool commands; records provenance | `seqguide/workflows.py`, `runner.py` |
| Plasmid comparison | orientation, rotation, alignment, structural changes | `seqguide/plasmid_align.py` |
| Effects | feature hit, codon/protein change, frameshift, stop gain, homopolymer flag | `seqguide/effects.py` |
| References | GenBank / SnapGene / FASTA reader | `seqguide/refmap.py` |
| Report | HTML report + figures | `seqguide/report.py`, `plots.py` |

### Tools used

| Task | Illumina | Nanopore |
|---|---|---|
| Trimming / QC | fastp | length filter (seqkit) |
| Subsampling to ~150× | seqkit | seqkit |
| Plasmid assembly | Unicycler (SPAdes) | Flye + Medaka polishing |
| Hybrid assembly | Unicycler with both read types | |
| Read alignment | minimap2 (`sr`) | minimap2 (`map-ont`) |
| Depth | samtools | samtools |
| Variant calling (amplicons) | bcftools (haploid or diploid, chosen by user) | bcftools (flagged as "confirm") |
| Comparison with the plasmid map | built-in, circular-aware | built-in |

The plasmid comparison is written in plain Python so it can be tested and explained. It:

1. finds the assembly's orientation and rotates it to start where your map starts (plasmids are circular, so assemblers start them anywhere);
2. aligns the two and lists SNPs, insertions, deletions, inversions and replaced regions;
3. places each difference in its feature and predicts the effect, e.g. "Frameshift in AmpR at amino acid 154" or "L10P (Leu → Pro)";
4. flags uncertain calls: low read depth, and 1-letter indels in runs like `AAAAAA` with Nanopore data.

It also writes `assembly_matched_to_reference.fasta`, the assembly in the same orientation and start position as your map, ready to open in SnapGene or Benchling.

### Built-in guidance rules (edit in `decide.py`)

* Plasmid → de novo assembly (with comparison when a map is given). Alignment-only is offered as an alternative with its drawbacks.
* Amplicon/targeted gene → alignment. Blocked without a reference, with a hint on where to get one.
* Paired files whose read names don't match → blocked.
* Depth below 30× (Illumina) / 20× (Nanopore) for plasmids, or 100× for amplicons → warning. Depth above 450× → automatic subsampling, explained to the user.
* Nanopore plasmid size from the read-length peak differs > 10% from the expected size → warning.
* Low average quality, very few reads or truncated files → warning.

---

## Validating on a new computer

`validation/` contains a validation set (8 APEX2 plasmids with known results, plus synthetic amplicon and
Nanopore samples) and one-click scripts that build SeqGuide, check every tool, run all samples and start the GUI:
`validate-mac.command` (Mac) and `validate-windows.bat` (Windows). Results go to `validation/results/`, the full
log to `validation/validation.log`. Expected results: no differences except C5196G (ori) in BA137 and BA137_align;
AMP1 G351T heterozygous (APEX2 V101F) and C552G homozygous (APEX2 H168D); NP1 C5196G and A1072C.

## For developers

```bash
pip install -r requirements-dev.txt
pytest tests                    # 20 tests, no bioinformatics tools needed
python tests/gui_smoke.py       # drives the GUI script through all five steps
streamlit run app/streamlit_app.py   # runs the GUI natively (tools must be on PATH)
```

`tests/fake_tools/` contains stand-ins for fastp, Unicycler, minimap2, samtools and bcftools. With them the whole workflow, including the report, can be tested on any machine.

Running without Docker on a Mac is possible with conda/mamba: `mamba env create -f environment.yml -n seqguide`, plus Medaka separately. Docker is the supported route because it behaves the same on Mac and Windows.

## Status and next steps

Version 0.1 is a working prototype.

**First real-data check (APEX2 plasmids, 8 samples, Illumina single-end, SPAdes assemblies):** the comparison
step agreed with the manual analysis on 8/8 plasmids, including the single real difference (C5196G in the ori of
BA137). That check also exposed one bug, now fixed: SPAdes leaves a k-mer-long (77 bp) duplicate at the ends of
circular contigs, and SeqGuide first reported it as an insertion. It now removes this overlap and says so in the report.
SeqGuide can also compare an **existing assembly** (your own, or a sequencing service's) with a map, with or without reads.
With reads, a built-in read check (`seqguide/readcheck.py`, pure Python) places every read on the map, measures depth and
confirms each difference. On the APEX2 FASTQ files it reproduced the manual bwa/bcftools results: same read-level variants
(C5196G in 728/728 reads vs 700/700 manually; none in the other 7), median depth within ~2-7%.

Still to do before relying on it:

- [ ] Build the Docker image on one Mac and one Windows laptop (the image has not been built yet).
- [ ] Run 3–5 plasmids whose sequence is already confirmed (e.g. by Sanger or a commercial whole-plasmid service), at least one Illumina and one Nanopore. Check that SeqGuide finds no false differences and catches known ones.
- [ ] Run an amplicon set with known variants (heterozygous and homozygous).
- [ ] Adjust thresholds in `decide.py` and wording in `explain.py` to lab conventions.

Ideas for later: batch mode for many plasmids at once, primer trimming for amplicons, feature transfer so the assembly is exported as an annotated GenBank file, and a Sanger trace comparison.
