"""Plain-language explanations shown next to every number and every choice.

Keep each entry short: what it is, then why it matters for *this* user.
Edit freely: this file is the "teaching voice" of the lab pipeline.
"""

GLOSSARY = {
    "fastq": (
        "FASTQ file",
        "Raw reads straight from the sequencer. Each read has its DNA letters plus a quality score for every letter.",
        "This is the normal starting point for analysis.",
    ),
    "fasta": (
        "FASTA file",
        "DNA sequences without quality scores.",
        "Usually a reference (e.g. your expected plasmid) or an assembled result, not raw reads.",
    ),
    "genbank": (
        "GenBank / SnapGene file",
        "A sequence with its features annotated (genes, promoters, resistance markers, primers).",
        "Best kind of reference: SeqGuide can tell you *which feature* a mutation falls in and what it does to the protein.",
    ),
    "paired": (
        "Paired-end reads (R1 / R2)",
        "Illumina read each DNA fragment from both ends, giving two files that belong together.",
        "Pairs carry extra information about distance between reads, which makes assembly and alignment more reliable. Both files must always be used together.",
    ),
    "single": (
        "Single-end reads",
        "Each DNA fragment was read from one end only.",
        "Works, but assemblies of repetitive regions are less reliable than with paired reads.",
    ),
    "long": (
        "Long reads (Nanopore / PacBio)",
        "Reads of thousands of bases; for whole-plasmid sequencing one read often covers the entire plasmid.",
        "Great for getting the full structure, but single reads have more errors, especially in runs of the same letter (AAAAAA).",
    ),
    "hybrid": (
        "Hybrid data",
        "Both short (Illumina) and long (Nanopore) reads for the same sample.",
        "The long reads give the structure, the short reads fix the small errors. Usually the most accurate option.",
    ),
    "read_length": (
        "Read length",
        "How many DNA letters each read contains.",
        "Illumina reads are typically 75–300 bp. Long reads peak at your plasmid size if the whole plasmid was sequenced.",
    ),
    "n_reads": (
        "Number of reads",
        "How many DNA fragments were sequenced.",
        "More reads means more depth; see 'Depth'.",
    ),
    "depth": (
        "Depth (coverage)",
        "How many times, on average, each position of your target was read. Calculated as total sequenced bases ÷ target size.",
        "Each position needs to be read many times to tell real mutations from sequencing errors. "
        "Rough guide: plasmid assembly ≥ 30× (Illumina) or ≥ 20× (Nanopore); variant calling in amplicons ≥ 100×.",
    ),
    "mean_q": (
        "Mean quality (Phred score)",
        "The sequencer's confidence in each letter. Q20 = 1 error in 100 letters, Q30 = 1 in 1,000.",
        "Illumina data is usually Q30+. Nanopore is typically Q10–Q20 per read, which is fine because many reads are combined.",
    ),
    "q30": (
        "% bases ≥ Q30",
        "Share of letters with at most a 1-in-1,000 error chance.",
        "For Illumina, above ~75% is good. Low values point to a problematic run.",
    ),
    "gc": (
        "GC content",
        "Percentage of G and C letters.",
        "Mostly informative. Very extreme values can cause uneven coverage.",
    ),
    "n50": (
        "N50",
        "Half of all sequenced bases are in reads at least this long.",
        "A standard way to summarise read lengths for long-read data.",
    ),
    "plasmid_size": (
        "Estimated plasmid size",
        "In whole-plasmid Nanopore data, most reads are full-length copies, so the read-length peak equals the plasmid size.",
        "If this is far from your expected size, the plasmid may have an insertion, a deletion, or be a different plasmid.",
    ),
    "dimer": (
        "Dimer peak",
        "A second read-length peak at twice the plasmid size.",
        "Usually harmless plasmid dimers or concatemers. Worth knowing if the assembly looks twice as long as expected.",
    ),
    "de_novo": (
        "De novo assembly",
        "Building the sequence from the reads alone, like a puzzle without the box picture.",
        "Finds insertions, deletions, rearrangements and unexpected sequence that a reference would hide. Recommended for plasmids.",
    ),
    "alignment": (
        "Alignment to a reference",
        "Placing each read onto a known sequence and looking for differences.",
        "Fast and sensitive for small changes (point mutations, small indels), but it cannot see sequence that is missing from the reference.",
    ),
    "comparison": (
        "Comparing an assembly to a reference",
        "Lining up your assembled plasmid with the plasmid map you expected.",
        "Shows every difference and which feature (gene, promoter, marker) it falls in.",
    ),
    "circular": (
        "Circular contig",
        "The assembler found that the sequence's end joins its start, as expected for a plasmid.",
        "A single circular contig of the expected size is the ideal result.",
    ),
    "frameshift": (
        "Frameshift",
        "An insertion or deletion whose length is not a multiple of 3 inside a protein-coding gene.",
        "Changes every amino acid after it and usually destroys the protein.",
    ),
    "homopolymer": (
        "Homopolymer",
        "A run of the same letter, e.g. AAAAAA.",
        "Nanopore often miscounts these, so a 1-letter insertion/deletion there may be a sequencing artifact rather than a real mutation.",
    ),
    "subsample": (
        "Subsampling",
        "Using only part of the reads.",
        "Plasmids are tiny, so runs often give 1000×+ depth. Too much depth can confuse assemblers; ~100–150× is plenty.",
    ),
    "trimming": (
        "Quality trimming",
        "Removing adapter sequences and low-quality ends from reads.",
        "Prevents sequencing artefacts from appearing as fake mutations.",
    ),
    "variant": (
        "Variant",
        "A position where your sample differs from the reference.",
        "SNP = one letter changed; insertion/deletion (indel) = letters added or removed.",
    ),
    "allele_fraction": (
        "Allele fraction",
        "Share of reads at that position that show the variant.",
        "~100% = all copies carry it. 30–70% may mean a mixed population (two clones) — consider re-picking a single colony.",
    ),
}


def explain(key: str) -> dict:
    title, what, why = GLOSSARY[key]
    return {"key": key, "title": title, "what": what, "why": why}
