"""Step 2 of the wizard: recommend an analysis path, and say *why*.

The rules are deliberately simple and readable so lab members can adjust them.
Every rule produces a sentence a non-expert can understand.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Optional

from .inspect import SampleGroup

# Thresholds (edit here to change lab policy)
MIN_DEPTH_PLASMID_SHORT = 30
MIN_DEPTH_PLASMID_LONG = 20
MIN_DEPTH_AMPLICON = 100
SUBSAMPLE_TARGET_DEPTH = 150
MIN_Q_ILLUMINA = 25
MIN_Q_NANOPORE = 10
SIZE_MISMATCH_FRACTION = 0.10

SAMPLE_TYPES = {
    "plasmid": "A plasmid (whole plasmid was sequenced)",
    "amplicon": "PCR amplicons or targeted genes",
    "other": "Something else (genome, transcriptome, metagenome…)",
}
REFERENCE_TYPES = {
    "annotated": "Yes — an annotated map (GenBank .gb or SnapGene .dna)",
    "sequence": "Yes — sequence only (FASTA)",
    "none": "No reference",
}
GOALS = {
    "verify": "Check that my construct / gene is correct",
    "find_mutations": "Find mutations compared to a known sequence",
    "unknown": "Find out what the sequence is",
}


@dataclass
class Message:
    level: str  # "ok", "info", "warn", "block"
    text: str
    learn_more: Optional[str] = None  # key in explain.GLOSSARY


@dataclass
class Option:
    workflow: str
    title: str
    summary: str
    pros: list
    cons: list


@dataclass
class Recommendation:
    workflow: Optional[str]          # "plasmid_denovo", "amplicon_align", "plasmid_align" or None
    title: str
    reasons: list                    # list[Message]
    messages: list                   # warnings and notes, list[Message]
    alternatives: list               # list[Option]
    params: dict = field(default_factory=dict)
    blocked: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


def _effective_bases(sample: SampleGroup, kind: str) -> int:
    total = 0
    for f in sample.files:
        is_long = f.platform in ("nanopore", "pacbio", "long")
        if (kind == "long") != is_long:
            continue
        if f.sampled and f.estimated_total_reads:
            total += int(f.mean_len * f.estimated_total_reads)
        else:
            total += f.total_bases
    return total


def _platform_of(sample: SampleGroup) -> str:
    if sample.layout in ("paired", "single"):
        return "illumina"
    if sample.layout == "long":
        return "nanopore"
    if sample.layout == "hybrid":
        return "hybrid"
    return "unknown"


def recommend(
    sample: SampleGroup,
    sample_type: str,
    reference: str,
    goal: str,
    expected_size: Optional[int] = None,
) -> Recommendation:
    reasons: list[Message] = []
    msgs: list[Message] = []
    alts: list[Option] = []
    params: dict = {"platform": _platform_of(sample), "layout": sample.layout}
    blocked = False

    # ---- universal data checks
    if sample.layout == "unknown":
        msgs.append(Message("block", "SeqGuide could not work out how the files of this sample belong together. "
                                     "Check that paired files end in _R1/_R2 (or _1/_2).", "paired"))
        blocked = True
    if sample.pair_check and not sample.pair_check["ok"]:
        msgs.append(Message("block", sample.pair_check["message"], "paired"))
        blocked = True
    elif sample.pair_check:
        msgs.append(Message("ok", sample.pair_check["message"], "paired"))
    if sample.layout == "single":
        msgs.append(Message("info", "Single-end reads: analysis works, but paired-end gives more reliable results.", "single"))
    for f in sample.files:
        for p in f.problems:
            msgs.append(Message("warn", f"{f.name}: {p}"))
        if f.mean_q is not None:
            is_long = f.platform in ("nanopore", "pacbio", "long")
            limit = MIN_Q_NANOPORE if is_long else MIN_Q_ILLUMINA
            if f.mean_q < limit:
                msgs.append(Message("warn", f"{f.name}: average quality is Q{f.mean_q:.0f}, below the usual Q{limit}. "
                                            "Trimming will help, but expect more uncertain positions.", "mean_q"))
        if f.n_reads and f.n_reads < 1000:
            msgs.append(Message("warn", f"{f.name}: only {f.n_reads} reads. This is very little data."))

    # ---- sample-type specific logic
    if sample_type == "other":
        return Recommendation(
            None, "Not supported yet",
            [Message("info", "This first version of SeqGuide handles plasmids and amplicons/targeted genes. "
                             "Whole genomes, RNA-seq and metagenomes need more computing power and different tools.")],
            msgs, [], params, blocked=True,
        )

    if sample_type == "plasmid":
        long_size = next((f.plasmid_size_estimate for f in sample.files if f.plasmid_size_estimate), None)
        size = expected_size or long_size
        params["expected_size"] = size
        params["compare"] = reference != "none"
        params["subsample_depth"] = SUBSAMPLE_TARGET_DEPTH
        if sample.layout in ("paired", "single"):
            params["assembler"] = "unicycler"
        elif sample.layout == "long":
            params["assembler"] = "flye+medaka"
        elif sample.layout == "hybrid":
            params["assembler"] = "unicycler-hybrid"

        reasons.append(Message("ok", "Plasmids are small, so a complete de novo assembly takes only minutes on a laptop.", "de_novo"))
        reasons.append(Message("ok", "Assembling from scratch reveals insertions, deletions, rearrangements or an unexpected "
                                     "backbone — changes that aligning reads to the expected map can miss.", "alignment"))
        if reference != "none":
            reasons.append(Message("ok", "Because you have a reference, the finished assembly will then be compared to it and "
                                         "every difference will be listed.", "comparison"))
            if reference == "annotated":
                reasons.append(Message("ok", "Your reference is annotated, so each difference will be placed in its feature "
                                             "(e.g. 'missense change in AmpR').", "genbank"))
            else:
                msgs.append(Message("info", "With a FASTA reference you will see *where* differences are but not which gene "
                                            "they hit. An annotated GenBank/SnapGene map gives more useful results.", "genbank"))
        else:
            reasons.append(Message("info", "Without a reference you will get the assembled sequence and its quality. "
                                           "Add the expected map later to compare.", "comparison"))
        if params.get("assembler") == "flye+medaka":
            msgs.append(Message("info", "Nanopore-only assemblies can contain 1-letter errors in runs like AAAAAA. "
                                        "These will be flagged in the report.", "homopolymer"))

        if long_size and expected_size:
            diff = abs(long_size - expected_size) / expected_size
            if diff > SIZE_MISMATCH_FRACTION:
                msgs.append(Message("warn", f"Reads suggest a plasmid of ~{long_size:,} bp, but you expected {expected_size:,} bp "
                                            f"({diff:.0%} different). The plasmid may carry an insertion/deletion, or this may be "
                                            "a different plasmid.", "plasmid_size"))
            else:
                msgs.append(Message("ok", f"Read-length peak (~{long_size:,} bp) matches your expected size.", "plasmid_size"))
        if any(f.dimer_peak for f in sample.files):
            msgs.append(Message("info", "A second peak at twice the plasmid size was seen (dimers). This is common and usually harmless.", "dimer"))

        if size:
            for kind, lim in (("short", MIN_DEPTH_PLASMID_SHORT), ("long", MIN_DEPTH_PLASMID_LONG)):
                bases = _effective_bases(sample, kind)
                if not bases:
                    continue
                d = bases / size
                params[f"depth_{kind}"] = round(d, 1)
                if d < lim:
                    msgs.append(Message("warn", f"Estimated {kind}-read depth is only {d:.0f}× (recommended ≥ {lim}×). "
                                                "The assembly may be incomplete.", "depth"))
                elif d > 3 * SUBSAMPLE_TARGET_DEPTH:
                    msgs.append(Message("info", f"Estimated {kind}-read depth is {d:,.0f}×, much more than needed. "
                                                f"Reads will be subsampled to ~{SUBSAMPLE_TARGET_DEPTH}× to keep the assembler happy and fast.",
                                        "subsample"))
                else:
                    msgs.append(Message("ok", f"Estimated {kind}-read depth is {d:,.0f}× — enough for assembly.", "depth"))
            if sample.layout in ("paired", "single"):
                msgs.append(Message("info", "For Illumina plasmid preps some reads come from host (E. coli) DNA, so real plasmid "
                                            "depth will be confirmed after assembly.", "depth"))
        else:
            msgs.append(Message("info", "Enter the expected plasmid size to get a depth estimate before running.", "depth"))

        title = "De novo plasmid assembly" + (" + comparison to your reference" if reference != "none" else "")
        if reference != "none":
            alts.append(Option(
                "plasmid_align", "Align reads to the reference only",
                "Map reads onto your expected plasmid and call differences.",
                ["Slightly faster", "Good for spotting single-letter changes"],
                ["Blind to inserted sequence that is not in your map", "Rearrangements and wrong backbones can go unnoticed"],
            ))
        return Recommendation("plasmid_denovo", title, reasons, msgs, alts, params, blocked)

    # ---- amplicons / targeted genes
    if reference == "none":
        msgs.append(Message("block", "Amplicon analysis needs the expected sequence(s) of your targets. Download the gene or "
                                     "region from NCBI (FASTA or GenBank) or export it from SnapGene, then add it here.", "alignment"))
        alts.append(Option(
            "plasmid_denovo", "Experimental: assemble amplicons de novo",
            "Try to build the amplicon sequence from reads alone.",
            ["Works without a reference"],
            ["Amplicon data has very uneven depth, which assemblers handle poorly", "Results need manual checking"],
        ))
        return Recommendation(None, "Reference needed", [], msgs, alts, params, blocked=True)

    params["preset"] = "map-ont" if sample.layout == "long" else "sr"
    params["min_depth"] = MIN_DEPTH_AMPLICON
    reasons.append(Message("ok", "Amplicons are short and their expected sequence is known, so aligning reads to it is the "
                                 "fastest and most sensitive way to find point mutations and small indels.", "alignment"))
    reasons.append(Message("ok", "Each target gets its own depth check, so you can see which amplicons worked.", "depth"))
    if reference == "annotated":
        reasons.append(Message("ok", "Your reference is annotated, so protein-level effects will be reported.", "genbank"))
    if sample.layout == "long":
        msgs.append(Message("warn", "Variant calls from Nanopore-only amplicon data should be confirmed (e.g. Sanger), "
                                    "particularly indels in runs like AAAAAA.", "homopolymer"))
    if expected_size:
        bases = _effective_bases(sample, "long" if sample.layout == "long" else "short")
        d = bases / expected_size
        params["depth_est"] = round(d, 1)
        level = "ok" if d >= MIN_DEPTH_AMPLICON else "warn"
        msgs.append(Message(level, f"Average depth across targets is roughly {d:,.0f}× "
                                   f"(recommended ≥ {MIN_DEPTH_AMPLICON}× for confident variant calls).", "depth"))
    alts.append(Option(
        "plasmid_denovo", "Assemble amplicons de novo instead",
        "Build the sequence without the reference, then compare.",
        ["Can reveal large unexpected insertions"],
        ["Uneven amplicon depth makes assembly unreliable", "Slower"],
    ))
    return Recommendation("amplicon_align", "Align reads to your target sequences", reasons, msgs, alts, params, blocked)
