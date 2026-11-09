"""Explain what each difference *does*: which feature it hits and the protein-level effect."""

from __future__ import annotations

from .plasmid_align import Variant
from .refmap import Feature, Reference
from .seqio import revcomp

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AA[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
AA_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His",
    "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Stop", "X": "?",
}

CODING_TYPES = {"CDS"}
# Non-coding feature types whose sequence usually matters for function
FUNCTIONAL_TYPES = {"promoter", "RBS", "rep_origin", "terminator", "regulatory", "protein_bind", "enhancer",
                    "polyA_signal", "LTR", "oriT", "sig_peptide", "5'UTR", "3'UTR", "intron", "exon", "mRNA",
                    "gene", "misc_RNA", "ncRNA", "primer_bind"}

SEVERITY_ORDER = {"high": 0, "moderate": 1, "low": 2, "": 3}


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def _cds_positions(f: Feature) -> list:
    pos = []
    for s, e, strand in f.parts:
        pos.extend(range(s, e) if strand >= 0 else range(e - 1, s - 1, -1))
    return pos


def _cds_seq(ref: str, f: Feature, positions: list, overrides: dict | None = None) -> str:
    overrides = overrides or {}
    strand = f.strand
    out = []
    for p in positions:
        b = overrides.get(p, ref[p])
        out.append(b if strand >= 0 else revcomp(b))
    return "".join(out)


def coding_effect(ref: Reference, f: Feature, v: Variant) -> tuple:
    """Return (effect sentence, severity) for a variant inside a CDS."""
    positions = _cds_positions(f)
    index = {p: i for i, p in enumerate(positions)}
    codon_start = int(f.qualifiers.get("codon_start", "1") or 1) - 1
    name = f.label

    if v.structural or v.kind in ("COMPLEX", "INVERSION") or len(v.ref) > 30 or len(v.alt) > 30:
        return f"Large change inside {name}: the protein is very likely disrupted.", "high"

    if v.kind in ("SNP", "MNP"):
        overrides = {v.ref_pos + i: v.alt[i] for i in range(len(v.ref))}
        idxs = [index[p] for p in overrides if p in index]
        if not idxs:
            return f"In {name} (outside its coding positions).", "low"
        ref_cds = _cds_seq(ref.seq, f, positions)[codon_start:]
        alt_cds = _cds_seq(ref.seq, f, positions, overrides)[codon_start:]
        ref_p, alt_p = translate(ref_cds), translate(alt_cds)
        changes = [(i, a, b) for i, (a, b) in enumerate(zip(ref_p, alt_p)) if a != b]
        first_codon = (min(idxs) - codon_start) // 3 + 1
        if not changes:
            return f"Silent change in {name} (codon {first_codon}): the protein is unchanged.", "low"
        parts, sev = [], "moderate"
        for i, a, b in changes:
            lbl = f"{a}{i + 1}{b}"
            if i == 0 and a == "M":
                parts.append(f"{lbl}: start codon lost, protein probably not made")
                sev = "high"
            elif b == "*":
                parts.append(f"{lbl}: premature stop codon, protein is truncated at amino acid {i + 1} of {len(ref_p)}")
                sev = "high"
            elif a == "*":
                parts.append(f"{lbl}: stop codon lost, protein is extended")
                sev = "high"
            else:
                parts.append(f"{lbl}: {AA_NAMES.get(a, a)} → {AA_NAMES.get(b, b)} (amino-acid change)")
        return f"In {name}: " + "; ".join(parts) + ".", sev

    # indels
    anchor = v.ref_pos if v.kind == "DEL" else v.ref_pos - 1
    i = index.get(anchor, index.get(v.ref_pos))
    aa = (i - codon_start) // 3 + 1 if i is not None else None
    delta = len(v.alt) - len(v.ref)
    where = f" at amino acid {aa}" if aa else ""
    if delta % 3:
        return (f"Frameshift in {name}{where}: {abs(delta)} bp {'inserted' if delta > 0 else 'deleted'}, "
                "so every amino acid after it is wrong. The protein is very likely non-functional."), "high"
    n = abs(delta) // 3
    return (f"In-frame {'insertion' if delta > 0 else 'deletion'} of {n} amino acid{'s' if n > 1 else ''} in {name}{where}. "
            "The protein keeps its reading frame but may be altered."), "moderate"


def annotate(ref: Reference, variants: list, platform: str = "illumina") -> list:
    for v in variants:
        end = v.ref_pos + len(v.ref)
        hits = [f for f in ref.features if f.overlaps(v.ref_pos, end)]
        v.features = [f"{f.label} ({f.type})" for f in hits]
        cds = [f for f in hits if f.type in CODING_TYPES]
        func = [f for f in hits if f.type in FUNCTIONAL_TYPES and f.type != "gene"]
        if cds:
            effects = [coding_effect(ref, f, v) for f in cds]
            effects.sort(key=lambda e: SEVERITY_ORDER[e[1]])
            v.effect = " ".join(e[0] for e in effects)
            v.severity = effects[0][1]
        elif func:
            f = func[0]
            v.effect = (f"Inside {f.label} ({f.type.replace('_', ' ')}). Non-coding, but this element can control "
                        "expression or replication, so the change may matter.")
            v.severity = "moderate"
        elif hits:
            v.effect = f"Inside {hits[0].label} ({hits[0].type}); effect unclear."
            v.severity = "low"
        elif ref.annotated:
            v.effect = "Outside all annotated features; usually harmless."
            v.severity = "low"
        else:
            v.effect = "Reference has no annotation, so the effect cannot be predicted. Use a GenBank/SnapGene map for this."
            v.severity = "moderate" if v.structural or v.size > 10 else "low"
        if v.structural and v.kind == "INVERSION":
            v.effect = "A block of sequence is flipped in orientation. " + v.effect
            v.severity = "high" if hits else "moderate"

        if platform in ("nanopore", "long") and v.kind in ("INS", "DEL") and v.homopolymer >= 4 and v.size <= 2:
            v.flags.append(f"In a run of {v.homopolymer} identical letters: Nanopore often miscounts these. "
                           "Possibly a sequencing artifact — confirm (e.g. Sanger) before acting on it.")
    variants.sort(key=lambda x: x.ref_pos)
    return variants


def verdict(variants: list, ref_covered: float, identity: float, notes: list) -> dict:
    """One-line overall verdict for the top of the report."""
    high = [v for v in variants if v.severity == "high"]
    mod = [v for v in variants if v.severity == "moderate"]
    if not variants and ref_covered > 99.5:
        return {"level": "ok", "text": "Perfect match: your plasmid is identical to the reference."}
    if high:
        return {"level": "block", "text": f"{len(high)} difference(s) likely break a feature (e.g. frameshift or stop codon). "
                                          "Check these before using the plasmid."}
    if mod:
        return {"level": "warn", "text": f"{len(variants)} difference(s) found; {len(mod)} may affect function. Review them below."}
    return {"level": "ok", "text": f"{len(variants)} minor difference(s), none expected to affect function."}
