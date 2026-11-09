"""Files for genome browsers (IGV) and plasmid editors (SnapGene, Benchling), written in pure Python.

Everything here also runs in the browser version, where samtools is not available:
  * FASTA + .fai index
  * GenBank with the map's features carried over to the sample's own sequence, plus every difference
  * GFF3 annotation and BED of differences (IGV tracks)
  * BAM + .bai (reads placed on the sample's sequence), written to the SAM/BAM specification
  * an IGV session file that opens everything at once
  * a ZIP of the whole results folder
"""

from __future__ import annotations

import io
import shutil
import struct
import zipfile
import zlib
from pathlib import Path
from typing import Optional

from .refmap import Feature

# ====================================================================== FASTA / FAI

def write_fasta_with_fai(records: list, path: Path, width: int = 60) -> None:
    """records: [(name, seq)]. Writes path and path + '.fai'."""
    path = Path(path)
    fai = []
    with open(path, "w", newline="\n") as fh:
        offset = 0
        for name, seq in records:
            header = f">{name}\n"
            fh.write(header)
            offset += len(header)
            fai.append(f"{name}\t{len(seq)}\t{offset}\t{width}\t{width + 1}\n")
            for i in range(0, len(seq), width):
                line = seq[i:i + width] + "\n"
                fh.write(line)
                offset += len(line)
    Path(str(path) + ".fai").write_text("".join(fai))


# ====================================================================== coordinates & features

def ref_to_sample(variants: list):
    """Function mapping a 0-based reference position to the sample sequence (after the comparison)."""
    vs = sorted(variants, key=lambda v: v.ref_pos)

    def f(r: int) -> int:
        shift = 0
        for v in vs:
            rs, re = v.ref_pos, v.ref_pos + len(v.ref)
            if not v.ref:                       # insertion before reference position rs
                if rs <= r:
                    shift += len(v.alt)
                    continue
                break
            if re <= r:
                shift += len(v.alt) - len(v.ref)
            elif rs <= r:                       # inside a deleted / replaced block
                return v.qry_pos + min(r - rs, max(len(v.alt) - 1, 0))
            else:
                break
        return r + shift
    return f


def transfer_features(features: list, variants: list, sample_len: int) -> list:
    f = ref_to_sample(variants)
    out = []
    for feat in features:
        parts = []
        for s, e, strand in feat.parts:
            ns, ne = f(s), f(e - 1) + 1 if e > s else f(s)
            ns, ne = max(0, min(ns, sample_len)), max(0, min(ne, sample_len))
            if ne > ns:
                parts.append((ns, ne, strand))
        if parts:
            q = dict(feat.qualifiers)
            q.pop("translation", None)       # the protein may have changed; SnapGene re-translates
            out.append(Feature(feat.type, feat.label, parts, q))
    return out


def difference_features(variants: list) -> list:
    feats = []
    for v in variants:
        s = v.qry_pos
        e = v.qry_pos + max(len(v.alt), 1)
        note = f"{v.describe()} vs reference. {v.effect}".strip()
        feats.append(Feature("misc_difference", f"SeqGuide: {v.describe()}"[:60], [(s, e, 1)],
                             {"note": note[:400], "severity": v.severity or ""}))
    return feats


# ====================================================================== GenBank / GFF3 / BED

def _loc(parts: list) -> str:
    def one(s, e):
        return f"{s + 1}..{e}" if e - s > 1 else f"{s + 1}"
    strand = parts[0][2]
    if strand < 0:
        bio = [(s, e) for s, e, _ in reversed(parts)]  # back to genomic order
        inner = ",".join(one(s, e) for s, e in bio)
        return f"complement(join({inner}))" if len(bio) > 1 else f"complement({inner})"
    inner = ",".join(one(s, e) for s, e, _ in parts)
    return f"join({inner})" if len(parts) > 1 else inner


def _wrap_qual(key: str, value: str) -> list:
    value = value.replace('"', "'")
    text = f'/{key}="{value}"'
    lines, width = [], 58
    while text:
        lines.append(" " * 21 + text[:width])
        text = text[width:]
    return lines


def write_genbank(path: Path, name: str, seq: str, features: list, circular: bool = True,
                  definition: str = "") -> None:
    name = "".join(ch for ch in name if ch.isalnum() or ch in "_-")[:16] or "sequence"
    topo = "circular" if circular else "linear"
    lines = [f"LOCUS       {name:<16}{len(seq):>11} bp    DNA     {topo:<8} SYN",
             f"DEFINITION  {definition or 'Assembled by SeqGuide'}.",
             "FEATURES             Location/Qualifiers",
             f"     source          1..{len(seq)}", '                     /mol_type="other DNA"']
    for f in sorted(features, key=lambda x: x.start):
        lines.append(f"     {f.type[:15]:<16}{_loc(f.parts)}")
        lines += _wrap_qual("label", f.label)
        for k in ("note", "gene", "product"):
            if f.qualifiers.get(k) and k != "label":
                lines += _wrap_qual(k, str(f.qualifiers[k]))
    lines.append("ORIGIN")
    low = seq.lower()
    for i in range(0, len(low), 60):
        chunk = " ".join(low[j:j + 10] for j in range(i, min(i + 60, len(low)), 10))
        lines.append(f"{i + 1:>9} {chunk}")
    lines.append("//")
    Path(path).write_text("\n".join(lines) + "\n")


def write_gff3(path: Path, seqid: str, features: list) -> None:
    out = ["##gff-version 3"]
    for i, f in enumerate(sorted(features, key=lambda x: x.start), 1):
        strand = "+" if f.strand >= 0 else "-"
        label = f.label.replace(";", ",").replace("=", ":")
        for s, e, _ in f.parts:
            out.append(f"{seqid}\tSeqGuide\t{f.type}\t{s + 1}\t{e}\t.\t{strand}\t.\tID=f{i};Name={label}")
    Path(path).write_text("\n".join(out) + "\n")


def write_bed(path: Path, seqid: str, variants: list) -> None:
    colors = {"high": "208,59,59", "moderate": "236,131,90", "low": "125,124,118"}
    out = ['track name="SeqGuide differences" itemRgb="On"']
    for v in variants:
        s = v.qry_pos
        e = s + max(len(v.alt), 1)
        out.append(f"{seqid}\t{s}\t{e}\t{v.describe().replace(' ', '_')[:60]}\t0\t.\t{s}\t{e}\t"
                   f"{colors.get(v.severity, '125,124,118')}")
    Path(path).write_text("\n".join(out) + "\n")


# ====================================================================== BAM + BAI

_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
_SEQ_CODE = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}
_CIGAR = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4}


def reg2bin(beg: int, end: int) -> int:
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


class _BGZFWriter:
    MAX = 0xff00

    def __init__(self, fh):
        self.fh = fh
        self.coffset = 0          # compressed offset where the current block will start
        self.buf = bytearray()

    def tell(self) -> int:
        return (self.coffset << 16) | len(self.buf)

    def _flush(self) -> None:
        if not self.buf:
            return
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        data = c.compress(bytes(self.buf)) + c.flush()
        bsize = 18 + len(data) + 8
        header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", bsize - 1)
        tail = struct.pack("<II", zlib.crc32(bytes(self.buf)) & 0xffffffff, len(self.buf))
        self.fh.write(header + data + tail)
        self.coffset += bsize
        self.buf = bytearray()

    def write(self, data: bytes, keep_together: bool = False) -> None:
        if keep_together and len(self.buf) + len(data) > self.MAX:
            self._flush()
        while data:
            room = self.MAX - len(self.buf)
            self.buf += data[:room]
            data = data[room:]
            if len(self.buf) >= self.MAX:
                self._flush()

    def close(self) -> None:
        self._flush()
        self.fh.write(_BGZF_EOF)


def _encode_record(ref_id: int, pos: int, name: str, mapq: int, flag: int, cigar: list, seq: str,
                   qual: Optional[str]) -> tuple:
    ref_span = sum(n for n, op in cigar if op in "MDN")
    bin_ = reg2bin(pos, pos + max(ref_span, 1))
    rn = name.encode()[:250] + b"\x00"
    cig = b"".join(struct.pack("<I", n << 4 | _CIGAR[op]) for n, op in cigar)
    l_seq = len(seq)
    packed = bytearray((l_seq + 1) // 2)
    for i, ch in enumerate(seq):
        code = _SEQ_CODE.get(ch, 15)
        if i % 2 == 0:
            packed[i // 2] = code << 4
        else:
            packed[i // 2] |= code
    qb = bytes(ord(c) - 33 for c in qual) if qual else b"\xff" * l_seq
    core = struct.pack("<iiBBHHHiiii", ref_id, pos, len(rn), mapq, bin_, len(cigar), flag, l_seq, -1, -1, 0)
    body = core + rn + cig + bytes(packed) + qb
    return struct.pack("<i", len(body)) + body, bin_, ref_span


def write_bam(path: Path, refs: list, alignments: list, program: str = "SeqGuide") -> int:
    """refs: [(name, length)]; alignments: [(ref_id, pos0, cigar[(n, op)], seq, qual, reverse, mapq, name)].
    Writes a coordinate-sorted BAM and its .bai index. Returns the number of records."""
    alignments = sorted(alignments, key=lambda a: (a[0], a[1]))
    path = Path(path)
    bins: list = [dict() for _ in refs]           # per ref: bin -> [[beg, end], ...]
    linear: list = [dict() for _ in refs]         # per ref: window -> min voffset
    meta = [[None, None, 0, 0] for _ in refs]     # first voff, last end voff, mapped, unmapped
    with open(path, "wb") as fh:
        w = _BGZFWriter(fh)
        text = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{n}\tLN:{L}\n" for n, L in refs) + \
               f"@PG\tID:{program}\tPN:{program}\n"
        tb = text.encode()
        hdr = b"BAM\x01" + struct.pack("<i", len(tb)) + tb + struct.pack("<i", len(refs))
        for n, L in refs:
            nb = n.encode() + b"\x00"
            hdr += struct.pack("<i", len(nb)) + nb + struct.pack("<i", L)
        w.write(hdr)
        w._flush()   # records start in a fresh block (simplifies the index)
        for ref_id, pos, cigar, seq, qual, rev, mapq, name in alignments:
            rec, bin_, span = _encode_record(ref_id, pos, name, mapq, 16 if rev else 0, cigar, seq, qual)
            if len(w.buf) + len(rec) > w.MAX:
                w._flush()
            beg = w.tell()
            w.write(rec)
            end = w.tell()
            chunks = bins[ref_id].setdefault(bin_, [])
            if chunks and chunks[-1][1] == beg:
                chunks[-1][1] = end
            else:
                chunks.append([beg, end])
            for win in range(pos >> 14, ((pos + max(span, 1) - 1) >> 14) + 1):
                if win not in linear[ref_id]:
                    linear[ref_id][win] = beg
            m = meta[ref_id]
            if m[0] is None:
                m[0] = beg
            m[1] = end
            m[2] += 1
        w.close()
    # ---- index
    out = bytearray(b"BAI\x01" + struct.pack("<i", len(refs)))
    for r in range(len(refs)):
        b = bins[r]
        n_bin = len(b) + (1 if meta[r][0] is not None else 0)
        out += struct.pack("<i", n_bin)
        for bin_id, chunks in sorted(b.items()):
            out += struct.pack("<Ii", bin_id, len(chunks))
            for beg, end in chunks:
                out += struct.pack("<QQ", beg, end)
        if meta[r][0] is not None:
            out += struct.pack("<Ii", 37450, 2) + struct.pack("<QQ", meta[r][0], meta[r][1]) + \
                   struct.pack("<QQ", meta[r][2], meta[r][3])
        lin = linear[r]
        n_intv = (max(lin) + 1) if lin else 0
        out += struct.pack("<i", n_intv)
        last = 0
        for win in range(n_intv):
            if win in lin:
                last = lin[win]
            out += struct.pack("<Q", last)
    out += struct.pack("<Q", 0)
    Path(str(path) + ".bai").write_bytes(bytes(out))
    return len(alignments)


def alignments_from_readcheck(placed: list, ref_len: int, ref_id: int = 0, max_bases: Optional[int] = None,
                              circular: bool = True) -> list:
    """Convert readcheck placements to BAM alignments (soft-clipping the unmatched ends and any
    part running past the end of the sequence). Optionally downsample to max_bases aligned bases."""
    import random
    rng = random.Random(7)
    if max_bases:
        total = sum(p[5] for p in placed)
        if total > max_bases:
            keep = max_bases / total
            placed = [p for p in placed if rng.random() < keep]
    out = []
    for pos, rev, seq, qual, lead, alen, mapq, name in placed:
        m = alen
        if pos + m > ref_len:            # circular wrap: show the part before the origin, clip the rest
            m = ref_len - pos
        right = len(seq) - lead - m
        if m <= 0 or right < 0:
            continue
        cigar = ([(lead, "S")] if lead else []) + [(m, "M")] + ([(right, "S")] if right else [])
        out.append((ref_id, pos, cigar, seq, qual, rev, mapq, name))
    return out


# ====================================================================== IGV session / bundle / zip

def write_igv_session(path: Path, genome_file: str, locus: str, tracks: list) -> None:
    res = "\n".join(f'    <Resource path="{t}"/>' for t in tracks)
    Path(path).write_text(f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{genome_file}" locus="{locus}" version="8">
  <Resources>
{res}
  </Resources>
</Session>
''')


HOW_TO = """How to open these files
=======================

SnapGene / Benchling / ApE
    Open {gb}. It is your sample's sequence, oriented and numbered like your map, with the map's
    features carried over and every difference marked as a "misc_difference" feature (prefixed "SeqGuide:").

IGV (free genome browser, https://igv.org)
    File > Open Session... > choose igv_session.xml. IGV loads:
      {fasta} (+ .fai index)       the sequence
      {gff} / differences.bed        features and differences as tracks
{bam_line}
    All files must stay together in this folder.
{extra}"""


def make_bundle(folder: Path, title: str, entries: list, placed: Optional[dict] = None,
                bam_from: Optional[Path] = None, note: str = "") -> Path:
    """Write the genome-browser folder.

    entries: [{"name", "seq", "features", "variants", "circular"}] (one per sequence; names are used as-is
             as sequence IDs, so they must match the names inside `bam_from` if that is given)
    placed:  {entry name: readcheck placements} -> SeqGuide writes reads.bam itself (browser version)
    bam_from: a coordinate-sorted, indexed BAM made by samtools against these exact sequences (Docker version)
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    stem = "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in title)[:40] or "sample"
    fasta, gff = f"{stem}.fasta", f"{stem}.gff3"
    gff_lines, bed_lines, gbs = ["##gff-version 3"], ['track name="SeqGuide differences" itemRgb="On"'], []
    for e in entries:
        feats = list(e["features"]) + difference_features(e["variants"])
        gb = f"{stem}.gb" if len(entries) == 1 else f"{stem}_{e['name']}.gb"
        write_genbank(folder / gb, e["name"], e["seq"], feats, e["circular"], f"{title}, annotated by SeqGuide")
        gbs.append(gb)
        write_gff3(folder / "_tmp.gff3", e["name"], feats)
        gff_lines += (folder / "_tmp.gff3").read_text().splitlines()[1:]
        write_bed(folder / "_tmp.bed", e["name"], e["variants"])
        bed_lines += (folder / "_tmp.bed").read_text().splitlines()[1:]
    for t in ("_tmp.gff3", "_tmp.bed"):
        (folder / t).unlink(missing_ok=True)
    (folder / gff).write_text("\n".join(gff_lines) + "\n")
    (folder / "differences.bed").write_text("\n".join(bed_lines) + "\n")
    write_fasta_with_fai([(e["name"], e["seq"]) for e in entries], folder / fasta)
    tracks = [gff, "differences.bed"]
    bam_line = ""
    if placed:
        aln = []
        for i, e in enumerate(entries):
            pl = placed.get(e["name"]) or []
            aln += alignments_from_readcheck(pl, len(e["seq"]), ref_id=i, max_bases=600 * len(e["seq"]),
                                             circular=e["circular"])
        if aln:
            write_bam(folder / "reads.bam", [(e["name"], len(e["seq"])) for e in entries], aln)
            tracks.append("reads.bam")
            bam_line = "      reads.bam (+ .bai index)      the reads placed on your sequence (up to ~600x shown)\n"
    elif bam_from and Path(bam_from).exists():
        shutil.copy(bam_from, folder / "reads.bam")
        if Path(str(bam_from) + ".bai").exists():
            shutil.copy(str(bam_from) + ".bai", folder / "reads.bam.bai")
        tracks.append("reads.bam")
        bam_line = "      reads.bam (+ .bai index)      the reads placed on your sequence\n"
    first = entries[0]
    write_igv_session(folder / "igv_session.xml", fasta, f"{first['name']}:1-{len(first['seq'])}", tracks)
    (folder / "HOW_TO_OPEN.txt").write_text(HOW_TO.format(gb=", ".join(gbs), fasta=fasta, gff=gff, bam_line=bam_line,
                                                         extra=("\n" + note + "\n") if note else ""))
    return folder


KEEP_SUFFIXES = (".html", ".csv", ".json", ".log", ".vcf", ".tsv")
KEEP_NAMES = {"final_assembly.fasta", "assembly_matched_to_reference.fasta", "ref.fa", "ref.fa.fai",
              "aligned.bam", "aligned.bam.bai"}


def zip_results(outdir: Path, sample: str) -> bytes:
    """ZIP of the results folder: reports, tables, sequences and the genome-browser folder
    (large intermediate files such as trimmed reads are left out)."""
    outdir = Path(outdir)
    buf = io.BytesIO()
    root = f"{sample}_SeqGuide_results"
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(outdir.rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(outdir)
            top = rel.parts[0]
            if top == "genome_browser" or (len(rel.parts) == 1 and (p.suffix in KEEP_SUFFIXES or p.name in KEEP_NAMES)):
                z.write(p, f"{root}/{rel.as_posix()}")
    return buf.getvalue()
