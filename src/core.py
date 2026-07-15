from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable
import json, re
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

IUPAC = {
    'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]',
    'K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'
}

@dataclass
class Variant:
    sample_id: str
    kind: str
    ref_start: int
    ref: str
    alt: str
    feature: str = 'intergenic'
    consequence: str = 'not_assessed'


def load_record(path: str):
    fmt = 'genbank' if Path(path).suffix.lower() in {'.gb','.gbk','.genbank'} else 'fasta'
    return SeqIO.read(path, fmt)


def _equal_length_variants(reference, sample) -> list[Variant]:
    """Call substitutions directly when no length change is present.

    This avoids dynamic-programming alignment for the common clone-QC case
    where the reference and sample have identical lengths.
    """
    reference_seq = str(reference.seq).upper()
    sample_seq = str(sample.seq).upper()
    variants: list[Variant] = []
    index = 0
    while index < len(reference_seq):
        if reference_seq[index] == sample_seq[index]:
            index += 1
            continue
        end = index + 1
        while end < len(reference_seq) and reference_seq[end] != sample_seq[end]:
            end += 1
        kind = 'SNV' if end - index == 1 else 'MNV'
        variants.append(
            Variant(
                sample.id,
                kind,
                index + 1,
                reference_seq[index:end],
                sample_seq[index:end],
            )
        )
        index = end
    return variants


def _global_alignment(reference_seq: str, sample_seq: str):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5
    return aligner.align(reference_seq, sample_seq)[0]


def align_and_call(reference, sample) -> list[Variant]:
    reference_seq = str(reference.seq).upper()
    sample_seq = str(sample.seq).upper()

    if len(reference_seq) == len(sample_seq):
        mismatch_count = sum(old != new for old, new in zip(reference_seq, sample_seq))
        fast_path_limit = max(5, int(len(reference_seq) * 0.05))
        if mismatch_count <= fast_path_limit:
            variants = _equal_length_variants(reference, sample)
            annotate(reference, variants)
            return variants

    aln = _global_alignment(reference_seq, sample_seq)
    ref_blocks, sample_blocks = aln.aligned
    variants: list[Variant] = []
    rp = sp = 0
    for (r0, r1), (s0, s1) in zip(ref_blocks, sample_blocks):
        if r0 > rp:
            variants.append(Variant(sample.id, 'deletion', int(rp) + 1, reference_seq[rp:r0], ''))
        if s0 > sp:
            variants.append(Variant(sample.id, 'insertion', int(rp) + 1, '', sample_seq[sp:s0]))
        rseg = reference_seq[r0:r1]
        sseg = sample_seq[s0:s1]
        i = 0
        while i < len(rseg):
            if rseg[i] != sseg[i]:
                j = i + 1
                while j < len(rseg) and rseg[j] != sseg[j]:
                    j += 1
                kind = 'SNV' if j - i == 1 else 'MNV'
                variants.append(Variant(sample.id, kind, int(r0 + i) + 1, rseg[i:j], sseg[i:j]))
                i = j
            else:
                i += 1
        rp, sp = r1, s1
    if rp < len(reference_seq):
        variants.append(Variant(sample.id, 'deletion', int(rp) + 1, reference_seq[rp:], ''))
    if sp < len(sample_seq):
        variants.append(Variant(sample.id, 'insertion', int(rp) + 1, '', sample_seq[sp:]))
    annotate(reference, variants)
    return variants


def annotate(reference, variants: list[Variant]):
    for v in variants:
        pos0 = max(v.ref_start-1, 0)
        overlapping = []
        for f in reference.features:
            if int(f.location.start) <= pos0 < int(f.location.end):
                label = f.qualifiers.get('label', f.qualifiers.get('gene',[f.type]))[0]
                overlapping.append((f,label))
        if overlapping:
            f,label = overlapping[0]
            v.feature = f'{f.type}:{label}'
            if f.type == 'CDS':
                if v.kind in {'insertion','deletion'}:
                    delta = len(v.alt)-len(v.ref)
                    v.consequence = 'in_frame_indel' if delta % 3 == 0 else 'frameshift'
                elif len(v.ref) == len(v.alt) == 1:
                    start = int(f.location.start); strand = f.location.strand or 1
                    cds = f.extract(reference.seq)
                    offset = pos0-start if strand == 1 else int(f.location.end)-1-pos0
                    codon_start = (offset//3)*3
                    old_codon = str(cds[codon_start:codon_start+3]).upper()
                    if len(old_codon)==3:
                        idx = offset%3
                        new_codon = old_codon[:idx]+v.alt+old_codon[idx+1:]
                        old_aa, new_aa = str(Seq(old_codon).translate()), str(Seq(new_codon).translate())
                        if old_aa == new_aa: v.consequence='synonymous'
                        elif new_aa == '*': v.consequence='nonsense'
                        elif old_aa == '*': v.consequence='stop_lost'
                        else: v.consequence=f'missense:{old_aa}{offset//3+1}{new_aa}'


def motif_hits(seq: str, motifs: Iterable[dict]) -> dict[str,list[int]]:
    seq = seq.upper().replace('U','T')
    out = {}
    for m in motifs:
        pattern = ''.join(IUPAC.get(c,c) for c in m['pattern'].upper())
        out[m['name']] = [x.start()+1 for x in re.finditer(f'(?=({pattern}))', seq)]
    return out


def run(reference_path: str, sample_paths: list[str], motifs: list[dict]) -> dict:
    ref = load_record(reference_path)
    result = {'reference': ref.id, 'samples': [], 'variants': []}
    ref_motifs = motif_hits(str(ref.seq), motifs)
    for p in sample_paths:
        sample = load_record(p)
        vars_ = align_and_call(ref, sample)
        sample_motifs = motif_hits(str(sample.seq), motifs)
        result['samples'].append({'sample_id':sample.id,'length':len(sample.seq),'variant_count':len(vars_),
            'motif_delta':{k:len(sample_motifs[k])-len(ref_motifs[k]) for k in ref_motifs}})
        result['variants'].extend(asdict(v) for v in vars_)
    return result
