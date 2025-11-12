#!/usr/bin/env python3
"""
build_hybrid_reference.py

Build a hybrid reference from contigs and a reference genome.

Features:
- Maps contigs onto reference using show-tiling coordinates.
- Fills gaps with reference sequence (lowercase).
- Handles overlaps with optional IUPAC blending.
- Outputs:
    - Single-line FASTA of hybrid reference
    - BED file marking source: ref_fill, contig, overlap_iupac
- Sanity check: all contigs in tiling must exist in FASTA.
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

IUPAC_CODES = {
    frozenset(['A','G']): 'R',
    frozenset(['A','T']): 'W',
    frozenset(['A','C']): 'M',
    frozenset(['G','C']): 'S',
    frozenset(['G','T']): 'K',
    frozenset(['C','T']): 'Y'
}

def iupac_merge(base1, base2):
    base1, base2 = base1.upper(), base2.upper()
    if base1 == base2:
        return base1
    if 'N' in (base1, base2):
        return 'N'
    return IUPAC_CODES.get(frozenset([base1, base2]), 'N')

def parse_show_tiling(tiling_file):
    blocks = []
    with open(tiling_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            ref_start, ref_end = int(parts[0]), int(parts[1])
            strand = parts[6]
            contig = parts[7]
            blocks.append({
                "ref_start": min(ref_start, ref_end),
                "ref_end": max(ref_start, ref_end),
                "strand": strand,
                "contig": contig
            })
    blocks.sort(key=lambda x: x["ref_start"])
    return blocks

def load_contigs(contig_fasta):
    contigs = {}
    for record in SeqIO.parse(contig_fasta, "fasta"):
        contigs[record.id] = str(record.seq)
    return contigs

def blend_sequences(seq1, seq2, use_iupac):
    if not use_iupac:
        mid = len(seq1) // 2
        return seq1[:mid] + seq2[mid:]
    else:
        return ''.join(iupac_merge(a,b) for a,b in zip(seq1, seq2))

def build_hybrid(ref_seq, contigs, blocks, sample_name, bed_out, use_iupac=False):
    hybrid_seq = []
    bed_lines = []

    prev_end = 0
    prev_contig_seq = None
    prev_ref_end = None
    prev_ref_start = None

    for b in blocks:
        if b["contig"] not in contigs:
            raise ValueError(f"Contig {b['contig']} not found in FASTA.")

        # Orient contig
        contig_seq_full = contigs[b["contig"]]
        if b["strand"] == "-":
            contig_seq_full = str(Seq(contig_seq_full).reverse_complement())

        # Fill gap before contig
        if prev_end < b["ref_start"] - 1:
            gap_seq = ref_seq[prev_end : b["ref_start"] - 1].lower()
            hybrid_seq.append(gap_seq)
            bed_lines.append(f"{sample_name}\t{prev_end}\t{b['ref_start']-1}\tref_fill\t.\t+")

        contig_seq = contig_seq_full

        # Handle overlap with previous contig
        if prev_contig_seq and b["ref_start"] - 1 <= prev_ref_end:
            overlap_start = b["ref_start"] - 1
            overlap_end = prev_ref_end
            overlap_len = overlap_end - overlap_start + 1

            seqA = prev_contig_seq[-overlap_len:]
            seqB = contig_seq[:overlap_len]
            blended = blend_sequences(seqA, seqB, use_iupac)
            hybrid_seq[-1] = hybrid_seq[-1][:-overlap_len] + blended

            # BED for overlap
            bed_lines.append(f"{sample_name}\t{overlap_start}\t{overlap_end}\toverlap_iupac\t.\t{b['strand']}")

            # Remaining post-overlap part of current contig
            post_overlap_seq = contig_seq[overlap_len:]
            if post_overlap_seq:
                hybrid_seq.append(post_overlap_seq)
                bed_lines.append(f"{sample_name}\t{overlap_end}\t{b['ref_end']-1}\tcontig\t{b['contig']}\t{b['strand']}")
        else:
            # No overlap
            hybrid_seq.append(contig_seq)
            bed_lines.append(f"{sample_name}\t{b['ref_start']-1}\t{b['ref_end']-1}\tcontig\t{b['contig']}\t{b['strand']}")

        prev_end = b["ref_end"]
        prev_contig_seq = contig_seq_full
        prev_ref_end = b["ref_end"] - 1
        prev_ref_start = b["ref_start"] - 1

    # Fill tail if reference extends beyond last contig
    if prev_end < len(ref_seq):
        tail_seq = ref_seq[prev_end:].lower()
        hybrid_seq.append(tail_seq)
        bed_lines.append(f"{sample_name}\t{prev_end}\t{len(ref_seq)}\tref_fill\t.\t+")

    # Write BED
    with open(bed_out, "w") as bedf:
        bedf.write("\n".join(bed_lines) + "\n")

    return ''.join(hybrid_seq)

def main():
    parser = argparse.ArgumentParser(description="Build hybrid reference with optional IUPAC blending.")
    parser.add_argument("ref_fa", help="Reference FASTA (single sequence)")
    parser.add_argument("contigs_fa", help="Assembly contigs FASTA")
    parser.add_argument("tiling_txt", help="Show-tiling output")
    parser.add_argument("sample_name", help="Sample name / FASTA header")
    parser.add_argument("--blend-iupac", action="store_true", help="Blend overlapping contigs using IUPAC")
    args = parser.parse_args()

    ref_record = next(SeqIO.parse(args.ref_fa, "fasta"))
    ref_seq = str(ref_record.seq)
    contigs = load_contigs(args.contigs_fa)
    blocks = parse_show_tiling(args.tiling_txt)

    bed_out = f"{args.sample_name}.hybrid.bed"
    hybrid_seq = build_hybrid(ref_seq, contigs, blocks, args.sample_name, bed_out, args.blend_iupac)

    out_fasta = f"{args.sample_name}.hybrid.fasta"
    with open(out_fasta, "w") as f:
        f.write(f">{args.sample_name}\n{hybrid_seq}\n")

    print(f"[done] Wrote {out_fasta} and {bed_out}")

if __name__ == "__main__":
    main()
