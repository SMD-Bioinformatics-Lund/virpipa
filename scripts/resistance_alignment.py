"""Conservative nucleotide projection for H77-numbered resistance sites."""
import hashlib
import json
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices

BASES = dict(zip("ACGTRYSWKMBDHVN", ("A", "C", "G", "T", "AG", "CT", "CG", "AT", "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT")))


def _aligner():
    matrix = substitution_matrices.Array(alphabet="".join(BASES), dims=2)
    for a in BASES:
        for b in BASES:
            matrix[a, b] = 0 if "N" in (a, b) else (2 if set(BASES[a]) & set(BASES[b]) else -3)
    aligner = PairwiseAligner(substitution_matrix=matrix, open_gap_score=-10, extend_gap_score=-1)
    aligner.end_gap_score = 0
    return aligner


def _segment(reference, query, left_terminal, right_terminal):
    aligner = _aligner()
    if not left_terminal:
        aligner.open_left_gap_score = -10
        aligner.extend_left_gap_score = -1
    if not right_terminal:
        aligner.open_right_gap_score = -10
        aligner.extend_right_gap_score = -1
    alignments = aligner.align(reference, query)
    shared = None
    unsafe = len(reference)
    count = 0
    exhausted = False
    first = None
    for alignment in alignments:
        count += 1
        if count > 256:
            exhausted = True
            break
        if first is None:
            first = alignment
        mapping = {}
        coords = alignment.coordinates
        for i in range(coords.shape[1] - 1):
            a, b = map(int, coords[0, i:i+2])
            c, d = map(int, coords[1, i:i+2])
            if b > a and d > c:
                mapping.update((a + j, c + j) for j in range(b - a))
            elif ((b-a) or (d-c)) % 3:
                terminal = (i == 0 and left_terminal) or (i == coords.shape[1]-2 and right_terminal)
                if not terminal:
                    unsafe = min(unsafe, a)
        if shared is None:
            shared = mapping
        else:
            shared = {p: q for p, q in shared.items() if mapping.get(p) == q}
    if exhausted:
        shared = {}
        unsafe = 0
    shared = shared or {}
    return {
        "mapping": shared, "unsafe_from": unsafe,
        "reference_aligned": str(first[0]), "sample_aligned": str(first[1]),
        "optimal_alignments_checked": min(count, 256), "alignment_limit_exceeded": exhausted,
    }


def align(reference, query):
    """Anchor on unique exact 9-mers; assess uncertainty within each interval.

    Anchors must agree with a whole-gene optimal alignment. Partitioning prevents
    a repetitive/unknown gene tail from hiding trustworthy upstream codons.
    """
    first = _aligner().align(reference, query)[0]
    indices = first.indices
    anchors = []
    last_r = last_q = 0
    for col in range(indices.shape[1]-8):
        r, q = map(int, indices[:, col])
        if r < last_r or q < last_q or r < 0 or q < 0:
            continue
        if list(indices[0, col:col+9]) != list(range(r, r+9)) or list(indices[1, col:col+9]) != list(range(q, q+9)):
            continue
        word = reference[r:r+9]
        if set(word) <= set("ACGT") and word == query[q:q+9] and reference.find(word) == reference.rfind(word) and query.find(word) == query.rfind(word):
            anchors.append((r, q))
            last_r, last_q = r+9, q+9
    mapping = {}
    a_parts, b_parts = [], []
    previous_r = previous_q = 0
    checked, exceeded = 1, False
    uncertain_from = len(reference)
    for r, q in anchors + [(len(reference), len(query))]:
        a, b = reference[previous_r:r], query[previous_q:q]
        if a and b:
            segment = _segment(a, b, previous_r == previous_q == 0, r == len(reference) and q == len(query))
            mapping.update((p+previous_r, s+previous_q) for p, s in segment["mapping"].items())
            a_parts.append(segment["reference_aligned"])
            b_parts.append(segment["sample_aligned"])
            checked = max(checked, segment["optimal_alignments_checked"])
            exceeded |= segment["alignment_limit_exceeded"]
            if segment["unsafe_from"] < len(a):
                uncertain_from = min(uncertain_from, previous_r+segment["unsafe_from"])
        else:
            a_parts.append(a or "-"*len(b))
            b_parts.append(b or "-"*len(a))
        if r < len(reference):
            mapping.update((r+i, q+i) for i in range(9))
            a_parts.append(reference[r:r+9]); b_parts.append(query[q:q+9])
        previous_r, previous_q = r+9, q+9
    ref_aligned, query_aligned = "".join(a_parts), "".join(b_parts)
    unsafe = uncertain_from
    ref_pos, column = 0, 0
    while column < len(ref_aligned):
        if "-" not in (ref_aligned[column], query_aligned[column]):
            ref_pos += 1; column += 1
            continue
        start = column
        deletion = query_aligned[column] == "-"
        while column < len(ref_aligned) and (query_aligned[column] == "-" if deletion else ref_aligned[column] == "-"):
            column += 1
        size = column-start
        if start > 0 and column < len(ref_aligned) and size % 3:
            unsafe = min(unsafe, ref_pos)
        if deletion:
            ref_pos += size
    deletions = set()
    for p in range(0, len(reference)-2, 3):
        if all(p+j not in mapping for j in range(3)) and p-1 in mapping and p+3 in mapping and mapping[p+3] == mapping[p-1]+1:
            deletions.add(p)
    return {"mapping": mapping, "deletions": deletions, "unsafe_from": unsafe,
            "reference_aligned": ref_aligned, "sample_aligned": query_aligned,
            "optimal_alignments_checked": checked, "alignment_limit_exceeded": exceeded,
            "method": "unique-9nt-anchored affine nucleotide alignment",
            "parameters": {"compatible": 2, "mismatch": -3, "unknown": 0, "gap_open": -10, "gap_extend": -1}}


def reference_genes(panel_path, subtype, h77, read_fasta):
    """Use exact assigned subtypes only; unknown subtypes fall back to H77."""
    if not panel_path:
        return {gene: {"nuc": item["nuc"], "accession": "AF011751.1",
                       "h77_map": dict(enumerate(range(len(item["nuc"]))))}
                for gene, item in h77.items()}
    path = Path(panel_path)
    panel = json.loads(path.read_text())
    if panel["numbering_reference"] != "AF011751.1":
        raise ValueError("Resistance panel must link to AF011751.1")
    result = reference_genes(None, subtype, h77, read_fasta)
    candidates = [r for r in panel["references"] if r.get("subtype") == subtype and any(c.isalpha() for c in subtype)]
    if len(candidates) > 1:
        raise ValueError(f"Multiple resistance references for {subtype}")
    if not candidates:
        return result
    selected = candidates[0]
    records = read_fasta(path.parent / selected["genome_fasta"])
    genome = records[selected["accession"]]
    if len(genome) != selected["length"] or hashlib.sha256(genome.encode()).hexdigest() != selected["sequence_sha256"]:
        raise ValueError("Resistance reference sequence does not match panel metadata")
    for gene, feature in selected["regions"].items():
        coding = genome[feature["start"]-1:feature["end"]]
        alignment = read_fasta(path.parent / panel["region_alignments"][gene])
        numbered, bridge = alignment["AF011751.1"], alignment[selected["accession"]]
        if (len(numbered) != len(bridge) or numbered.replace("-", "") != h77[gene]["nuc"]
                or bridge.replace("-", "") != coding or len(coding) % 3 or feature["strand"] != "+"):
            raise ValueError(f"Invalid resistance linking alignment for {gene}")
        h, s = -1, -1
        mapping = {}
        for a, b in zip(numbered, bridge):
            h += a != "-"
            s += b != "-"
            if a != "-" and b != "-":
                mapping[h] = s
        result[gene] = {"nuc": coding, "accession": selected["accession"], "h77_map": mapping}
    return result
