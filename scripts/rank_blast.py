#!/usr/bin/env python3
"""Reference-level BLAST evidence; scores are descriptive, not confidence values."""

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


METHOD = "ordered-nonoverlapping-hsp-chain-v1"


def read_hits(path):
    hits = {}
    query = None
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith(("#", "query acc.ver")):
            continue
        fields = line.split("\t")
        try:
            if len(fields) not in (12, 13):
                raise ValueError("expected 12 or 13 BLAST columns")
            qid, subject = fields[:2]
            if query is not None and query != qid:
                raise ValueError("multiple queries are not supported")
            query = qid
            identity = float(fields[2])
            length, mismatches, gaps, qs, qe, ss, se = map(int, fields[3:10])
            evalue, score = map(float, fields[10:12])
            if (not all(math.isfinite(v) for v in (identity, evalue, score))
                    or not 0 <= identity <= 100 or min(length, qs, qe, ss, se) < 1
                    or min(mismatches, gaps, evalue) < 0 or score <= 0
                    or max(abs(qe - qs) + 1, abs(se - ss) + 1) > length):
                raise ValueError("invalid BLAST metrics")
            # Normalize query direction; subject coordinates then encode orientation.
            if qs > qe:
                qs, qe, ss, se = qe, qs, se, ss
            direction = 1 if ss <= se else -1
            nident = int(fields[12]) if len(fields) == 13 else identity * length / 100
            if not 0 <= nident <= length:
                raise ValueError("invalid identical-base count")
            hit = dict(query_start=qs, query_end=qe, subject_start=ss,
                       subject_end=se, orientation=direction, aligned_bases=length,
                       identical_bases=nident, bit_score=score)
            hits.setdefault(subject, []).append(hit)
        except ValueError as error:
            raise ValueError(f"{path}:{number}: {error}") from error
    return hits


def best_chain(hits):
    """Dynamic programming over strictly ordered, disjoint query/subject intervals."""
    ordered = sorted(hits, key=lambda h: (h["query_end"], h["query_start"],
                                        h["orientation"], h["subject_start"],
                                        h["subject_end"], -h["bit_score"],
                                        h["aligned_bases"], -h["identical_bases"]))
    chains = []
    scores = []
    for index, hit in enumerate(ordered):
        direction = hit["orientation"]
        predecessors = [i for i, previous in enumerate(ordered[:index])
                        if previous["orientation"] == direction
                        and previous["query_end"] < hit["query_start"]
                        and direction * previous["subject_end"] < direction * hit["subject_start"]]
        best = max(predecessors, key=lambda i: scores[i], default=None)
        chains.append((chains[best] if best is not None else []) + [hit])
        scores.append((scores[best] if best is not None else 0) + hit["bit_score"])
    return chains[max(range(len(scores)), key=lambda i: scores[i])]


def rank(blast, sequence):
    references = []
    for subject, hits in sorted(read_hits(blast).items()):
        chain = best_chain(hits)
        if any(h["query_end"] > len(sequence) for h in hits):
            raise ValueError(f"BLAST coordinates exceed query length: {subject}")
        label = subject.split("_", 1)[0]
        subtype = label if re.fullmatch(r"[1-9][0-9]*[a-z]*", label) else None
        covered = sum(h["query_end"] - h["query_start"] + 1 for h in chain)
        observed = sum(base != "N" for h in chain
                       for base in sequence[h["query_start"] - 1:h["query_end"]])
        total_observed = sum(base != "N" for base in sequence)
        aligned = sum(h["aligned_bases"] for h in chain)
        identical = sum(h["identical_bases"] for h in chain)
        references.append(dict(reference=subject, subtype=subtype,
                               bit_score=sum(h["bit_score"] for h in chain),
                               identity=100 * identical / aligned,
                               aligned_bases=aligned, query_aligned_bases=covered,
                               identical_bases=identical,
                               query_coverage_pct=100 * covered / len(sequence),
                               observed_coverage_pct=100 * observed / total_observed if total_observed else None,
                               segments=chain))
    references.sort(key=lambda r: (-r["bit_score"], r["reference"]))
    groups = []
    for label in sorted({r["subtype"] for r in references if r["subtype"]}):
        members = [r for r in references if r["subtype"] == label]
        groups.append(dict(subtype=label, bit_score=members[0]["bit_score"], references=members))
    groups.sort(key=lambda g: (-g["bit_score"], g["subtype"]))
    # Unknown labels must not silently disappear and reveal a lower-scoring assignment.
    top = references[0]["bit_score"] if references else None
    winners = [g for g in groups if g["bit_score"] == top]
    unknown_top = any(r["subtype"] is None and r["bit_score"] == top for r in references)
    status = "unassigned"
    selected = None
    if len(winners) == 1 and not unknown_top:
        selected = winners[0]["subtype"]
        status = "genotype_only" if selected.isdigit() else "assigned"
    elif winners:
        status = "ambiguous"
    return dict(schema_version="1.0", method=METHOD, status=status, subtype=selected,
                selected_reference=winners[0]["references"][0]["reference"] if selected else None,
                groups=groups, unrecognized_references=[r for r in references if r["subtype"] is None],
                score_margin=groups[0]["bit_score"] - groups[1]["bit_score"] if len(groups) > 1 else None,
                query_length=len(sequence), observed_query_bases=sum(b != "N" for b in sequence))


def read_sequence(path):
    lines = Path(path).read_text().splitlines()
    if sum(line.startswith(">") for line in lines) != 1 or not lines or not lines[0].startswith(">"):
        raise ValueError("typing requires one FASTA record")
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
        raise ValueError("invalid analytical FASTA sequence")
    return sequence


def summarize(blast, fasta, database_files=()):
    sequence = read_sequence(fasta)
    result = rank(blast, sequence)
    result["input"] = dict(role="analytical_15pct_iupac", sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest())
    result["database"] = {Path(path).name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                          for path in sorted(database_files)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--database-file", action="append", default=[], type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(summarize(args.blast, args.fasta, args.database_file), indent=2) + "\n")


if __name__ == "__main__":
    main()
