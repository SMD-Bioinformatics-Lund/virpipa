#!/usr/bin/env python3
"""Call HCV resistance from a 15% IUPAC consensus and optional read evidence."""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from resistance_alignment import align, reference_genes
from consensus_metadata import load as load_consensus

try:
    import pysam
except ImportError:
    pysam = None
from update_geno2pheno_rules import (
    load_rules_json,
    load_rules_rows_from_csv,
    normalize_rows,
)

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
TABLE = {
    a + b + c: x
    for a, row in zip(
        "TCAG",
        (
            "FFLLSSSSYY**CC*W",
            "LLLLPPPPHHQQRRRR",
            "IIIMTTTTNNKKSSRR",
            "VVVVAAAADDEEGGGG",
        ),
    )
    for b, part in zip("TCAG", (row[0:4], row[4:8], row[8:12], row[12:16]))
    for c, x in zip("TCAG", part)
}
H77 = {"NS3": (3420, 5312), "NS5A": (6258, 7601), "NS5B": (7602, 9374)}
RANK = {"resistant": 0, "reduced susceptibility": 1}


def fasta(path):
    out = {}
    name = None
    parts = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if name is not None:
        out[name] = "".join(parts).upper()
    if not out:
        raise ValueError(f"No sequence found in {path}")
    return out


def gff(path):
    out = {}
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        fields = raw.rstrip().split("\t")
        if raw.startswith("#") or len(fields) != 9 or fields[2] != "gene":
            continue
        attrs = {}
        for item in fields[8].split(";"):
            separator = "=" if "=" in item else (":" if ":" in item else None)
            if separator:
                key, value = item.split(separator, 1)
                attrs[key] = value
        if attrs.get("gene"):
            if attrs["gene"] in out:
                raise ValueError(f"Multiple annotations for gene {attrs['gene']}")
            out[attrs["gene"]] = {
                "chrom": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "phase": fields[7],
                "attributes": attrs,
            }
    return out


def rc(seq):
    return seq.translate(str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN"))[::-1]


def codons(codon):
    result = [""]
    for char in codon.upper():
        result = [prefix + base for prefix in result for base in IUPAC.get(char, "N")]
    return result


def aas(codon):
    return (
        {TABLE.get(item, "X") for item in codons(codon)} if len(codon) == 3 else set()
    )


def protein(seq):
    result = []
    for offset in range(0, len(seq) - 2, 3):
        choices = sorted(aas(seq[offset : offset + 3]) - {"X"})
        result.append(choices[0] if len(choices) == 1 else "X")
    return "".join(result)


def parts(definition):
    result = []
    for raw in re.split(r"\s+and\s+", (definition or "").strip(), flags=re.I):
        alternatives = [re.fullmatch(r"(\d+)\s*(del|any|[A-Za-z*])", item.strip()) for item in re.split(r"\s+or\s+", raw, flags=re.I)]
        if not all(alternatives) or len({m.group(1) for m in alternatives}) != 1:
            raise ValueError(f"Unsupported resistance rule: {definition}")
        options = sorted({m.group(2) for m in alternatives})
        result.append(
            {"position": int(alternatives[0].group(1)), "aa": "/".join(options),
             "aa_options": options, "raw": raw.strip()}
        )
    return result


def matching_amino_acids(part, site):
    observed = set(site.get("amino_acids", []))
    if "any" in part["aa_options"]:
        # The catalogue's 168any denotes a substitution, not wild-type D168.
        return observed - {site["h77_aa"], "del", "*", "X"}
    return observed & set(part["aa_options"])


def genotype(subtype):
    match = re.match(r"\d+", subtype.strip())
    return match.group() if match else ""


def subtype_match(subtype, selector):
    return any(
        item == subtype.lower() or (item.isdigit() and item == genotype(subtype))
        for item in (x.strip().lower() for x in (selector or "").split(","))
    )


def licensed(subtype, selector):
    return genotype(subtype) in {x.strip() for x in (selector or "").split(",")}


def rules(path):
    if Path(path).suffix.lower() == ".json":
        return load_rules_json(path)["rules"]
    _, rows = load_rules_rows_from_csv(path)
    data = normalize_rows(rows)
    return [dict(zip(data["columns"], row)) for row in data["rules"]]


def h77_genes(path):
    genome = next(iter(fasta(path).values()))
    out = {}
    for gene, (start, end) in H77.items():
        nuc = genome[start - 1 : end]
        out[gene] = {"nuc": nuc, "protein": protein(nuc)}
    return out


def unassessed_site(gene, pos, reference_codon, strand, reason):
    return {"gene": gene, "h77_position": pos, "sample_aa_position": None,
            "genomic_start": None, "genomic_end": None, "strand": strand,
            "h77_codon": reference_codon, "sample_codon": "",
            "h77_aa": next(iter(aas(reference_codon))), "amino_acids": [],
            "assessed": False, "assessment_reason": reason}


def make_sites(sequence, genes, h77, applicable, references=None, alignments=None):
    sites = {}
    for gene in sorted({r["region"] for r in applicable} & H77.keys()):
        positions = sorted({p["position"] for r in applicable if r["region"] == gene for p in parts(r["rule_definition"])})
        if gene not in genes:
            for pos in positions:
                sites[(gene, pos)] = unassessed_site(gene, pos, h77[gene]["nuc"][(pos-1)*3:pos*3], ".", "missing_coverage")
            continue
        info = genes[gene]
        if info["strand"] not in {"+", "-"} or not 1 <= info["start"] <= info["end"] <= len(sequence):
            raise ValueError(f"Invalid gene coordinates for {gene}")
        nuc = sequence[info["start"] - 1 : info["end"]]
        if info["strand"] == "-":
            nuc = rc(nuc)
        if not nuc:
            continue
        reference = (references or {}).get(gene, {"nuc": h77[gene]["nuc"], "accession": "AF011751.1", "h77_map": dict(enumerate(range(len(h77[gene]["nuc"]))))})
        projection = align(reference["nuc"], nuc)
        mapping = projection["mapping"]
        bridge = reference["h77_map"]
        unsafe = projection["unsafe_from"]
        warnings = []
        partial_start = info.get("attributes", {}).get("trunc5") == "1"
        query_phase = 0
        if mapping:
            initial = min(mapping)
            delta = (mapping[initial]-initial) % 3
            if partial_start:
                query_phase = delta
            elif delta:
                unsafe = min(unsafe, initial)
                warnings.append("Uncertain 5-prime coding boundary/phase; downstream codons are unassessed")
        if len(nuc) % 3:
            warnings.append("Non-triplet annotated interval; inspect terminal boundaries and internal indels")
        if info.get("attributes", {}).get("annotation_inferred") == "true":
            warnings.append("VADR inferred this gene from a misc_feature; boundaries require review")
        if projection["alignment_limit_exceeded"]:
            warnings.append("Alignment ambiguity exceeded the evaluation limit; affected interval and downstream sites are unassessed")
        if projection["unsafe_from"] < len(reference["nuc"]) and not projection["alignment_limit_exceeded"]:
            warnings.append(f"Unresolved internal frameshift from reference nucleotide {projection['unsafe_from']+1}")
        for offset in range(query_phase, len(nuc)-5, 3):
            observed = aas(nuc[offset:offset+3])
            if observed == {"*"}:
                unsafe = min(unsafe, min((p for p, q in mapping.items() if q >= offset), default=unsafe))
                warnings.append(f"Internal stop at annotated interval codon {(offset-query_phase)//3+1}")
            elif "*" in observed and "N" not in nuc[offset:offset+3]:
                warnings.append(f"Possible IUPAC stop at annotated interval codon {(offset-query_phase)//3+1}")
        if alignments is not None:
            alignments[gene] = {k: v for k, v in projection.items() if k not in {"mapping", "deletions"}}
            alignments[gene].update(reference=reference["accession"], warnings=warnings,
                                    unsafe_from=unsafe,
                                    genomic_start=info["start"], genomic_end=info["end"],
                                    annotation=info.get("attributes", {}),
                                    reference_translation=protein(reference["nuc"]),
                                    sample_translation=protein(nuc[query_phase:]), sample_translation_phase=query_phase)
        for pos in positions:
            numbered = (pos-1)*3
            reference_positions = [bridge.get(numbered+j) for j in range(3)]
            query_positions = [mapping.get(p) for p in reference_positions]
            reason = None
            deletion = False
            if None in reference_positions:
                reason = "numbering_reference_gap"
            elif min(reference_positions) >= unsafe or max(reference_positions) >= unsafe:
                reason = "coding_integrity"
            elif reference_positions != list(range(reference_positions[0], reference_positions[0]+3)):
                reason = "uncertain_codon_mapping"
            elif reference_positions[0] in projection["deletions"]:
                deletion = True
            elif None in query_positions:
                reason = "missing_coverage" if not mapping or reference_positions[-1] < min(mapping) or reference_positions[0] > max(mapping) else "ambiguous_alignment"
            elif query_positions != list(range(query_positions[0], query_positions[0]+3)):
                reason = "disrupted_codon"
            sample_pos = None if reason or deletion or info.get("attributes", {}).get("trunc5") == "1" else query_positions[0]//3+1
            ref_codon = h77[gene]["nuc"][(pos - 1) * 3 : pos * 3]
            ref_aa = next(iter(aas(ref_codon)))
            if reason:
                sites[(gene, pos)] = unassessed_site(gene, pos, ref_codon, info["strand"], reason)
                continue
            offset = mapping[reference_positions[0]-1] if deletion else query_positions[0]
            sample_codon = "---" if deletion else nuc[offset : offset + 3]
            if info["strand"] == "+":
                start = info["start"] + offset
                end = start if deletion else start + 2
            else:
                end = info["end"] - offset
                start = end if deletion else end - 2
            possible = ["del"] if deletion else sorted(aas(sample_codon))
            assessed = deletion or (
                len(sample_codon) == 3
                and "N" not in sample_codon
                and bool(possible)
                and not ({"X", "*"} & set(possible))
            )
            sites[(gene, pos)] = {
                "gene": gene,
                "h77_position": pos,
                "sample_aa_position": sample_pos,
                "sample_nt_offset": offset+1,
                "genomic_start": start,
                "genomic_end": end,
                "strand": info["strand"],
                "h77_codon": ref_codon,
                "sample_codon": sample_codon,
                "h77_aa": ref_aa,
                "amino_acids": possible,
                "assessed": assessed,
                "assessment_reason": None if assessed else ("missing_coverage" if "N" in sample_codon else "possible_stop_codon"),
                "deletion_anchor": deletion,
            }
    return sites


def read_codon(read, positions, min_bq):
    found = {}
    query_positions = {}
    qualities = read.query_qualities or []
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
        if (
            ref_pos in positions
            and query_pos is not None
            and query_pos < len(read.query_sequence or "")
            and query_pos < len(qualities)
            and qualities[query_pos] >= min_bq
        ):
            found[ref_pos] = read.query_sequence[query_pos].upper()
            query_positions[ref_pos] = query_pos
    if len(found) != 3 or any(found[p] not in "ACGT" for p in positions):
        return None
    ordered = [query_positions[p] for p in sorted(positions)]
    if ordered != list(range(ordered[0], ordered[0]+3)):
        return None
    return "".join(found[p] for p in positions)


def evidence(sites, cram, reference, min_bq, min_depth):
    if not cram:
        return
    if pysam is None:
        raise RuntimeError("pysam is required with --cram")
    with pysam.AlignmentFile(cram, "rc", reference_filename=reference) as handle:
        name, sequence = next(iter(fasta(reference).items()))
        if handle.references != (name,) or handle.lengths != (len(sequence),):
            raise ValueError("CRAM reference identity/length does not match analytical consensus")
        sq = handle.header.to_dict()["SQ"][0]
        if sq.get("M5") and sq["M5"].lower() != hashlib.md5(sequence.encode()).hexdigest():
            raise ValueError("CRAM reference checksum does not match analytical consensus")
        for site in sites.values():
            start, end = site["genomic_start"], site["genomic_end"]
            if start is None or site.get("deletion_anchor"):
                site["codon_evidence"] = {
                    "depth": 0,
                    "amino_acid_frequencies": {},
                    "assessed": False,
                }
                continue
            positions = list(range(start - 1, end))
            ordered = positions if site["strand"] == "+" else list(reversed(positions))
            templates = defaultdict(set)
            for read in handle.fetch(
                site.get("chrom", handle.references[0]), start - 1, end
            ):
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_qcfail
                    or read.is_duplicate
                ):
                    continue
                codon = read_codon(read, ordered, min_bq)
                if codon:
                    templates[read.query_name].add(
                        codon.translate(str.maketrans("ACGT", "TGCA")) if site["strand"] == "-" else codon
                    )
            counts = Counter(
                next(iter(value)) for value in templates.values() if len(value) == 1
            )
            aa_counts = Counter()
            for codon, count in counts.items():
                aa_counts[TABLE.get(codon, "X")] += count
            depth = sum(counts.values())
            site["codon_evidence"] = {
                "depth": depth,
                "codon_counts": dict(sorted(counts.items())),
                "amino_acid_frequencies": (
                    {
                        aa: round(count / depth, 6)
                        for aa, count in sorted(aa_counts.items())
                    }
                    if depth
                    else {}
                ),
                "assessed": depth >= min_depth,
            }


def indel_evidence(cram, reference, alignments, min_depth, min_bq):
    """Summarize read/assembly phase disagreement, not separate strain calls.

    Use net indel length in 30-nt windows so nearby alternative CIGAR placements
    contribute to the same signal. Conflicting mates are excluded.
    """
    if not cram:
        return
    with pysam.AlignmentFile(cram, "rc", reference_filename=reference) as handle:
        chrom = handle.references[0]
        for gene, alignment in alignments.items():
            start, end = alignment["genomic_start"]-1, alignment["genomic_end"]
            windows = defaultdict(lambda: defaultdict(set))
            for read in handle.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate or any(op == 3 for op, _ in read.cigartuples or []):
                    continue
                net, position, query_position = Counter(), read.reference_start, 0
                invalid = set()
                qualities = read.query_qualities or []
                for operation, length in read.cigartuples or []:
                    if operation in {1, 2}:
                        right = query_position+length if operation == 1 else query_position
                        if query_position < 1 or right >= len(qualities) or min(qualities[query_position-1:right+1]) < min_bq:
                            invalid.add(position//30)
                        else:
                            net[position//30] += length if operation == 1 else -length
                    if operation in {0, 2, 3, 7, 8}:
                        position += length
                    if operation in {0, 1, 4, 7, 8}:
                        query_position += length
                for bucket in range(max(start, read.reference_start)//30, min(end, read.reference_end)//30+1):
                    left, right = bucket*30, (bucket+1)*30
                    if bucket not in invalid and left >= start and right <= end and read.reference_start < left and read.reference_end > right:
                        windows[bucket][read.query_name].add(net[bucket])
            findings = []
            for bucket, templates in sorted(windows.items()):
                counts = Counter(next(iter(values)) for values in templates.values() if len(values) == 1)
                total = sum(counts.values())
                affected = sum(count for net, count in counts.items() if net % 3)
                if total >= min_depth and affected/total >= 0.15:
                    findings.append({"start": bucket*30+1, "end": (bucket+1)*30,
                                     "depth": total, "non_triplet_fraction": affected/total,
                                     "net_indel_counts": dict(sorted(counts.items()))})
            alignment["read_indel_evidence"] = findings
            if findings:
                alignment["warnings"].append("Reads support non-triplet indel differences from the assembly; inspect alignment and mixed-population signals (no automatic repair)")


def evaluate(applicable, sites, subtype):
    out = []
    for index, rule in enumerate(applicable, 1):
        rule_parts = parts(rule["rule_definition"])
        assessed = []
        matched = []
        reasons = []
        for part in rule_parts:
            site = sites.get((rule["region"], part["position"]))
            ok = bool(site and site["assessed"])
            assessed.append(ok)
            matched.append(ok and bool(matching_amino_acids(part, site)))
            reasons.append(None if ok else (site.get("assessment_reason", "missing_coverage") if site else "missing_coverage"))
        state = "not_assessed"
        if rule_parts and any(ok and not match for ok, match in zip(assessed, matched)):
            state = "none"
        elif rule_parts and all(assessed):
            state = "full" if all(matched) else "none"
        out.append(
            {
                "id": f"rule-{index}",
                **rule,
                "parts": rule_parts,
                "component_states": ["matched" if match else ("disproved" if ok else "unknown") for ok, match in zip(assessed, matched)],
                "component_reasons": reasons,
                "match_state": state,
                "licensed": licensed(
                    subtype, rule.get("drug_licensed_for_genotype", "")
                ),
            }
        )
    return out


def compound_warning(drug_rules):
    if len(drug_rules) != 1:
        return None
    rule = drug_rules[0]
    states = rule.get("component_states", [])
    missing = [reason for state, reason in zip(states, rule.get("component_reasons", [])) if state == "unknown"]
    if len(states) > 1 and "matched" in states and missing and set(missing) == {"missing_coverage"} and "disproved" not in states:
        return f"Incomplete compound rule: {rule['rule_definition']}; missing coverage at required positions"
    return None


def outcome(drug_rules, sites, drug_licensed, catalogue):
    if not drug_licensed:
        return "not licensed", None
    if not drug_rules:
        positions = {
            (r["region"], p["position"])
            for r in catalogue
            for p in parts(r["rule_definition"])
        }
        assessed = bool(positions) and all(
            sites.get(key, {}).get("assessed", False) for key in positions
        )
        return ("susceptible" if assessed else "not assessed"), None
    full = [r for r in drug_rules if r["match_state"] == "full"]
    if full:
        prediction = min(
            (r["prediction"] for r in full), key=lambda x: RANK.get(x.lower(), 9)
        )
        return prediction, prediction
    if any(r["match_state"] == "not_assessed" for r in drug_rules):
        return "not assessed", None
    changed = any(
        (site := sites.get((r["region"], p["position"])))
        and site["assessed"]
        and set(site["amino_acids"]) != {site["h77_aa"]}
        for r in drug_rules
        for p in r["parts"]
    )
    return ("susceptible with scored substitutions" if changed else "susceptible"), None


def escape(value):
    out = []
    for char in str(value):
        out.extend(
            [f"%{byte:02X}" for byte in char.encode()]
            if char in "%\t\n\r;=,&"
            else [char]
        )
    return "".join(out)


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def qualify_sites(sites, metadata):
    """Qualify analytical calls using the evidence that generated the FASTA."""
    for site in sites.values():
        site["nucleotide_depths"] = []
        site["low_depth"] = False
        if not metadata or site.get("genomic_start") is None:
            continue
        start, end = site["genomic_start"], site["genomic_end"]
        if site.get("deletion_anchor"):
            positions = [start-1, start+1]
        else:
            positions = list(range(start, end+1))
        depths = [metadata["usable_depths"][p-1] for p in positions if 1 <= p <= metadata["length"]]
        if site.get("strand") == "-":
            depths.reverse()
        site["nucleotide_depths"] = depths
        site["minimum_nucleotide_depth"] = min(depths) if depths else None
        site["depth_basis"] = "deletion_flanks" if site.get("deletion_anchor") else "consensus_allele_counts"
        site["low_depth"] = site["assessed"] and any(1 <= d <= 6 for d in depths)
        if site.get("assessment_reason") == "missing_coverage" and "N" in site.get("sample_codon", ""):
            reasons = [metadata["masking_reasons"].get(str(p)) for p in positions]
            if "all_four_alleles" in reasons and "missing_coverage" not in reasons:
                site["assessment_reason"] = "all_four_alleles"
            elif any(r and r != "missing_coverage" for r in reasons) and "missing_coverage" not in reasons:
                site["assessment_reason"] = "uninterpretable_consensus_evidence"


def write(output, sample, subtype, sites, evaluations, all_rules, rules_path, args):
    output.mkdir(parents=True, exist_ok=True)
    chrom = getattr(args, "track_chrom", next(iter(fasta(args.fasta))))
    matched = defaultdict(list)
    relevant = defaultdict(list)
    for rule in evaluations:
        for part in rule["parts"]:
            key = (rule["region"], part["position"])
            relevant[key].append(rule)
            if rule["match_state"] == "full" and matching_amino_acids(part, sites.get(key, {})):
                matched[key].append(rule)
    rows = []
    for key, site_rules in sorted(matched.items()):
        site = sites[key]
        supporting = {aa for rule in site_rules for part in rule["parts"]
                      if (rule["region"], part["position"]) == key
                      for aa in matching_amino_acids(part, site)}
        alt = ",".join(sorted(supporting))
        rows.append(
            {
                "sample": sample,
                "gene": site["gene"],
                "genomic_start": site["genomic_start"],
                "genomic_end": site["genomic_end"],
                "ref_nuc": site["h77_codon"],
                "alt_nuc": site["sample_codon"],
                "aa_pos": site["h77_position"],
                "ref_aa": site["h77_aa"],
                "alt_aa": alt,
                "rule_definition": "; ".join(
                    sorted({r["rule_definition"] for r in site_rules})
                ),
                "drugs": ", ".join(sorted({r["drug"] for r in site_rules})),
                "prediction": "; ".join(sorted({r["prediction"] for r in site_rules})),
                "reference": "; ".join(sorted({r["reference"] for r in site_rules})),
                "strand": site["strand"],
                "low_depth": site.get("low_depth", False),
                "nucleotide_depths": "/".join(map(str, site.get("nucleotide_depths", []))),
            }
        )
    fields = "sample gene genomic_start genomic_end ref_nuc alt_nuc aa_pos ref_aa alt_aa rule_definition drugs prediction reference strand low_depth nucleotide_depths".split()
    with open(
        output / f"{sample}_resistance.tsv", "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with open(output / f"{sample}_resistance.bed", "w", encoding="utf-8") as handle:
        handle.write("#chrom\tstart\tend\tname\tscore\tstrand\n")
        for r in rows:
            handle.write(
                f"{chrom}\t{r['genomic_start']-1}\t{r['genomic_end']}\t{r['gene']}:{r['aa_pos']}:{r['ref_aa']}>{r['alt_aa']}\t0\t{r['strand']}\n"
            )
    with open(output / f"{sample}_resistance.gff", "w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for r in rows:
            attrs = {
                "ID": f"{r['gene']}:{r['aa_pos']}",
                "Name": f"{r['gene']}:{r['ref_aa']}{r['aa_pos']}{r['alt_aa']}",
                "gene": r["gene"],
                "aa_pos": r["aa_pos"],
                "aa_change": f"{r['ref_aa']}>{r['alt_aa']}",
                "ref_codon": r["ref_nuc"],
                "sample_codon": r["alt_nuc"],
                "drugs": r["drugs"],
                "prediction": r["prediction"],
                "rule_definition": r["rule_definition"],
                "reference": r["reference"],
                "low_depth": r["low_depth"],
                "nucleotide_depths": r["nucleotide_depths"],
            }
            handle.write(
                f"{chrom}\tgeno2pheno\tresistance_mutation\t{r['genomic_start']}\t{r['genomic_end']}\t.\t{r['strand']}\t.\t"
                + ";".join(f"{k}={escape(v)}" for k, v in attrs.items())
                + "\n"
            )
    by_drug = defaultdict(list)
    for rule in evaluations:
        by_drug[rule["drug"]].append(rule)
    catalogue = defaultdict(list)
    for rule in all_rules:
        catalogue[rule["drug"]].append(rule)
    drug_outcomes = []
    with open(
        output / f"{sample}_resistance_by_drug.tsv", "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "drug",
                "outcome",
                "prediction",
                "licensed",
                "matched_rules",
                "partial_rules",
                "warnings",
            ]
        )
        for drug, drug_catalogue in sorted(catalogue.items()):
            drug_rules = by_drug[drug]
            drug_licensed = any(
                licensed(subtype, r.get("drug_licensed_for_genotype", ""))
                for r in drug_catalogue
            )
            state, prediction = outcome(drug_rules, sites, drug_licensed, drug_catalogue)
            scored = sorted({f"{s['gene']}:{s['h77_aa']}{s['h77_position']}{aa}"
                             for r in (drug_rules or drug_catalogue) for p in parts(r["rule_definition"])
                             if (s := sites.get((r["region"], p["position"]))) and s["assessed"]
                             for aa in s["amino_acids"] if aa != s["h77_aa"]})
            supporting = sorted({f"{r['region']}:{s['h77_aa']}{p['position']}{aa}"
                                 for r in drug_rules if r["match_state"] == "full" for p in r["parts"]
                                 if (s := sites.get((r["region"], p["position"])))
                                 for aa in matching_amino_acids(p, s)})
            record = {
                "drug": drug,
                "outcome": state,
                "prediction": prediction,
                "licensed": drug_licensed,
                "applicable_rule_count": len(drug_rules),
                "warnings": [compound_warning(drug_rules)] if drug_licensed and compound_warning(drug_rules) else [],
                "mutations": supporting if supporting else scored,
                "scored_substitutions": scored,
                "matched_rules": sorted(
                    {
                        r["rule_definition"]
                        for r in drug_rules
                        if r["match_state"] == "full"
                    }
                ),
                "partial_rules": sorted(
                    {
                        r["rule_definition"]
                        for r in drug_rules
                        if compound_warning(drug_rules) and r["match_state"] == "not_assessed"
                    }
                ),
            }
            low_sites = sorted({(r["region"], p["position"]) for r in (drug_rules or drug_catalogue)
                                for p in parts(r["rule_definition"])
                                if sites.get((r["region"], p["position"]), {}).get("low_depth")}) if drug_licensed else []
            record["low_depth_positions"] = [f"{gene}:{pos}" for gene, pos in low_sites]
            if low_sites:
                record["warnings"].append("Low-depth support (1-6 usable observations); curator review required: "+", ".join(f"{gene}:{pos} depths="+"/".join(map(str, sites[(gene,pos)]["nucleotide_depths"])) for gene,pos in low_sites))
            drug_outcomes.append(record)
            writer.writerow(
                [
                    drug,
                    state,
                    prediction or "",
                    str(record["licensed"]).lower(),
                    "; ".join(record["matched_rules"]),
                    "; ".join(record["partial_rules"]),
                    "; ".join(record["warnings"]),
                ]
            )
    with open(
        output / f"{sample}_resistance_sites.gff3", "w", encoding="utf-8"
    ) as handle:
        handle.write("##gff-version 3\n")
        for gene, alignment in getattr(args, "gene_alignments", {}).items():
            if alignment["warnings"]:
                attrs = {"ID": f"integrity:{gene}", "Name": f"{gene}:coding review", "warnings": "; ".join(alignment["warnings"])}
                handle.write(f"{chrom}\tVirPipa\tsequence_integrity\t{alignment['genomic_start']}\t{alignment['genomic_end']}\t.\t.\t.\t" + ";".join(f"{k}={escape(v)}" for k, v in attrs.items()) + "\n")
        for key, site in sorted(sites.items()):
            if site["genomic_start"] is None:
                continue
            rs = relevant[key]
            ev = site.get("codon_evidence", {})
            frequencies = ",".join(
                f"{aa}:{freq:.1%}"
                for aa, freq in ev.get("amino_acid_frequencies", {}).items()
            )
            state = (
                "not_assessed"
                if not site["assessed"]
                else ("resistance" if key in matched else ("scored_substitution" if set(site["amino_acids"]) != {site["h77_aa"]} else "reference"))
            )
            attrs = {
                "ID": f"site:{site['gene']}:{site['h77_position']}",
                "Name": f"{site['gene']}:{site['h77_aa']}{site['h77_position']}{'/'.join(site['amino_acids']) or '?'}",
                "gene": site["gene"],
                "h77_position": site["h77_position"],
                "sample_aa_position": site["sample_aa_position"],
                "h77_aa": site["h77_aa"],
                "observed_aa": ",".join(site["amino_acids"]),
                "resistance_aa": ",".join(sorted({aa for rule in matched.get(key, []) for part in rule["parts"] if (rule["region"], part["position"]) == key for aa in matching_amino_acids(part, site)})),
                "h77_codon": site["h77_codon"],
                "sample_codon": site["sample_codon"],
                "state": state,
                "assessment_reason": site.get("assessment_reason") or "",
                "deletion_anchor": site.get("deletion_anchor", False),
                "drugs": ", ".join(sorted({r["drug"] for r in rs})),
                "rules": "; ".join(sorted({r["rule_definition"] for r in rs})),
                "codon_depth": ev.get("depth", ""),
                "nucleotide_depths": "/".join(map(str, site.get("nucleotide_depths", []))),
                "low_depth": site.get("low_depth", False),
                "aa_frequencies": frequencies,
            }
            handle.write(
                f"{chrom}\tgeno2pheno\tresistance_site\t{site['genomic_start']}\t{site['genomic_end']}\t.\t{site['strand']}\t.\t"
                + ";".join(f"{k}={escape(v)}" for k, v in attrs.items())
                + "\n"
            )
    payload = {
        "schema_version": "3.1",
        "caller": "annotate_vcf_resistance.py",
        "caller_version": "3.1",
        "sample": sample,
        "subtype": subtype,
        "baseline_cutoff": 0.15,
        "baseline_source": "15% IUPAC consensus",
        "rules_sha256": checksum(rules_path),
        "h77_sha256": checksum(args.h77_fasta),
        "h77_reference": next(iter(fasta(args.h77_fasta))),
        "reference_panel_sha256": checksum(args.reference_panel),
        "consensus_provenance_verified": getattr(args, "consensus_provenance_verified", False),
        "exploratory_cutoff_range": {"minimum": 0.02, "maximum": 0.50},
        "codon_evidence_parameters": {
            "minimum_base_quality": args.minimum_base_quality,
            "minimum_depth": args.minimum_depth,
            "primary_alignments_only": True,
            "overlapping_mates_count_once": True,
            "mate_conflicts_discarded": True,
            "minimum_mapping_quality": 0,
            "indel_window_nt": 30,
            "indel_warning_fraction": 0.15,
        },
        "warnings": list(getattr(args, "consensus_warnings", [])),
        "gene_alignments": getattr(args, "gene_alignments", {}),
        "integrity_findings": [{"gene": gene, "message": message} for gene, alignment in getattr(args, "gene_alignments", {}).items() for message in alignment["warnings"]],
        "coordinate_reference": {"contig": chrom, "fasta_sha256": checksum(getattr(args, "track_reference_fasta", None) or args.fasta), "coordinate_source_sha256": getattr(args, "coordinate_source_sha256", None)},
        "consensus_sha256": checksum(args.fasta),
        "limitations": ["IUPAC mixtures do not establish linkage between compound-rule mutations", "No automatic assembly repair or strain separation"],
        "drug_outcomes": drug_outcomes,
        "sites": [site for _, site in sorted(sites.items())],
        "rules": [
            {k: v for k, v in rule.items() if k != "reference"} for rule in evaluations
        ],
    }
    for site in payload["sites"]:
        ev = site.get("codon_evidence")
        if ev and ev["assessed"]:
            called = {
                aa for aa, freq in ev["amino_acid_frequencies"].items() if freq >= 0.15
            }
            if set(site["amino_acids"]) != called:
                payload["warnings"].append(
                    f"15% IUPAC/codon evidence disagreement at {site['gene']} {site['h77_position']}: {site['amino_acids']} vs {sorted(called)}"
                )
    with open(output / f"{sample}_resistance.json", "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", "-f", required=True)
    parser.add_argument("--gff", "-g", required=True)
    parser.add_argument("--subtype", "-s", required=True)
    parser.add_argument("--rules", "-r", required=True)
    parser.add_argument(
        "--h77-fasta",
        default=str(Path(__file__).resolve().parents[1] / "refgenomes/resistance/AF011751.1.fasta"),
    )
    parser.add_argument("--cram")
    parser.add_argument("--consensus-metadata")
    parser.add_argument("--track-reference-metadata")
    parser.add_argument("--allow-legacy-consensus", action="store_true", help="Comparison only: accept unverified FASTA without a consensus sidecar")
    parser.add_argument("--track-reference-fasta", help="Coordinate-preserving main scaffold for combined IGV tracks")
    parser.add_argument("--reference-panel", default=str(Path(__file__).resolve().parents[1] / "refgenomes/resistance/panel.json"))
    parser.add_argument("--vcf", help="Deprecated and ignored")
    parser.add_argument("--sample-name")
    parser.add_argument("--output-dir", "-o", default="results")
    parser.add_argument("--minimum-base-quality", type=int, default=13)
    parser.add_argument("--minimum-depth", type=int, default=7)
    args = parser.parse_args()
    records = fasta(args.fasta)
    if len(records) != 1:
        raise ValueError("Resistance calling requires exactly one consensus contig")
    name, sequence = next(iter(records.items()))
    args.track_chrom = name
    metadata = load_consensus(args.fasta, args.consensus_metadata, required=not args.allow_legacy_consensus)
    if metadata and (metadata["mode"] != "iupac" or metadata["min_af"] != 0.15 or metadata["minimum_usable_depth"] != 1 or metadata["reference_fill"]):
        raise ValueError("Baseline resistance requires versioned 15% analytical IUPAC consensus")
    args.consensus_warnings = [] if metadata else ["Legacy comparison consensus: allele selection, coverage and low-depth support are unverified; regenerate before clinical interpretation"]
    args.consensus_provenance_verified = metadata is not None
    args.coordinate_source_sha256 = metadata.get("coordinate_source_sha256") if metadata else None
    if args.track_reference_fasta:
        track_records = fasta(args.track_reference_fasta)
        track_metadata = load_consensus(args.track_reference_fasta, args.track_reference_metadata)
        if len(track_records) != 1 or track_metadata["mode"] != "scaffold" or args.consensus_warnings:
            raise ValueError("Track rebasing requires versioned coordinate-preserving consensus and scaffold")
        if not re.fullmatch(r"[0-9a-f]{64}", args.coordinate_source_sha256 or "") or args.coordinate_source_sha256 != track_metadata.get("coordinate_source_sha256"):
            raise ValueError("Track rebasing requires a shared coordinate-source checksum")
        args.track_chrom, scaffold = next(iter(track_records.items()))
        if len(sequence) != len(scaffold):
            raise ValueError("Analytical consensus and track reference have different coordinate lengths")
    name = args.sample_name or name
    all_rules = rules(args.rules)
    applicable = [
        r for r in all_rules if subtype_match(args.subtype, r["subtype_pattern"])
    ]
    if not applicable:
        raise SystemExit(f"No geno2pheno rules apply to subtype {args.subtype}")
    h77 = h77_genes(args.h77_fasta)
    h77_records = fasta(args.h77_fasta)
    if len(h77_records) != 1 or next(iter(h77_records)) != "AF011751.1" or hashlib.sha256(next(iter(h77_records.values())).encode()).hexdigest() != "72fb77d5ba74c2d865e542b48e27f8662d3173752dbf217a6e5f08dc1537aa73":
        raise ValueError("Resistance numbering requires the pinned AF011751.1 sequence")
    references = reference_genes(args.reference_panel, args.subtype, h77, fasta)
    if args.subtype != "1a" and all(item["accession"] == "AF011751.1" for item in references.values()):
        args.consensus_warnings.append("No verified exact-subtype alignment bridge is available; using H77 directly and leaving uncertain positions unassessed")
    args.gene_alignments = {}
    sites = make_sites(sequence, gff(args.gff), h77, all_rules, references, args.gene_alignments)
    qualify_sites(sites, metadata)
    evidence(
        sites, args.cram, args.fasta, args.minimum_base_quality, args.minimum_depth
    )
    indel_evidence(args.cram, args.fasta, args.gene_alignments, args.minimum_depth, args.minimum_base_quality)
    evaluations = evaluate(applicable, sites, args.subtype)
    write(
        Path(args.output_dir), name, args.subtype, sites, evaluations, all_rules, args.rules, args
    )
    print(
        f"Called {sum(r['match_state']=='full' for r in evaluations)} matching rules across {len(sites)} resistance sites"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        raise SystemExit(1)
