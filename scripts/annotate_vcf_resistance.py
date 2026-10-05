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
    for raw in open(path, encoding="utf-8"):
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
    for raw in open(path, encoding="utf-8"):
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
            out[attrs["gene"]] = {
                "chrom": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
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
        result.append(choices[0] if choices else "X")
    return "".join(result)


def position_map(reference, query):
    """Needleman-Wunsch map from 1-based H77 residues to sample residues."""
    rows, cols = len(reference) + 1, len(query) + 1
    score = [[0] * cols for _ in range(rows)]
    trace = [[0] * cols for _ in range(rows)]
    for i in range(1, rows):
        score[i][0] = -2 * i
        trace[i][0] = 1
    for j in range(1, cols):
        score[0][j] = -2 * j
        trace[0][j] = 2
    for i in range(1, rows):
        for j in range(1, cols):
            values = (
                score[i - 1][j - 1] + (2 if reference[i - 1] == query[j - 1] else -1),
                score[i - 1][j] - 2,
                score[i][j - 1] - 2,
            )
            score[i][j] = max(values)
            trace[i][j] = values.index(score[i][j])
    result = {}
    i, j = len(reference), len(query)
    while i or j:
        direction = trace[i][j]
        if i and j and direction == 0:
            result[i] = j
            i -= 1
            j -= 1
        elif i and (not j or direction == 1):
            result[i] = None
            i -= 1
        else:
            j -= 1
    return result


def parts(definition):
    result = []
    for raw in re.split(r"\s+and\s+", (definition or "").strip(), flags=re.I):
        match = re.fullmatch(r"(\d+)\s*(del|[A-Za-z*]+)?", raw.strip())
        if match:
            result.append(
                {
                    "position": int(match.group(1)),
                    "aa": match.group(2) or "",
                    "raw": raw.strip(),
                }
            )
    return result


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


def make_sites(sequence, genes, h77, applicable):
    sites = {}
    for gene in sorted({r["region"] for r in applicable} & H77.keys()):
        if gene not in genes:
            continue
        info = genes[gene]
        nuc = sequence[info["start"] - 1 : info["end"]]
        if info["strand"] == "-":
            nuc = rc(nuc)
        mapping = position_map(h77[gene]["protein"], protein(nuc))
        positions = sorted(
            {
                p["position"]
                for r in applicable
                if r["region"] == gene
                for p in parts(r["rule_definition"])
            }
        )
        for pos in positions:
            sample_pos = mapping.get(pos)
            ref_codon = h77[gene]["nuc"][(pos - 1) * 3 : pos * 3]
            ref_aa = next(iter(aas(ref_codon)))
            if sample_pos is None:
                sites[(gene, pos)] = {
                    "gene": gene,
                    "h77_position": pos,
                    "sample_aa_position": None,
                    "genomic_start": None,
                    "genomic_end": None,
                    "strand": info["strand"],
                    "h77_codon": ref_codon,
                    "sample_codon": "---",
                    "h77_aa": ref_aa,
                    "amino_acids": ["del"],
                    "assessed": True,
                }
                continue
            offset = (sample_pos - 1) * 3
            sample_codon = nuc[offset : offset + 3]
            if info["strand"] == "+":
                start = info["start"] + offset
                end = start + 2
            else:
                end = info["end"] - offset
                start = end - 2
            possible = sorted(aas(sample_codon))
            assessed = (
                len(sample_codon) == 3
                and bool(possible)
                and not ({"X", "*"} & set(possible))
            )
            sites[(gene, pos)] = {
                "gene": gene,
                "h77_position": pos,
                "sample_aa_position": sample_pos,
                "genomic_start": start,
                "genomic_end": end,
                "strand": info["strand"],
                "h77_codon": ref_codon,
                "sample_codon": sample_codon,
                "h77_aa": ref_aa,
                "amino_acids": possible,
                "assessed": assessed,
            }
    return sites


def read_codon(read, positions, min_bq):
    found = {}
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
    if len(found) != 3 or any(found[p] not in "ACGT" for p in positions):
        return None
    return "".join(found[p] for p in positions)


def evidence(sites, cram, reference, min_bq, min_depth):
    if not cram:
        return
    if pysam is None:
        raise RuntimeError("pysam is required with --cram")
    with pysam.AlignmentFile(cram, "rc", reference_filename=reference) as handle:
        for site in sites.values():
            start, end = site["genomic_start"], site["genomic_end"]
            if start is None:
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
                        rc(codon) if site["strand"] == "-" else codon
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


def evaluate(applicable, sites, subtype):
    out = []
    for index, rule in enumerate(applicable, 1):
        rule_parts = parts(rule["rule_definition"])
        assessed = []
        matched = []
        for part in rule_parts:
            site = sites.get((rule["region"], part["position"]))
            ok = bool(site and site["assessed"])
            assessed.append(ok)
            matched.append(ok and part["aa"] in site["amino_acids"])
        state = "not_assessed"
        if rule_parts and all(assessed):
            state = (
                "full"
                if all(matched)
                else ("partial" if len(rule_parts) > 1 and any(matched) else "none")
            )
        out.append(
            {
                "id": f"rule-{index}",
                **rule,
                "parts": rule_parts,
                "match_state": state,
                "licensed": licensed(
                    subtype, rule.get("drug_licensed_for_genotype", "")
                ),
            }
        )
    return out


def outcome(drug_rules, sites):
    if not any(r["licensed"] for r in drug_rules):
        return "not licensed", None
    full = [r for r in drug_rules if r["match_state"] == "full"]
    if full:
        prediction = min(
            (r["prediction"] for r in full), key=lambda x: RANK.get(x.lower(), 9)
        )
        return prediction, prediction
    if not any(r["match_state"] != "not_assessed" for r in drug_rules):
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


def write(output, sample, subtype, sites, evaluations, rules_path, args):
    output.mkdir(parents=True, exist_ok=True)
    matched = defaultdict(list)
    relevant = defaultdict(list)
    for rule in evaluations:
        for part in rule["parts"]:
            key = (rule["region"], part["position"])
            relevant[key].append(rule)
            if rule["match_state"] == "full" and part["aa"] in sites.get(key, {}).get(
                "amino_acids", []
            ):
                matched[key].append(rule)
    rows = []
    for key, site_rules in sorted(matched.items()):
        site = sites[key]
        alt = (
            ",".join(sorted(set(site["amino_acids"]) - {site["h77_aa"]}))
            or site["h77_aa"]
        )
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
            }
        )
    fields = "sample gene genomic_start genomic_end ref_nuc alt_nuc aa_pos ref_aa alt_aa rule_definition drugs prediction reference strand".split()
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
                f"{sample}\t{r['genomic_start']-1}\t{r['genomic_end']}\t{r['gene']}:{r['aa_pos']}:{r['ref_aa']}>{r['alt_aa']}\t0\t{r['strand']}\n"
            )
    with open(output / f"{sample}_resistance.gff", "w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for r in rows:
            attrs = {
                "ID": f"{r['gene']}:{r['aa_pos']}",
                "gene": r["gene"],
                "aa_pos": r["aa_pos"],
                "aa_change": f"{r['ref_aa']}>{r['alt_aa']}",
                "ref_codon": r["ref_nuc"],
                "sample_codon": r["alt_nuc"],
                "drugs": r["drugs"],
                "prediction": r["prediction"],
                "rule_definition": r["rule_definition"],
                "reference": r["reference"],
            }
            handle.write(
                f"{sample}\tgeno2pheno\tresistance_mutation\t{r['genomic_start']}\t{r['genomic_end']}\t.\t{r['strand']}\t.\t"
                + ";".join(f"{k}={escape(v)}" for k, v in attrs.items())
                + "\n"
            )
    by_drug = defaultdict(list)
    for rule in evaluations:
        by_drug[rule["drug"]].append(rule)
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
            ]
        )
        for drug, drug_rules in sorted(by_drug.items()):
            state, prediction = outcome(drug_rules, sites)
            record = {
                "drug": drug,
                "outcome": state,
                "prediction": prediction,
                "licensed": any(r["licensed"] for r in drug_rules),
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
                        if r["match_state"] == "partial"
                    }
                ),
            }
            drug_outcomes.append(record)
            writer.writerow(
                [
                    drug,
                    state,
                    prediction or "",
                    str(record["licensed"]).lower(),
                    "; ".join(record["matched_rules"]),
                    "; ".join(record["partial_rules"]),
                ]
            )
    with open(
        output / f"{sample}_resistance_sites.gff3", "w", encoding="utf-8"
    ) as handle:
        handle.write("##gff-version 3\n")
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
                else ("resistance" if key in matched else "observed")
            )
            attrs = {
                "ID": f"site:{site['gene']}:{site['h77_position']}",
                "gene": site["gene"],
                "h77_position": site["h77_position"],
                "sample_aa_position": site["sample_aa_position"],
                "h77_aa": site["h77_aa"],
                "observed_aa": ",".join(site["amino_acids"]),
                "h77_codon": site["h77_codon"],
                "sample_codon": site["sample_codon"],
                "state": state,
                "drugs": ", ".join(sorted({r["drug"] for r in rs})),
                "rules": "; ".join(sorted({r["rule_definition"] for r in rs})),
                "codon_depth": ev.get("depth", ""),
                "aa_frequencies": frequencies,
            }
            handle.write(
                f"{sample}\tgeno2pheno\tresistance_site\t{site['genomic_start']}\t{site['genomic_end']}\t.\t{site['strand']}\t.\t"
                + ";".join(f"{k}={escape(v)}" for k, v in attrs.items())
                + "\n"
            )
    payload = {
        "schema_version": "2.0",
        "caller": "annotate_vcf_resistance.py",
        "caller_version": "2.0",
        "sample": sample,
        "subtype": subtype,
        "baseline_cutoff": 0.15,
        "baseline_source": "15% IUPAC consensus",
        "rules_sha256": checksum(rules_path),
        "h77_sha256": checksum(args.h77_fasta),
        "exploratory_cutoff_range": {"minimum": 0.02, "maximum": 0.50},
        "codon_evidence_parameters": {
            "minimum_base_quality": args.minimum_base_quality,
            "minimum_depth": args.minimum_depth,
            "primary_alignments_only": True,
            "overlapping_mates_count_once": True,
            "mate_conflicts_discarded": True,
        },
        "warnings": [],
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
    parser.add_argument("--h77-fasta", required=True)
    parser.add_argument("--cram")
    parser.add_argument("--vcf", help="Deprecated and ignored")
    parser.add_argument("--sample-name")
    parser.add_argument("--output-dir", "-o", default="results")
    parser.add_argument("--minimum-base-quality", type=int, default=13)
    parser.add_argument("--minimum-depth", type=int, default=7)
    args = parser.parse_args()
    records = fasta(args.fasta)
    name, sequence = next(iter(records.items()))
    name = args.sample_name or name
    all_rules = rules(args.rules)
    applicable = [
        r for r in all_rules if subtype_match(args.subtype, r["subtype_pattern"])
    ]
    if not applicable:
        raise SystemExit(f"No geno2pheno rules apply to subtype {args.subtype}")
    sites = make_sites(sequence, gff(args.gff), h77_genes(args.h77_fasta), all_rules)
    evidence(
        sites, args.cram, args.fasta, args.minimum_base_quality, args.minimum_depth
    )
    evaluations = evaluate(applicable, sites, args.subtype)
    write(
        Path(args.output_dir), name, args.subtype, sites, evaluations, args.rules, args
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
