#!/usr/bin/env python3
"""
Parse geno2pheno HCV resistance rules into structured JSON format.
Handles subtype matching (1 matches all 1a,1b,1c; "1a,4a" matches either) and
compound rules (e.g., "445F and 451S" means ALL must be present).
"""

import json
import csv
import re
import argparse
import os
from collections import defaultdict


IUPAC_CODE = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def parse_subtype_pattern(pattern):
    """
    Parse subtype pattern and return list of subtypes it applies to.
    Examples:
        "1"      -> ["1a", "1b", "1c", ...] (all genotype 1 subtypes)
        "1a"     -> ["1a"]
        "1a, 4a" -> ["1a", "4a"]
        "3a"     -> ["3a"]
        "1,2,3,4,5,6" -> ["1a", "1b", "2a", "2b", ... all]
    """
    pattern = pattern.strip()
    
    if not pattern:
        return []
    
    subtypes = set()
    
    parts = [p.strip() for p in pattern.split(',')]
    
    for part in parts:
        if part == '1':
            subtypes.update([f'1{suffix}' for suffix in ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']])
        elif part == '2':
            subtypes.update([f'2{suffix}' for suffix in ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'k', 'l', 'm']])
        elif part == '3':
            subtypes.update([f'3{suffix}' for suffix in ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm']])
        elif part == '4':
            subtypes.update([f'4{suffix}' for suffix in ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r']])
        elif part == '5':
            subtypes.update([f'5{suffix}' for suffix in ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q']])
        elif part == '6':
            subtypes.update([f'6{suffix}' for suffix in ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 't', 'u', 'v', 'w', 'x']])
        elif part.startswith('1') or part.startswith('2') or part.startswith('3') or \
             part.startswith('4') or part.startswith('5') or part.startswith('6'):
            subtypes.add(part)
        else:
            subtypes.add(part)
    
    return list(subtypes)


def parse_rule_definition(rule_def):
    """
    Parse rule definition like "155K" or "41R and 80R" or "32 del" into list of dicts.
    Returns list of dicts: [{"position": 155, "aa": "K", "raw": "155K"}, ...]
    """
    rule_def = rule_def.strip()
    if not rule_def:
        return []
    
    rules = []
    parts = re.split(r'\s+and\s+', rule_def, flags=re.IGNORECASE)
    
    for part in parts:
        part = part.strip()
        match = re.match(r'(\d+)\s*(del|[A-Za-z*]+)?', part)
        if match:
            pos = int(match.group(1))
            aa = match.group(2) if match.group(2) else ''
            if aa == 'del':
                aa = 'del'
            rules.append({
                'position': pos,
                'aa': aa,
                'raw': part
            })
    
    return rules


def load_rules_csv(csv_file):
    """Load rules from CSV file."""
    rules = []
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    header = None
    data_lines = []
    for line in lines:
        if line.startswith('#'):
            continue
        if not header:
            reader = csv.reader([line])
            header = next(reader)
            continue
        data_lines.append(line)
    
    reader = csv.DictReader(data_lines, fieldnames=header)
    for row in reader:
        drugs_val = row.get('Drugs')
        if not drugs_val or drugs_val.startswith('#'):
            continue
            
        rule_defs = parse_rule_definition(row['Rule Definitions'])
        subtypes = parse_subtype_pattern(row['Rule for Subtype'])
        
        for rule in rule_defs:
            rules.append({
                'drug': row['Drugs'].strip(),
                'region': row['Region'].strip(),
                'position': rule['position'],
                'aa': rule['aa'],
                'rule_definition': row['Rule Definitions'].strip(),
                'subtypes': subtypes,
                'original_subtype': row['Rule for Subtype'].strip(),
                'drug_licensed': row['Drug licensed for genotype'].strip(),
                'prediction': row['Predictions'].strip(),
                'reference': row['Reference'].strip(),
                'is_compound': len(rule_defs) > 1
            })
    
    return rules


def build_rules_index(rules):
    """Build indexable structure from rules for fast lookup."""
    index = defaultdict(list)
    
    for rule in rules:
        key = (rule['region'], rule['position'], rule['aa'])
        index[key].append(rule)
    
    return index


def get_aa_from_codon(codon):
    """Get amino acid from codon, handling ambiguous bases."""
    if not codon or len(codon) != 3:
        return set()
    
    codon = codon.upper()
    
    bases = [IUPAC_CODE.get(b, [b]) for b in codon]
    
    possible_codons = ['', '', '', '']
    for b0 in bases[0]:
        for b1 in bases[1]:
            for b2 in bases[2]:
                full_codon = b0 + b1 + b2
                aa = CODON_TABLE.get(full_codon)
                if aa:
                    possible_codons.append(aa)
    
    return set(possible_codons)


def match_rules(region, aa_pos, alt_aa, subtype, rules_index):
    """
    Find matching rules for a variant.
    Returns list of matching rule dicts.
    """
    matches = []
    
    key = (region, aa_pos, alt_aa)
    
    for rule in rules_index.get(key, []):
        if subtype in rule['subtypes']:
            matches.append(rule)
    
    return matches


def check_compound_rule(found_variants, compound_rules):
    """
    Check if all required variants for a compound rule are present.
    found_variants: dict of {(region, pos): {aa, ...}}
    compound_rules: list of rule dicts that are part of compound rule
    Returns: list of rule dicts where all parts are found
    """
    if not compound_rules:
        return []
    
    grouped = defaultdict(list)
    for rule in compound_rules:
        rule_def = rule['rule_definition']
        grouped[rule_def].append(rule)
    
    results = []
    for rule_def, rule_parts in grouped.items():
        parsed = parse_rule_definition(rule_def)
        
        all_found = True
        for part in parsed:
            region = rule_parts[0]['region']
            var_key = (region, part['position'])
            if var_key not in found_variants:
                all_found = False
                break
            if part['aa'] not in found_variants[var_key]:
                all_found = False
                break
        
        if all_found:
            results.extend(rule_parts)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Parse geno2pheno HCV resistance rules into JSON'
    )
    parser.add_argument(
        '-i', '--input',
        default='hbv_result_rules.csv',
        help='Input CSV file (default: hbv_result_rules.csv)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output JSON file (default: <input>.json)'
    )
    parser.add_argument(
        '--subtype',
        help='Test subtype matching (e.g., 3a, 1b)'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        return 1
    
    print(f"Loading rules from {args.input}...")
    rules = load_rules_csv(args.input)
    print(f"Loaded {len(rules)} rule entries")
    
    rules_index = build_rules_index(rules)
    print(f"Built index with {len(rules_index)} unique (region, pos, aa) combinations")
    
    drugs = set(r['drug'] for r in rules)
    regions = set(r['region'] for r in rules)
    print(f"Drugs: {len(drugs)}")
    print(f"Regions: {regions}")
    
    if args.subtype:
        print(f"\nTesting subtype '{args.subtype}':")
        matching_rules = []
        for rule in rules:
            if args.subtype in rule['subtypes']:
                matching_rules.append(rule)
        print(f"Found {len(matching_rules)} rules matching subtype")
        
        by_drug = defaultdict(list)
        for rule in matching_rules:
            by_drug[rule['drug']].append(rule)
        
        for drug, drug_rules in sorted(by_drug.items()):
            print(f"  {drug}: {len(drug_rules)} rules")
    
    output_file = args.output
    if not output_file:
        output_file = args.input.replace('.csv', '.json')
    
    data = {
        'rules': rules,
        'index': {f"{k[0]}|{k[1]}|{k[2]}": len(v) for k, v in rules_index.items()},
        'metadata': {
            'total_rules': len(rules),
            'unique_positions': len(rules_index),
            'drugs': sorted(list(drugs)),
            'regions': sorted(list(regions))
        }
    }
    
    print(f"\nWriting to {output_file}...")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2)
    
    print("Done!")
    return 0


if __name__ == '__main__':
    exit(main())
