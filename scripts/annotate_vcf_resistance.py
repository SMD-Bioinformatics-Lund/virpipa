#!/usr/bin/env python3
"""
Annotate VCF variants with HCV drug resistance information from geno2pheno rules.

This script:
1. Parses VCF file to get variants
2. Uses VADR GFF to map genomic positions to genes (NS3, NS5A, NS5B, etc.)
3. Extracts codons from IUPAC FASTA at variant positions
4. Translates codons to amino acids (handling IUPAC ambiguity)
5. Matches against resistance rules for the given subtype
6. Outputs TSV results and BED file for IGV
"""

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def parse_vcf(vcf_file):
    """Parse VCF file using bcftools query."""
    variants = []
    
    cmd = [
        'bcftools', 'query',
        '-f', '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/DP\t%FORMAT/AD\t%INFO/AF\n',
        vcf_file
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running bcftools query: {e}\n")
        sys.stderr.write(f"stderr: {e.stderr}\n")
        return []
    
    for line in lines:
        if not line.strip():
            continue
        
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        
        chrom = parts[0]
        pos = int(parts[1])
        ref = parts[2]
        alt = parts[3]
        qual = float(parts[4]) if parts[4] != '.' else 0
        filter_status = parts[5] if parts[5] != '.' else 'PASS'
        dp = int(parts[6]) if parts[6] != '.' else None
        ad_str = parts[7] if parts[7] != '.' else None
        af_str = parts[8] if parts[8] != '.' else None
        
        freq = None
        if ad_str and dp and dp > 0:
            ad_vals = [int(x) for x in ad_str.split(',')]
            total_alt = sum(ad_vals[1:]) if len(ad_vals) > 1 else 0
            freq = total_alt / dp
        elif af_str:
            freq = float(af_str)
        
        for alt_allele in alt.split(','):
            variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt_allele,
                'qual': qual,
                'filter': filter_status,
                'dp': dp,
                'af': freq,
                'info': {}
            })
    
    return variants


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


def parse_gff(gff_file):
    """Parse VADR GFF file to get gene coordinates."""
    genes = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            
            if feature != 'gene':
                continue
            
            gene_id = None
            gene_name = None
            product = None
            
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    if key == 'ID':
                        gene_id = value
                    elif key == 'gene':
                        gene_name = value
                    elif key == 'product':
                        product = value
            
            if gene_name:
                genes[gene_name] = {
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'id': gene_id,
                    'product': product
                }
    
    return genes


def parse_fasta(fasta_file):
    """Parse FASTA file, returning sequence as string."""
    seq = []
    name = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is None:
                    name = line[1:].split()[0]
                continue
            seq.append(line.upper())
    
    return name, ''.join(seq)


def expand_iupac_codon(codon):
    """Expand ambiguous codon to list of all possible codons."""
    if not codon or len(codon) != 3:
        return []
    
    codon = codon.upper()
    bases = [IUPAC_CODE.get(b, [b]) for b in codon]
    
    codons = []
    for b0 in bases[0]:
        for b1 in bases[1]:
            for b2 in bases[2]:
                codons.append(b0 + b1 + b2)
    
    return codons


def translate_codon(codon):
    """Translate codon to amino acid."""
    return CODON_TABLE.get(codon.upper(), 'X')


def get_possible_aa(codon):
    """Get all possible amino acids from an ambiguous codon."""
    codons = expand_iupac_codon(codon)
    aa_set = set()
    for c in codons:
        aa = translate_codon(c)
        if aa:
            aa_set.add(aa)
    return aa_set


def find_gene_for_position(pos, genes):
    """Find which gene contains the given position."""
    for gene_name, gene_info in genes.items():
        if gene_info['start'] <= pos <= gene_info['end']:
            return gene_name, gene_info
    return None, None


def calculate_aa_position(genomic_pos, gene_start, gene_end, strand):
    """Calculate amino acid position from genomic position."""
    if strand == '+':
        offset = genomic_pos - gene_start
    else:
        offset = gene_end - genomic_pos
    
    aa_pos = (offset // 3) + 1
    codon_offset = offset % 3
    
    return aa_pos, codon_offset


def extract_codon_from_fasta(genomic_pos, gene_start, gene_end, fasta_seq, strand, codon_offset):
    """Extract codon from FASTA at the given genomic position."""
    codons = []
    
    if strand == '+':
        codon_start = genomic_pos - codon_offset
    else:
        codon_start = genomic_pos + codon_offset
    
    for i in range(3):
        pos = codon_start + i
        if strand == '-':
            pos = codon_start - i
        
        if pos < 1 or pos > len(fasta_seq):
            codons.append('N')
        else:
            codons.append(fasta_seq[pos - 1])
    
    ref_codon = ''.join(codons)
    
    return ref_codon, codon_start


def generate_alt_codons(ref_codon, ref_bases, alt_bases):
    """Generate alternative codons with the variant base(s) substituted."""
    alt_codons = set()
    
    for ref_base, alt_base in zip(ref_bases, alt_bases):
        if ref_base == alt_base:
            continue
        
        alt_codon_list = list(ref_codon)
        
        for i, (rb, ab) in enumerate(zip(ref_bases, alt_bases)):
            if rb != ab:
                alt_codon_list[i] = ab
        
        alt_codons.add(''.join(alt_codon_list))
    
    return list(alt_codons)


def parse_rule_definition(rule_def):
    """Parse rule definition into list of (position, amino_acid) tuples."""
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


def parse_subtype_pattern(pattern):
    """Parse subtype pattern and return list of subtypes it applies to."""
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
            import csv as csv_module
            reader = csv_module.reader([line])
            header = next(reader)
            continue
        data_lines.append(line)
    
    import csv as csv_module
    reader = csv_module.DictReader(data_lines, fieldnames=header)
    for row in reader:
        if not row.get('Drugs') or row['Drugs'].startswith('#'):
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
                'prediction': row['Predictions'].strip(),
                'reference': row['Reference'].strip(),
            })
    
    return rules


def build_rules_index(rules):
    """Build indexable structure from rules for fast lookup."""
    index = defaultdict(list)
    
    for rule in rules:
        key = (rule['region'], rule['position'], rule['aa'])
        index[key].append(rule)
    
    return index


def match_variant_to_rules(region, aa_pos, possible_aa, subtype, rules_index):
    """Match a variant against resistance rules."""
    matches = []
    
    for aa in possible_aa:
        key = (region, aa_pos, aa)
        for rule in rules_index.get(key, []):
            if subtype in rule['subtypes']:
                matches.append(rule)
    
    return matches


def main():
    parser = argparse.ArgumentParser(
        description='Annotate VCF variants with HCV drug resistance information'
    )
    parser.add_argument(
        '--vcf', '-v',
        required=True,
        help='Input VCF file'
    )
    parser.add_argument(
        '--gff', '-g',
        required=True,
        help='VADR GFF file with gene coordinates'
    )
    parser.add_argument(
        '--fasta', '-f',
        required=True,
        help='IUPAC FASTA file'
    )
    parser.add_argument(
        '--subtype', '-s',
        required=True,
        help='HCV subtype (e.g., 3a, 1b)'
    )
    parser.add_argument(
        '--rules', '-r',
        default='hbv_result_rules.csv',
        help='Rules CSV file (default: hbv_result_rules.csv)'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='results',
        help='Output directory (default: results subfolder of VCF location)'
    )
    parser.add_argument(
        '--sample-name',
        help='Sample name (default: derived from VCF)'
    )
    parser.add_argument(
        '--ref-bed',
        action='store_true',
        help='Also generate reference BED with all resistance positions from rules'
    )
    parser.add_argument(
        '--assets-dir',
        default='assets',
        help='Directory for reference files (default: assets)'
    )
    
    args = parser.parse_args()
    
    print(f"Loading rules from {args.rules}...")
    rules = load_rules_csv(args.rules)
    rules_index = build_rules_index(rules)
    print(f"Loaded {len(rules)} rules")
    
    print(f"Parsing GFF: {args.gff}")
    genes = parse_gff(args.gff)
    print(f"Found {len(genes)} genes: {', '.join(genes.keys())}")
    
    print(f"Parsing FASTA: {args.fasta}")
    fasta_name, fasta_seq = parse_fasta(args.fasta)
    print(f"FASTA length: {len(fasta_seq)} bp")
    
    sample_name = args.sample_name
    if not sample_name:
        sample_name = Path(args.vcf).stem.replace('.vcf', '').replace('.gz', '')
    
    vcf_parent = Path(args.vcf).parent
    sample_folder = vcf_parent.parent
    output_dir = sample_folder / args.output_dir
    
    assets_dir = Path(args.assets_dir)
    
    if args.ref_bed:
        ref_bed_file = assets_dir / "resistance_reference.bed"
        print(f"\nGenerating reference BED: {ref_bed_file}")
        assets_dir.mkdir(parents=True, exist_ok=True)
        with open(ref_bed_file, 'w') as f:
            f.write("#chrom\tstart\tend\tname\tscore\tstrand\tgene\tdrugs\tprediction\n")
            for rule in rules:
                gene = rule['region']
                pos = rule['position']
                aa = rule['aa']
                drug = rule['drug']
                pred = rule['prediction']
                strand = '+'
                if gene in genes:
                    gene_info = genes.get(gene)
                    if gene_info:
                        strand = gene_info.get('strand', '+')
                        gene_start = gene_info['start']
                        codon_start = gene_start + (pos - 1) * 3
                        codon_end = codon_start + 3
                        f.write(f"REF\t{codon_start-1}\t{codon_end}\t{gene}:{pos}:{aa}\t0\t{strand}\t{gene}\t{drug}\t{pred}\n")
        print(f"Reference BED written with {len(rules)} entries")
    
    print(f"Parsing VCF: {args.vcf}")
    variants = parse_vcf(args.vcf)
    print(f"Found {len(variants)} variants with ALT alleles")
    
    print(f"\nAnalyzing variants for subtype {args.subtype}...")
    
    results = []
    bed_entries = []
    
    for var in variants:
        chrom = var['chrom']
        pos = var['pos']
        ref = var['ref']
        alt = var['alt']
        
        gene_name, gene_info = find_gene_for_position(pos, genes)
        
        if not gene_name:
            continue
        
        if gene_name not in ['NS3', 'NS5A', 'NS5B', 'NS2', 'E1', 'E2', 'core', 'p7', 'NS4A', 'NS4B']:
            continue
        
        aa_pos, codon_offset = calculate_aa_position(
            pos, gene_info['start'], gene_info['end'], gene_info['strand']
        )
        
        ref_codon, codon_start = extract_codon_from_fasta(
            pos, gene_info['start'], gene_info['end'], fasta_seq,
            gene_info['strand'], codon_offset
        )
        
        codon_end = codon_start + 2
        
        possible_ref_aa = get_possible_aa(ref_codon)
        
        ref_base = ref[0] if ref else ''
        alt_base = alt[0] if alt else ''
        
        if len(ref) > 1 or len(alt) > 1:
            continue
        
        alt_codons = generate_alt_codons(ref_codon, ref_base, alt_base)
        possible_alt_aa = set()
        for ac in alt_codons:
            possible_alt_aa.update(get_possible_aa(ac))
        
        matches = match_variant_to_rules(
            gene_name, aa_pos, possible_alt_aa, args.subtype, rules_index
        )
        
        for match in matches:
            ref_aa_str = ','.join(sorted(possible_ref_aa)) if possible_ref_aa else '-'
            alt_aa_str = ','.join(sorted(possible_alt_aa)) if possible_alt_aa else '-'
            
            results.append({
                'sample': sample_name,
                'gene': gene_name,
                'genomic_pos': pos,
                'codon_start': codon_start,
                'codon_end': codon_end,
                'ref_nuc': ref,
                'alt_nuc': alt,
                'aa_pos': aa_pos,
                'ref_aa': ref_aa_str,
                'alt_aa': alt_aa_str,
                'rule_definition': match['rule_definition'],
                'drug': match['drug'],
                'prediction': match['prediction'],
                'reference': match['reference'],
                'strand': gene_info['strand']
            })
    
    print(f"Found {len(results)} raw rule matches")
    
    grouped = defaultdict(lambda: {
        'codon_start': None,
        'codon_end': None,
        'ref_nuc': '',
        'alt_nuc': '',
        'drugs': set(),
        'predictions': set(),
        'references': set(),
        'rule_definitions': set()
    })
    
    for r in results:
        key = (r['gene'], r['aa_pos'], r['ref_aa'], r['alt_aa'])
        if grouped[key]['codon_start'] is None:
            grouped[key]['codon_start'] = r['codon_start']
            grouped[key]['codon_end'] = r['codon_end']
        grouped[key]['ref_nuc'] = r['ref_nuc']
        grouped[key]['alt_nuc'] = r['alt_nuc']
        grouped[key]['drugs'].add(r['drug'])
        grouped[key]['predictions'].add(r['prediction'])
        grouped[key]['references'].add(r['reference'])
        grouped[key]['rule_definitions'].add(r['rule_definition'])
        grouped[key]['strand'] = r['strand']
    
    consolidated = []
    for (gene, aa_pos, ref_aa, alt_aa), data in sorted(grouped.items()):
        consolidated.append({
            'sample': sample_name,
            'gene': gene,
            'genomic_start': data['codon_start'],
            'genomic_end': data['codon_end'],
            'ref_nuc': data['ref_nuc'],
            'alt_nuc': data['alt_nuc'],
            'aa_pos': aa_pos,
            'ref_aa': ref_aa,
            'alt_aa': alt_aa,
            'rule_definition': '; '.join(sorted(data['rule_definitions'])),
            'drugs': ', '.join(sorted(data['drugs'])),
            'prediction': '; '.join(sorted(data['predictions'])),
            'reference': '; '.join(sorted(data['references'])),
            'strand': data['strand']
        })
    
    print(f"Consolidated to {len(consolidated)} unique amino acid changes")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tsv_file = output_dir / f"{sample_name}_resistance.tsv"
    print(f"\nWriting TSV: {tsv_file}")
    with open(tsv_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'sample', 'gene', 'genomic_start', 'genomic_end', 'ref_nuc', 'alt_nuc',
            'aa_pos', 'ref_aa', 'alt_aa', 'rule_definition',
            'drugs', 'prediction', 'reference', 'strand'
        ], delimiter='\t')
        writer.writeheader()
        for r in consolidated:
            writer.writerow(r)
    
    bed_file = output_dir / f"{sample_name}_resistance.bed"
    print(f"Writing BED: {bed_file}")
    with open(bed_file, 'w') as f:
        f.write("#chrom\tstart\tend\tname\tscore\tstrand\n")
        for r in consolidated:
            bed_name = f"{r['gene']}:{r['aa_pos']}:{r['ref_aa']}>{r['alt_aa']}"
            f.write(f"{sample_name}\t{r['genomic_start']-1}\t{r['genomic_end']}\t{bed_name}\t0\t{r['strand']}\n")
    
    drug_focused = defaultdict(list)
    for r in consolidated:
        for drug in r['drugs'].split(', '):
            drug_focused[drug.strip()].append(r)
    
    drug_tsv_file = output_dir / f"{sample_name}_resistance_by_drug.tsv"
    print(f"Writing drug-focused TSV: {drug_tsv_file}")
    with open(drug_tsv_file, 'w') as f:
        f.write(f"#Resistance analysis for {sample_name} (subtype: {args.subtype})\n")
        f.write(f"#Generated by annotate_vcf_resistance.py\n\n")
        
        for drug in sorted(drug_focused.keys()):
            f.write(f"## {drug}\n")
            writer = csv.DictWriter(f, fieldnames=[
                'sample', 'gene', 'genomic_start', 'genomic_end', 'ref_nuc', 'alt_nuc',
                'aa_pos', 'ref_aa', 'alt_aa', 'rule_definition',
                'prediction', 'reference'
            ], delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for r in drug_focused[drug]:
                writer.writerow(r)
            f.write("\n")
    
    print("\nDone!")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
