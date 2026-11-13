#!/usr/bin/env python3
"""Experimental read-codon scorecard using existing numbering and rule code.

This writes benchmark evidence only, never a production resistance artifact.
Typing and catalogue are held constant. Unquantified deletions and existing
coding-integrity exclusions remain unassessed in the read-codon comparison.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

import annotate_vcf_resistance as caller
from benchmark_consensus import digest, write_json
from consensus_metadata import sequence
from resistance_alignment import align, reference_genes

ROOT = Path(__file__).resolve().parents[1]


def qualify_read_sites(sites, minimum=7):
    for site in sites.values():
        evidence = site.get('codon_evidence', {})
        counts = Counter()
        for codon, count in evidence.get('codon_counts', {}).items():
            counts[caller.TABLE[codon]] += count
        depth = sum(counts.values())
        site['read_amino_acid_counts'] = dict(counts)
        site['sequence_amino_acids'] = site['amino_acids']
        site['sequence_assessed'] = site['assessed']
        eligible = (site['assessed'] or site.get('assessment_reason') == 'missing_coverage')
        if site.get('deletion_anchor'):
            site['assessed'] = False
            site['assessment_reason'] = 'deletion_frequency_not_quantified'
        elif depth < minimum:
            site['assessed'] = False
            site['assessment_reason'] = 'insufficient_complete_codons'
        elif eligible and site.get('genomic_start') is not None:
            called = sorted(aa for aa, count in counts.items() if count*100 >= 15*depth)
            site['amino_acids'] = called
            site['assessed'] = '*' not in called
            site['assessment_reason'] = None if site['assessed'] else 'qualifying_stop_codon'


def score(native_dir, archived, arm, all_rules, h77, references):
    ref = native_dir/(arm+'.fasta')
    _, scaffold = sequence(native_dir/'scaffold_control.fasta')
    _, dna = sequence(ref)
    if arm in ('scaffold_control', 'baseline_iupac', 'baseline_majority',
               'awk_baqoff_iupac', 'awk_baqoff_majority'):
        if len(scaffold) != len(dna):
            raise ValueError('Coordinate-preserving control changed length')
        projection = dict(mapping=dict(enumerate(range(len(dna)))), alignment_limit_exceeded=False)
    else:
        projection = align(scaffold, dna)
    genes = {}
    for gene, info in archived['gene_alignments'].items():
        start, end = info['genomic_start']-1, info['genomic_end']-1
        mapping = projection['mapping']
        strands = {s['strand'] for s in archived['sites'] if s['gene'] == gene and s.get('strand') in ('+', '-')}
        if start in mapping and end in mapping and mapping[start] <= mapping[end] and len(strands) == 1:
            genes[gene] = dict(start=mapping[start]+1, end=mapping[end]+1,
                               strand=next(iter(strands)), attributes=info['annotation'])
    alignments = {}
    sites = caller.make_sites(dna, genes, h77, all_rules, references, alignments)
    caller.evidence(sites, str(native_dir/(arm+'.cram')), str(ref), 13, 7)
    caller.indel_evidence(str(native_dir/(arm+'.cram')), str(ref), alignments, 7, 13)
    qualify_read_sites(sites)
    subtype = archived['subtype'] or ''
    applicable = [r for r in all_rules if caller.subtype_match(subtype, r['subtype_pattern'])]
    evaluations = caller.evaluate(applicable, sites, subtype)
    outcomes = []
    typing_status = (archived.get('typing_evidence') or {}).get('status')
    for drug in sorted({r['drug'] for r in all_rules}):
        catalogue = [r for r in all_rules if r['drug'] == drug]
        rules = [r for r in evaluations if r['drug'] == drug]
        licensed = any(caller.licensed(subtype, r.get('drug_licensed_for_genotype', '')) for r in catalogue)
        state, prediction = caller.outcome(rules, sites, licensed, catalogue)
        if typing_status in ('ambiguous', 'unassigned') or (typing_status == 'genotype_only' and licensed and not rules):
            state, prediction = 'not assessed', None
        outcomes.append(dict(drug=drug, outcome=state, prediction=prediction))
    return dict(arm=arm, length=len(dna), unknown_bases=dna.count('N'), fasta_sha256=digest(ref),
                coordinate_projection=dict(mapped_bases=len(projection['mapping']),
                    alignment_limit_exceeded=projection['alignment_limit_exceeded'], genes_retained=sorted(genes)),
                gene_alignments=alignments, sites=list(sites.values()), drug_outcomes=outcomes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-dir', type=Path, required=True)
    parser.add_argument('--archived-resistance', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    archived = json.loads(args.archived_resistance.read_text())
    rules_path = ROOT/'assets/hcv_geno2pheno_rules.csv'
    all_rules = caller.rules(rules_path)
    h77_path = ROOT/'refgenomes/resistance/AF011751.1.fasta'
    h77 = caller.h77_genes(h77_path)
    references = reference_genes(ROOT/'refgenomes/resistance/panel.json', archived['subtype'] or '', h77, caller.fasta)
    results = [score(args.native_dir, archived, arm, all_rules, h77, references) for arm in
        ('scaffold_control', 'baseline_iupac', 'baseline_majority', 'bcftools_genotype',
         'ivar_majority', 'ivar_085', 'ivar_090', 'awk_baqoff_iupac', 'awk_baqoff_majority')]
    write_json(args.output, dict(experimental=True, baseline_source='observed complete read codons >=15%',
        minimum_depth=7, minimum_base_quality=13, typing_held_constant=archived.get('typing_evidence'),
        rules_sha256=digest(rules_path), h77_sha256=digest(h77_path),
        archived_resistance_sha256=digest(args.archived_resistance),
        limitations=['Archived VADR boundaries projected conservatively; no fresh VADR gene annotation',
                     'Deletion frequencies are not quantified by the existing read-codon extractor',
                     'Compound-rule mutations are not phased across distant codons',
                     'Low-depth observations remain in counts but are not assessed below seven complete codons'],
        arms=results))


if __name__ == '__main__':
    main()
