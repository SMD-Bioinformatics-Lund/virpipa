#!/usr/bin/env python3
"""Validate native benchmark artifacts and compare like-for-like read-codon calls."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import pysam
from benchmark_consensus import digest, write_json
from consensus_metadata import sequence


def calls(sites):
    return {(s['gene'], s['h77_position']): s for s in sites}


def compare(before, after):
    old, new = calls(before['sites']), calls(after['sites'])
    changes = []
    for key in sorted(old.keys() | new.keys()):
        a, b = old.get(key), new.get(key)
        fields = ('assessed', 'assessment_reason', 'amino_acids')
        if a is None or b is None or any(a[k] != b[k] for k in fields):
            changes.append(dict(gene=key[0], h77_position=key[1],
                before={k: a.get(k) for k in fields} if a else None,
                after={k: b.get(k) for k in fields} if b else None,
                before_counts=a.get('read_amino_acid_counts') if a else None,
                after_counts=b.get('read_amino_acid_counts') if b else None))
    old_drugs = {r['drug']: r for r in before['drug_outcomes']}
    new_drugs = {r['drug']: r for r in after['drug_outcomes']}
    drugs = [dict(drug=drug, before=old_drugs.get(drug), after=new_drugs.get(drug))
             for drug in sorted(old_drugs.keys() | new_drugs.keys())
             if any((old_drugs.get(drug) or {}).get(k) != (new_drugs.get(drug) or {}).get(k)
                    for k in ('outcome', 'prediction'))]
    return dict(site_changes=changes, drug_changes=drugs,
        changed_assessed_amino_acids=sum(s['assessed'] and new.get(k, {}).get('assessed', False)
            and s['amino_acids'] != new[k]['amino_acids'] for k, s in old.items()),
        lost_assessments=sum(s['assessed'] and not new.get(k, {}).get('assessed', False) for k, s in old.items()),
        gained_assessments=sum(s['assessed'] and not old.get(k, {}).get('assessed', False) for k, s in new.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--native-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--laptop-paths', action='store_true')
    parser.add_argument('--sample', action='append', default=[])
    args = parser.parse_args()
    if not (args.inventory.parent/'COMPLETE').is_file() or not (args.native_root.parent/'NATIVE_COMPLETE').is_file():
        raise ValueError('Inventory and native campaign must both be complete')
    records = json.loads(args.inventory.read_text())
    if args.sample:
        records = [r for r in records if r['sample_id'] in args.sample]
        if {r['sample_id'] for r in records} != set(args.sample):
            raise ValueError('Requested sample is absent from inventory')
    args.output.mkdir(parents=True, exist_ok=False)
    results, totals, containers = [], {}, set()
    for index, record in enumerate(records):
        directory = args.native_root/record['run_name']/record['sample_id']/'native'
        if not (directory/'COMPLETE').is_file():
            raise ValueError(f'Incomplete native experiment: {directory}')
        data = json.loads((directory/'scores.json').read_text())
        expected_hashes = sorted(record['sha256'][k] for k in ('scaffold', 'scaffold_cram', 'iupac', 'majority'))
        if sorted(line.split()[0] for line in (directory/'inputs.sha256').read_text().splitlines()) != expected_hashes:
            raise ValueError('Native task did not use the frozen baseline inputs')
        containers.add((directory/'containers.sha256').read_text())
        arms = {a['arm']: a for a in data['arms']}
        if len(arms) != 9 or 'baseline_iupac' not in arms:
            raise ValueError('Expected all nine native benchmark arms')
        summaries = []
        for name, arm in arms.items():
            fasta_path = directory/(name+'.fasta')
            identifier, dna = sequence(fasta_path)
            if arm['fasta_sha256'] != digest(fasta_path):
                raise ValueError('FASTA changed after scoring')
            with pysam.AlignmentFile(str(directory/(name+'.cram')), 'rc', reference_filename=str(fasta_path)) as handle:
                if handle.references != (identifier,) or handle.lengths != (len(dna),):
                    raise ValueError('Native CRAM reference identity/length mismatch')
                if handle.header.to_dict()['SQ'][0].get('M5') != hashlib.md5(dna.encode()).hexdigest():
                    raise ValueError('Native CRAM reference checksum mismatch')
                handle.check_index()
                next(handle.fetch(identifier, 0, min(len(dna), 120)), None)
            delta = compare(arms['baseline_iupac'], arm)
            counts = Counter(s.get('assessment_reason') or 'assessed' for s in arm['sites'])
            warnings = {gene: info['warnings'] for gene, info in arm['gene_alignments'].items() if info['warnings']}
            summaries.append(dict(arm=name, length=len(dna), unknown_bases=dna.count('N'),
                assessed=sum(s['assessed'] for s in arm['sites']), assessment_reasons=dict(counts),
                coordinate_projection=arm['coordinate_projection'], integrity_warnings=warnings, **delta))
            total = totals.setdefault(name, Counter())
            total.update(samples=1, samples_with_site_changes=bool(delta['site_changes']),
                         samples_with_drug_changes=bool(delta['drug_changes']),
                         site_changes=len(delta['site_changes']), drug_changes=len(delta['drug_changes']),
                         changed_assessed_amino_acids=delta['changed_assessed_amino_acids'],
                         lost_assessments=delta['lost_assessments'], gained_assessments=delta['gained_assessments'])
        resistance = Path(record['paths']['resistance'])
        if args.laptop_paths and str(resistance).startswith('/fs1/'):
            resistance = Path('/mnt'+str(resistance))
        production = json.loads(resistance.read_text())
        results.append(dict(run_name=record['run_name'], sample_id=record['sample_id'], arms=summaries,
            production_to_read_codons=compare(production, arms['baseline_iupac'])))
        print(f'{index+1}/{len(records)} validated', flush=True)
    write_json(args.output/'comparisons.json', results)
    write_json(args.output/'summary.json', dict(samples=len(records), totals=totals,
        baseline='baseline_iupac with observed read amino acids >=15%, minimum seven complete codons',
        inventory_sha256=digest(args.inventory), container_manifests_identical=len(containers) == 1,
        limitations=data['limitations']))
    if len(containers) != 1:
        raise ValueError('Container provenance differs between native tasks')
    (args.output/'COMPLETE').touch()


if __name__ == '__main__':
    main()
