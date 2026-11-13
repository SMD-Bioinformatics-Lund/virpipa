#!/usr/bin/env python3
"""Freeze a cohort inventory and compare archived scaffold/IUPAC codon evidence.

Run directories are supplied explicitly, one per line. No recursive discovery of
unrelated campaigns, archive writes, database imports or marker publication.
"""
import argparse
import copy
import json
from pathlib import Path

from benchmark_consensus import digest, write_json
from annotate_vcf_resistance import evidence as production_evidence


def inventory(run_dirs, expected):
    samples = []
    identities = set()
    for directory in run_dirs:
        run = Path(directory)
        for qc in sorted(run.glob('*/*_qc_summary.json')):
            sample = qc.name.removesuffix('_qc_summary.json')
            payload = json.loads(qc.read_text())
            if payload['sample_id'] != sample or payload['run_name'] != run.name:
                raise ValueError(f'QC identity differs from archive path: {qc}')
            identity = (run.name, sample)
            if identity in identities:
                raise ValueError(f'Duplicate sample/run: {identity}')
            identities.add(identity)
            paths = {name: qc.parent/(sample+suffix) for name, suffix in {
                'scaffold': '.fasta', 'scaffold_cram': '.cram',
                'scaffold_index': '.cram.crai',
                'iupac': '-0.15-iupac.fasta', 'iupac_cram': '-0.15-iupac.cram',
                'iupac_index': '-0.15-iupac.cram.crai',
                'majority': '-majority.fasta', 'resistance': '_resistance.json',
                'coverage': '-coverage-1x.bed', 'qc': '_qc_summary.json'}.items()}
            missing = [str(p) for p in paths.values() if not p.is_file()]
            if missing:
                raise ValueError(f'Missing cohort inputs: {missing}')
            samples.append(dict(run_name=run.name, sample_id=sample,
                                paths={k: str(v.resolve()) for k, v in paths.items()},
                                sha256={k: digest(v) for k, v in paths.items()}))
    if len(samples) != expected:
        raise ValueError(f'Expected {expected} analyses, found {len(samples)}')
    return samples


def evidence(sample):
    paths = sample['paths']
    archived = json.loads(Path(paths['resistance']).read_text())
    eligible = []
    for site in archived['sites']:
        start, end = site.get('genomic_start'), site.get('genomic_end')
        if start is None or end is None or end-start != 2 or site.get('deletion_anchor'):
            continue
        eligible.append(site)
    modes = {}
    for mode in ('scaffold', 'iupac'):
        mapped = {index: copy.deepcopy(site) for index, site in enumerate(eligible)}
        production_evidence(mapped, paths[mode+'_cram'], paths[mode], 13, 7)
        modes[mode] = mapped
    sites = []
    for index, site in enumerate(eligible):
        results = {mode: modes[mode][index]['codon_evidence'] for mode in modes}
        if results['scaffold'] != results['iupac']:
            sites.append(dict(gene=site['gene'], position=site['h77_position'],
                              archived_assessed=site['assessed'], **results))
    return dict(run_name=sample['run_name'], sample_id=sample['sample_id'],
                changed_site_evidence=sites)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dirs', type=Path, required=True)
    parser.add_argument('--expected-samples', type=int, default=240)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--evidence', action='store_true', help='Compare complete codon evidence at catalogued sites')
    args = parser.parse_args()
    directories = [x.strip() for x in args.run_dirs.read_text().splitlines() if x.strip()]
    if len(directories) != len(set(directories)):
        raise ValueError('Duplicate run directories')
    args.output.mkdir(parents=True, exist_ok=False)
    samples = inventory(directories, args.expected_samples)
    write_json(args.output/'cohort.json', samples)
    if args.evidence:
        results = []
        for index, sample in enumerate(samples):
            results.append(evidence(sample))
            write_json(args.output/'codon-comparison.json', results)
            print(f'{index+1}/{len(samples)} {sample["sample_id"]}', flush=True)
    # Verify the archived files have not changed while measuring them.
    for sample in samples:
        for key, path in sample['paths'].items():
            if digest(path) != sample['sha256'][key]:
                raise ValueError(f'Input changed during benchmark: {path}')
    (args.output/'COMPLETE').touch()


if __name__ == '__main__':
    main()
