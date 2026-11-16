#!/usr/bin/env python3
"""Prepare native-reference coverage and Virtitta trimming for cohort distances."""
import argparse
import json
from pathlib import Path
import sys

from benchmark_consensus import digest, write_json
from consensus_metadata import sequence


def depth_mask(path, identifier, length):
    mask = []
    for line in path.read_text().splitlines():
        chrom, position, depth = line.split('\t')
        if chrom != identifier or int(position) != len(mask)+1 or int(depth) < 0:
            raise ValueError(f'Invalid native-reference depth: {path}')
        mask.append(int(depth) > 0)
    if len(mask) != length:
        raise ValueError(f'Incomplete native-reference depth: {path}')
    return mask


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--native-root', type=Path, required=True)
    parser.add_argument('--virtitta-source', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--laptop-paths', action='store_true')
    parser.add_argument('--sample', action='append', default=[])
    args = parser.parse_args()
    sys.path.insert(0, str(args.virtitta_source))
    from virtitta.config import load_config
    from virtitta.distance import _trim_coordinates, parse_coverage_bed, distance_config_snapshot

    config = load_config(args.config)
    records = json.loads(args.inventory.read_text())
    if not (args.inventory.parent/'COMPLETE').is_file():
        raise ValueError('Frozen inventory is not complete')
    if args.sample:
        records = [r for r in records if r['sample_id'] in args.sample]
        if {r['sample_id'] for r in records} != set(args.sample):
            raise ValueError('Requested sample is absent from the inventory')
    args.output.mkdir(parents=True, exist_ok=False)
    rows, provenance, incompatible = [], [], []
    methods = {'baseline': ('baseline_iupac', 'baseline_majority'),
               'baseline_remap': ('baseline_iupac', 'baseline_majority'),
               'bcftools': ('bcftools_genotype', 'bcftools_genotype'),
               'ivar085': ('ivar_085', 'ivar_majority'),
               'ivar090': ('ivar_090', 'ivar_majority'),
               'awk_baqoff': ('awk_baqoff_iupac', 'awk_baqoff_majority')}
    for index, record in enumerate(records):
        directory = args.native_root/record['run_name']/record['sample_id']/'native'
        if not (directory/'COMPLETE').is_file():
            raise ValueError(f'Native experiment is not complete: {directory}')
        resistance = Path(record['paths']['resistance'])
        if args.laptop_paths and str(resistance).startswith('/fs1/'):
            resistance = Path('/mnt'+str(resistance))
        archived = json.loads(resistance.read_text())
        subtype = archived['subtype'] if (archived.get('typing_evidence') or {}).get('status') == 'assigned' else 'uncertain'
        identity = record['sample_id']+'_'+record['run_name']
        for method, arms in methods.items():
            sequences = [sequence(directory/(arm+'.fasta')) for arm in arms]
            coordinates = [_trim_coordinates(config, dna) for _, dna in sequences]
            # Match Virtitta's shared trimming only where its paired interface
            # is usable. Length-changing pairs remain available for pooled tests.
            compatible = len(sequences[0][1]) == len(sequences[1][1])
            if compatible:
                coordinates[1] = coordinates[0]
            else:
                incompatible.append(dict(sample_run_id=identity, method=method,
                    reason='Length-incompatible native IUPAC/majority; independently trimmed for pooled comparison'))
            for mode, arm, (identifier, dna), (start, end) in zip(
                    ('iupac', 'majority'), arms, sequences, coordinates):
                fasta_path = directory/(arm+'.fasta')
                if method == 'baseline':
                    coverage = Path(record['paths']['coverage'])
                    scaffold = Path(record['paths']['scaffold'])
                    if args.laptop_paths:
                        coverage = Path('/mnt'+str(coverage)) if str(coverage).startswith('/fs1/') else coverage
                        scaffold = Path('/mnt'+str(scaffold)) if str(scaffold).startswith('/fs1/') else scaffold
                    scaffold_id, scaffold_dna = sequence(scaffold)
                    if len(scaffold_dna) != len(dna):
                        raise ValueError('Baseline coverage and consensus lengths differ')
                    mask = parse_coverage_bed(coverage, scaffold_id, len(dna))
                else:
                    coverage = directory/(arm+'.depth.tsv')
                    mask = depth_mask(coverage, identifier, len(dna))
                token = f's{index}_{method}_{mode}'
                output_fasta, output_mask = args.output/(token+'.fasta'), args.output/(token+'.mask.json')
                output_fasta.write_text('>'+token+'\n'+dna[start:end]+'\n')
                write_json(output_mask, mask[start:end])
                rows.append(dict(sample_run_id=identity, method=method, mode=mode,
                    subtype=subtype, fasta=str(output_fasta.resolve()), mask=str(output_mask.resolve())))
                provenance.append(dict(token=token, fasta=str(fasta_path), fasta_sha256=digest(fasta_path),
                    coverage=str(coverage), coverage_sha256=digest(coverage),
                    coverage_policy='published scaffold >=1x' if method == 'baseline' else 'native remapping samtools depth -aa >=1x',
                    source_length=len(dna), trim_start=start, trim_end=end,
                    covered_after_trimming=sum(mask[start:end]), operational_pair_compatible=compatible))
    write_json(args.output/'manifest.json', rows)
    write_json(args.output/'preparation.json', dict(inventory_sha256=digest(args.inventory),
        config_sha256=digest(args.config), distance_config=distance_config_snapshot(config),
        virtitta_distance_sha256=digest(args.virtitta_source/'virtitta/distance.py'),
        limitations=['Equal length follows current Virtitta pairing rules; does not prove indel-coordinate equivalence',
                    'Uncertain typing is a diagnostic group, not a supported subtype cluster'],
        incompatible=incompatible, inputs=provenance))
    (args.output/'COMPLETE').touch()


if __name__ == '__main__':
    main()
