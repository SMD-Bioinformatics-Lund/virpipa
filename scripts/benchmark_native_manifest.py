#!/usr/bin/env python3
"""Export native-tool inputs from an already frozen cohort inventory."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sample', action='append', default=[])
    parser.add_argument('--hopper-paths', action='store_true')
    args = parser.parse_args()
    if not (args.inventory.parent/'COMPLETE').is_file():
        raise ValueError('Inventory source-integrity verification is not complete')
    samples = json.loads(args.inventory.read_text())
    if args.sample:
        wanted = set(args.sample)
        samples = [s for s in samples if s['sample_id'] in wanted]
        if {s['sample_id'] for s in samples} != wanted:
            raise ValueError('Requested sample was not found in frozen inventory')
    fields = ['run_name', 'sample_id', 'scaffold', 'scaffold_cram', 'scaffold_index', 'iupac', 'majority', 'resistance']
    with args.output.open('x') as output:
        writer = csv.DictWriter(output, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for sample in samples:
            row = {k: sample[k] for k in fields[:2]}
            for key in fields[2:]:
                path = sample['paths'][key]
                if args.hopper_paths and path.startswith(('/mnt/fs1/', '/mnt/fs2/')):
                    path = path[4:]
                row[key] = path
            writer.writerow(row)


if __name__ == '__main__':
    main()
