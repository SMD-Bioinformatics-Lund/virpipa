#!/usr/bin/env python3
"""Package completed known-molecule probes for the native Nextflow workflow."""
import argparse
import csv
import json
from pathlib import Path
import shutil

import pysam
from benchmark_consensus import write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probes', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--execution-root', type=Path, required=True,
                        help='Absolute package path as seen by the workflow')
    args = parser.parse_args()
    results = json.loads((args.probes/'results.json').read_text())
    cases = json.loads((args.probes/'fixtures/fixtures.json').read_text())
    if {r['case'] for r in results if r['truth_pass']} != {c['case'] for c in cases}:
        raise ValueError('All input probes must pass known-molecule extraction')
    args.output.mkdir(parents=True, exist_ok=False)
    fields = ['run_name', 'sample_id', 'scaffold', 'scaffold_cram', 'scaffold_index',
              'iupac', 'majority', 'resistance']
    with (args.output/'inputs.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for case in cases:
            name = case['case']
            directory = args.output/name
            directory.mkdir()
            fasta = directory/'reference.fasta'
            shutil.copyfile(case['reference'], fasta)
            pysam.faidx(str(fasta))
            cram = directory/'source.cram'
            with pysam.AlignmentFile(case['bam']) as reads:
                with pysam.AlignmentFile(str(cram), 'wc', header=reads.header,
                        reference_filename=str(fasta), format_options=[b'embed_ref=1']) as output:
                    for read in reads.fetch(until_eof=True):
                        output.write(read)
            pysam.index(str(cram))
            for source, target in [('baseline.fasta', 'iupac.fasta'),
                                   ('baseline-majority.fasta', 'majority.fasta')]:
                shutil.copyfile(args.probes/'results'/name/source, directory/target)
            write_json(directory/'resistance.json', dict(subtype='1a',
                typing_evidence={'status': 'assigned'}, gene_alignments={}, sites=[]))
            write_json(directory/'truth.json', case)
            paths = dict(scaffold='reference.fasta', scaffold_cram='source.cram',
                scaffold_index='source.cram.crai', iupac='iupac.fasta',
                majority='majority.fasta', resistance='resistance.json')
            writer.writerow(dict(run_name='SYNTHETIC', sample_id=name,
                **{k: str(args.execution_root/name/v) for k, v in paths.items()}))
    write_json(args.output/'truth.json', cases)
    (args.output/'COMPLETE').touch()


if __name__ == '__main__':
    main()
