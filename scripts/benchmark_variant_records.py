#!/usr/bin/env python3
"""Probe equivalent indel representations and deliberately conflicting records."""
import argparse
import json
from pathlib import Path

from benchmark_consensus import Runner, digest, write_json
from consensus_metadata import sequence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tools', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True,
                        help='Reference from benchmark_consensus.py fixtures')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    identifier, dna = sequence(args.reference)
    if dna[90:101] != 'GCCCCCCCCCC':
        raise ValueError('Expected the defined homopolymer fixture')
    ref = args.output/'reference.fasta'
    ref.write_text('>'+identifier+'\n'+dna+'\n')
    runner = Runner(json.loads(args.tools.read_text()), args.output)
    runner.run('samtools', ['faidx', ref])
    expected = dna[:91]+dna[92:]
    # Each valid record removes one C from the same run of ten Cs.
    records = {'left': [(91, 'GC', 'G')], 'right': [(100, 'CC', 'C')],
               'padded': [(91, dna[90:101], dna[90:100])],
               # The SNP refers to the deleted base: intentionally inconsistent.
               'conflicting_overlap': [(91, 'GC', 'G'), (92, 'C', 'T')]}
    header = (f'##fileformat=VCFv4.2\n##contig=<ID={identifier},length={len(dna)}>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruth\n')
    results = []
    for name, variants in records.items():
        directory = args.output/name
        directory.mkdir()
        vcf = directory/'raw.vcf'
        vcf.write_text(header+''.join(f'{identifier}\t{pos}\t.\t{a}\t{b}\t100\tPASS\t.\tGT:AD\t1:0,100\n'
                                     for pos, a, b in variants))
        normalized = directory/'normalized.vcf.gz'
        runner.run('bcftools', ['norm', '-f', ref, '-Oz', '-o', normalized, vcf])
        runner.run('bcftools', ['index', normalized])
        runner.run('bcftools', ['view', normalized], stdout=directory/'normalized.vcf')
        log = args.output/'tools.stderr.log'
        previous_log_length = log.stat().st_size
        runner.run('bcftools', ['consensus', '-f', ref, '-s', 'truth', '-H', '1', normalized],
                   stdout=directory/'bcftools.fasta')
        warnings = log.read_bytes()[previous_log_length:].decode()
        outputs = {path.name: dict(length=len(sequence(path)[1]), unknown=sequence(path)[1].count('N'),
                    matches_known_deletion=sequence(path)[1] == expected) for path in directory.glob('*.fasta')}
        results.append(dict(case=name, consistent_callset=name != 'conflicting_overlap',
                            tool_messages=warnings, outputs=outputs))
    write_json(args.output/'results.json', dict(reference_sha256=digest(ref), cases=results,
        bcftools_version=runner.run('bcftools', ['--version']).decode(),
        interpretation='Conflicting overlap has no consistent haploid truth; inspect tool warnings, never silently merge records'))
    if not all(r['outputs']['bcftools.fasta']['matches_known_deletion']
               for r in results if r['consistent_callset']):
        raise SystemExit('Equivalent valid representations failed known-deletion truth')
    (args.output/'COMPLETE').touch()


if __name__ == '__main__':
    main()
