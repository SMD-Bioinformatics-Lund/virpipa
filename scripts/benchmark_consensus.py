#!/usr/bin/env python3
"""Isolated consensus capability benchmark; never writes published archives.

Tool configuration is JSON mapping bcftools, samtools and ivar to argv lists.
All commands, versions and input hashes are retained in the new output directory.
"""
import argparse
import copy
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import subprocess

import pysam

from annotate_vcf_resistance import IUPAC, TABLE, read_codon
from consensus_metadata import generate, sequence

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    sha = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(chunk)
    return sha.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


class Runner:
    def __init__(self, tools, output):
        self.tools, self.output = tools, Path(output)
        self.commands = []

    def run(self, tool, args, *, stdout=None, stdin=None):
        command = self.tools[tool] + list(map(str, args))
        with (self.output / 'tools.stderr.log').open('a') as log:
            result = subprocess.run(command, input=stdin, stdout=subprocess.PIPE,
                                    stderr=log)
        entry = dict(argv=command, returncode=result.returncode,
                     stdout=str(stdout) if stdout else None)
        self.commands.append(entry)
        write_json(self.output / 'commands.json', self.commands)
        if result.returncode:
            raise RuntimeError(f'{tool} failed; see {self.output}/tools.stderr.log')
        if stdout:
            Path(stdout).write_bytes(result.stdout)
        return result.stdout


def codon_counts(bam, reference, positions, min_bq=13):
    """Use production codon extraction; one count per agreeing template."""
    templates = {}
    with pysam.AlignmentFile(str(bam), reference_filename=str(reference)) as reads:
        for read in reads.fetch(reads.references[0], min(positions), max(positions)+1):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_qcfail or read.is_duplicate):
                continue
            codon = read_codon(read, positions, min_bq)
            if codon:
                templates.setdefault(read.query_name, set()).add(codon)
    codons = Counter(next(iter(values)) for values in templates.values() if len(values) == 1)
    amino = Counter()
    for codon, count in codons.items():
        amino[TABLE[codon]] += count
    total = sum(amino.values())
    return dict(codons=dict(codons), amino_acids=dict(amino), depth=total,
                qualifying_amino_acids=sorted(aa for aa, n in amino.items()
                                             if n * 100 >= 15 * total))


def fixtures(output):
    """Known molecules and explicit CIGARs, not truth inferred by a caller."""
    output = Path(output)
    # A deterministic HCV reference segment supplies realistic flanking context.
    _, h77 = sequence(ROOT / 'refgenomes/resistance/AF011751.1.fasta')
    ref = h77[3400:3640]
    ref = ref[:60] + 'TTT' + ref[63:90] + 'GCCCCCCCCCC' + ref[101:]
    cases = []
    for fraction in (14, 15, 16):
        cases.append((f'snv_{fraction}', [('TTT', 100-fraction), ('TCT', fraction)], None))
    cases += [('synonymous_sum', [('TTT', 84), ('GCT', 8), ('GCC', 8)], None),
              ('unobserved_codon', [('ATA', 50), ('TTT', 50)], None),
              ('multiallelic', [('TTT', 60), ('TCT', 20), ('TAT', 20)], None),
              ('low_depth', [('TCT', 2)], None),
              ('poor_quality', [('TCT', 20)], None),
              ('agreeing_mates', [('TCT', 20)], None),
              ('conflicting_mates', [('TCT', 20)], None),
              ('partial_codon', [('TCT', 20)], None),
              ('no_coverage', [], None)]
    for length in (1, 2, 3):
        for fraction in (8, 15, 50, 90):
            for kind in ('insertion', 'deletion'):
                cases.append((f'{kind}_{length}_{fraction}', [('TTT', 100)],
                              (kind, length, fraction)))
    manifest = []
    for name, molecules, indel in cases:
        dest = output / name
        dest.mkdir(parents=True)
        fasta = dest / 'reference.fasta'
        fasta.write_text('>truth\n'+ref+'\n')
        pysam.faidx(str(fasta))
        bam = dest / 'reads.bam'
        expected = Counter()
        with pysam.AlignmentFile(str(bam), 'wb', header={
                'HD': {'SO': 'coordinate'}, 'SQ': [{'SN': 'truth', 'LN': len(ref)}]}) as handle:
            index = 0
            for codon, count in molecules:
                for _ in range(count):
                    seq = ref[:60]+codon+ref[63:]
                    cigar = [(0, len(ref))]
                    if indel and index < indel[2]:
                        kind, length, _ = indel
                        if kind == 'insertion':
                            seq = seq[:91]+'C'*length+seq[91:]
                            cigar = [(0, 91), (1, length), (0, len(ref)-91)]
                        else:
                            seq = seq[:91]+seq[91+length:]
                            cigar = [(0, 91), (2, length), (0, len(ref)-91-length)]
                    read = pysam.AlignedSegment()
                    read.query_name = f'molecule_{index}'
                    read.query_sequence, read.flag = seq, 0
                    read.reference_id, read.reference_start, read.mapping_quality = 0, 0, 60
                    read.cigartuples = cigar
                    read.query_qualities = pysam.qualitystring_to_array('I'*len(seq))
                    if name == 'poor_quality':
                        qualities = read.query_qualities
                        qualities[61] = 2
                        read.query_qualities = qualities
                    if name == 'partial_codon':
                        read.query_sequence = seq[:62]
                        read.cigartuples = [(0, 62)]
                        read.query_qualities = pysam.qualitystring_to_array('I'*62)
                    handle.write(read)
                    if name in ('agreeing_mates', 'conflicting_mates'):
                        mate = copy.copy(read)
                        if name == 'conflicting_mates':
                            mate.query_sequence = seq[:60]+'TTT'+seq[63:]
                            mate.query_qualities = pysam.qualitystring_to_array('I'*len(seq))
                        handle.write(mate)
                    if name not in ('poor_quality', 'partial_codon', 'conflicting_mates'):
                        expected[codon] += 1
                    index += 1
        pysam.index(str(bam))
        aa = Counter()
        for codon, count in expected.items():
            aa[TABLE[codon]] += count
        total = sum(expected.values())
        truth = dict(case=name, reference=str(fasta), bam=str(bam),
                     codon_positions_zero_based=[60, 61, 62], codons=dict(expected),
                     amino_acids=dict(aa), depth=total, indel=indel,
                     qualifying_amino_acids=sorted(a for a, n in aa.items()
                                                  if n*100 >= 15*total))
        write_json(dest/'truth.json', truth)
        manifest.append(truth)
    write_json(output/'fixtures.json', manifest)
    return manifest


def run_case(case, runner, output):
    dest = output / case['case']
    dest.mkdir()
    ref, bam = Path(case['reference']), Path(case['bam'])
    bcf, vcf, normalized = dest/'calls.bcf', dest/'calls.vcf.gz', dest/'normalized.vcf.gz'
    # BAQ is disabled in this controlled probe, consistently for both engines.
    runner.run('bcftools', ['mpileup', '-B', '-Q', '13', '-d', '1000000',
                           '-a', 'FORMAT/AD,FORMAT/DP', '-f', ref, '-Ob', '-o', bcf, bam])
    runner.run('bcftools', ['call', '-m', '-A', '--ploidy', '1', '-Oz', '-o', vcf, bcf])
    runner.run('bcftools', ['view', vcf], stdout=dest/'calls.vcf')
    runner.run('bcftools', ['norm', '-f', ref, '-Oz', '-o', normalized, vcf])
    runner.run('bcftools', ['index', normalized])
    runner.run('bcftools', ['view', normalized], stdout=dest/'normalized.vcf')
    runner.run('bcftools', ['consensus', '-f', ref, '-s', '-', '-H', 'I', normalized],
               stdout=dest/'bcftools-alt-diagnostic.fasta')
    # The sample's called haploid genotype cannot express minority mixtures.
    # Using its actual ID, rather than -s -, makes that limitation measurable.
    sample = runner.run('bcftools', ['query', '-l', normalized]).decode().strip()
    runner.run('bcftools', ['consensus', '-f', ref, '-s', sample, '-H', 'I', normalized],
               stdout=dest/'bcftools-genotype.fasta')
    generate(dest/'calls.vcf', ref, dest/'baseline.fasta', 'truth', 'iupac')
    generate(dest/'normalized.vcf', ref, dest/'baseline-normalized-diagnostic.fasta', 'truth', 'iupac')
    if 'ivar' in runner.tools:
        pileup = runner.run('samtools', ['mpileup', '-aa', '-A', '-B', '-Q', '0',
                                         '-d', '1000000', '-f', ref, bam])
        (dest/'reads.pileup').write_bytes(pileup)
        for threshold in ('0', '0.15', '0.85', '0.90'):
            runner.run('ivar', ['consensus', '-q', '13', '-m', '1', '-n', 'N',
                               '-c', '0.8', '-t', threshold,
                               '-p', dest/f'ivar-{threshold}'], stdin=pileup)
    actual = codon_counts(bam, ref, case['codon_positions_zero_based'])
    matches = all(actual[key] == case[key] for key in
                  ('codons', 'amino_acids', 'depth', 'qualifying_amino_acids'))
    sequences = {}
    for fasta in sorted(list(dest.glob('*.fasta'))+list(dest.glob('*.fa'))):
        _, seq = sequence(fasta)
        codon = seq[60:63]
        possibilities = (sorted({TABLE[''.join(bases)] for bases in
                                 itertools.product(*(IUPAC[b] for b in codon))})
                         if len(codon) == 3 and 'N' not in codon else [])
        sequences[fasta.name] = dict(length=len(seq), unknown=seq.count('N'),
                                    codon_61_63=codon, amino_acid_possibilities=possibilities,
                                    extra_amino_acids=sorted(set(possibilities)-set(case['qualifying_amino_acids'])),
                                    missing_amino_acids=sorted(set(case['qualifying_amino_acids'])-set(possibilities)),
                                    sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
    return dict(case=case['case'], truth_pass=matches, codon_evidence=actual,
                outputs=sequences, ivar_tested='ivar' in runner.tools)


def prepare_probe_distances(output):
    """Create coverage-backed distance inputs for coordinate-compatible probes.

    A single simulated indel changes length if applied; such outputs are excluded
    explicitly until a tool-provided coordinate map is available.
    """
    output = Path(output)
    manifest, excluded = [], []
    destination = output/'distance-inputs'
    destination.mkdir()
    for case in json.loads((output/'fixtures/fixtures.json').read_text()):
        root = output/'results'/case['case']
        generate(root/'calls.vcf', case['reference'], root/'baseline-majority.fasta',
                 'truth', 'analytical_majority')
        _, ref = sequence(case['reference'])
        with pysam.AlignmentFile(case['bam']) as reads:
            counts = reads.count_coverage(reads.references[0], quality_threshold=13)
        mask = [sum(values) > 0 for values in zip(*counts)]
        mask_path = destination/(case['case']+'.mask.json')
        write_json(mask_path, mask)
        methods = {'baseline': ('baseline.fasta', 'baseline-majority.fasta'),
                   'bcftools_genotype': ('bcftools-genotype.fasta', 'bcftools-genotype.fasta')}
        for cutoff in ('0.15', '0.85', '0.90'):
            if (root/f'ivar-{cutoff}.fa').exists():
                methods['ivar_'+cutoff.replace('.', '')] = (f'ivar-{cutoff}.fa', 'ivar-0.fa')
        for method, variants in methods.items():
            for mode, name in zip(('iupac', 'majority'), variants):
                path = root/name
                if len(sequence(path)[1]) != len(ref):
                    excluded.append(dict(case=case['case'], method=method, mode=mode,
                                         reason='Length-changing output needs explicit coverage projection'))
                    continue
                manifest.append(dict(sample_run_id=case['case'], method=method, mode=mode,
                                     subtype='synthetic', fasta=str(path.resolve()), mask=str(mask_path.resolve())))
    write_json(destination/'manifest.json', manifest)
    write_json(destination/'excluded.json', excluded)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--tools', type=Path, help='JSON command argv arrays; omit to generate fixtures only')
    parser.add_argument('--prepare-distances', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    cases = fixtures(args.output/'fixtures')
    if not args.tools:
        return
    tools = json.loads(args.tools.read_text())
    if any(not isinstance(v, list) or not v or any(not isinstance(x, str) for x in v)
           for v in tools.values()):
        raise ValueError('Tool commands must be nonempty string arrays')
    runner = Runner(tools, args.output)
    versions = {name: runner.run(name, ['version' if name == 'ivar' else '--version']).decode()
                for name in tools}
    write_json(args.output/'provenance.json', dict(versions=versions, tools=tools,
        commit=subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD']).decode().strip(),
        scripts={str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), ROOT/'scripts/vcf_to_iupac.awk',
                    ROOT/'scripts/consensus_metadata.py', ROOT/'scripts/annotate_vcf_resistance.py']},
        inputs={c['case']: {k: digest(c[k]) for k in ('reference', 'bam')} for c in cases}))
    output = args.output/'results'
    output.mkdir()
    results = [run_case(case, runner, output) for case in cases]
    write_json(args.output/'results.json', results)
    if args.prepare_distances:
        prepare_probe_distances(args.output)
    if not all(r['truth_pass'] for r in results):
        raise SystemExit('Codon evidence disagrees with known molecules')


if __name__ == '__main__':
    main()
