#!/usr/bin/env python3
"""Check fresh baseline publication and compare biological payloads to inventory."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import runpy
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', required=True, type=Path)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--scripts', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--hopper-paths', action='store_true')
    args = parser.parse_args()
    sys.path.insert(0, str(args.scripts))
    import pysam
    from consensus_metadata import sequence
    from benchmark_consensus import digest

    records = json.loads(args.inventory.read_text())
    expected = {(r['run_name'], r['sample_id']) for r in records}
    results = args.baseline/'results'
    published = {(p.parent.parent.name, p.name.removesuffix('_qc_summary.json'))
                 for p in results.glob('*/*/*_qc_summary.json')}
    report = dict(expected=len(expected), published=len(published),
                  missing=sorted(expected-published), extra=sorted(published-expected),
                  samples=[], failures=[], task_statuses={})
    traces = list(results.glob('pipeline_info_*/trace.txt'))
    if len(traces) != 1:
        report['failures'].append('Expected one baseline task trace')
    else:
        with traces[0].open() as stream:
            tasks = list(csv.DictReader(stream, delimiter='\t'))
        report['task_statuses'] = dict(Counter(r['status'] for r in tasks))
        report['failed_tasks'] = [r for r in tasks if r['status'] not in ('COMPLETED', 'CACHED')
                                  or r.get('exit') not in ('0', '-', '')]
        if report['failed_tasks']:
            report['failures'].append('Task trace contains unsuccessful tasks')
    layout = runpy.run_path(str(args.scripts/'check_publish_layout.py'))['main']
    fasta_names = {'scaffold': '.fasta', 'iupac': '-0.15-iupac.fasta', 'majority': '-majority.fasta'}
    for record in records:
        run, sample = record['run_name'], record['sample_id']
        directory = results/run/sample
        try:
            sys.argv = ['check_publish_layout.py', '--sample-dir', str(directory),
                        '--sample-id', sample, '--mode', 'debug']
            layout()
            differences = []
            for key, suffix in fasta_names.items():
                old = record['paths'][key]
                if args.hopper_paths and old.startswith(('/mnt/fs1/', '/mnt/fs2/')):
                    old = old[4:]
                _, previous = sequence(old)
                identifier, current = sequence(directory/(sample+suffix))
                if current != previous:
                    differences.append(key)
                if key in ('scaffold', 'iupac'):
                    cram_name = sample+('' if key == 'scaffold' else '-0.15-iupac')+'.cram'
                    with pysam.AlignmentFile(str(directory/cram_name), 'rc',
                            reference_filename=str(directory/(sample+suffix))) as handle:
                        if handle.references != (identifier,) or handle.lengths != (len(current),):
                            raise ValueError('CRAM reference identity/length mismatch')
                        if handle.header.to_dict()['SQ'][0].get('M5') != hashlib.md5(current.encode()).hexdigest():
                            raise ValueError('CRAM reference MD5 mismatch')
                        handle.check_index()
                        next(handle.fetch(identifier, 0, min(120, len(current))), None)
            old_resistance = record['paths']['resistance']
            if args.hopper_paths and old_resistance.startswith('/mnt/fs1/'):
                old_resistance = old_resistance[4:]
            before = json.loads(Path(old_resistance).read_text())
            after = json.loads((directory/(sample+'_resistance.json')).read_text())
            def outcomes(data):
                return sorted((r['drug'], r['outcome'], r['prediction']) for r in data['drug_outcomes'])
            report['samples'].append(dict(run_name=run, sample_id=sample, sequence_changes=differences,
                subtype_changed=before['subtype'] != after['subtype'],
                drug_outcomes_changed=outcomes(before) != outcomes(after)))
        except (Exception, SystemExit) as exc:
            report['failures'].append(dict(run_name=run, sample_id=sample, error=str(exc)))
    report['summary'] = {key: sum(bool(r[key]) for r in report['samples'])
                         for key in ('sequence_changes', 'subtype_changed', 'drug_outcomes_changed')}
    report['inventory_sha256'] = digest(args.inventory)
    report['baseline_commit'] = (args.baseline/'source_commit.txt').read_text().strip()
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('samples', 'failed_tasks')}, indent=2))
    if report['failures'] or report['missing'] or report['extra']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
