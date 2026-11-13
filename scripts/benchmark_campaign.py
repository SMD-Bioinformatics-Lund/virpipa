#!/usr/bin/env python3
"""Validate and freeze an explicit campaign CSV before a fresh Hopper run."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path


def fingerprint(path):
    path = Path(path)
    before = path.stat()
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'Input changed while hashing: {path}')
    return dict(path=str(path), size=after.st_size, mtime_ns=after.st_mtime_ns, sha256=sha.hexdigest())


def validate(csv_path, expected):
    with Path(csv_path).open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected:
        raise ValueError(f'Expected {expected} rows, found {len(rows)}')
    samples = [r['clarity_sample_id'] for r in rows]
    # The current pipeline keys several joins by sample_id alone.
    if any(not s for s in samples) or len(set(samples)) != len(samples):
        raise ValueError('Campaign requires unique, nonempty sample IDs')
    for row in rows:
        if not row['sequencing_run']:
            raise ValueError('Missing sequencing run')
        for key in ('read1', 'read2'):
            if not Path(row[key]).is_file():
                raise ValueError(f'Missing FASTQ: {row[key]}')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected', type=int, default=240)
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    rows = validate(args.csv, args.expected)
    paths = sorted({r[key] for r in rows for key in ('read1', 'read2')})
    if args.check_only:
        payload = dict(rows=len(rows), fastqs=len(paths), samplesheet=fingerprint(args.csv))
    else:
        with ThreadPoolExecutor(max_workers=4) as pool:
            inputs = list(pool.map(fingerprint, paths))
        payload = dict(samplesheet=fingerprint(args.csv), samples=rows, inputs=inputs)
    with args.output.open('x') as stream:
        json.dump(payload, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
