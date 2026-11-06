#!/usr/bin/env python3
"""Compare paired FASTQs by ordered content and an order-independent pair multiset."""

import argparse
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile


def records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            record = [header] + [handle.readline() for _ in range(3)]
            record = [line.rstrip(b"\r\n") for line in record]
            if (not record[0].startswith(b"@") or not record[2].startswith(b"+")
                    or not record[1] or len(record[1]) != len(record[3])):
                raise ValueError(f"Malformed four-line FASTQ record in {path}")
            yield record


def pair_id(header):
    name = header.split()[0][1:]
    return name[:-2] if name.endswith((b"/1", b"/2")) else name


def fingerprint(read1, read2, directory, label):
    digests = directory / f"{label}.digests"
    ordered, identifiers = hashlib.sha256(), hashlib.sha256()
    count = bases = 0
    with digests.open("wb") as output:
        for left, right in itertools.zip_longest(records(read1), records(read2)):
            if left is None or right is None:
                raise ValueError("R1/R2 have different record counts")
            name = pair_id(left[0])
            if name != pair_id(right[0]):
                raise ValueError(f"R1/R2 read identifiers differ at pair {count + 1}")
            # Keep complete headers, sequences, plus lines and qualities. Framing
            # by normalized newlines cannot conflate fields of a FASTQ record.
            digest = hashlib.sha256(b"\n".join(left + right) + b"\n").digest()
            output.write(digest.hex().encode() + b"\n")
            ordered.update(digest)
            identifiers.update(name + b"\n")
            count += 1
            bases += len(left[1]) + len(right[1])
    sorted_path = directory / f"{label}.sorted"
    with sorted_path.open("wb") as output:
        subprocess.run(["sort", "-S", "32M", "--parallel=1", "-T", str(directory), str(digests)],
                       stdout=output, check=True, env={**os.environ, "LC_ALL": "C"})
    return dict(read1=str(read1), read2=str(read2), pairs=count, bases=bases,
                ordered_pairs_sha256=ordered.hexdigest(),
                ordered_pair_ids_sha256=identifiers.hexdigest()), sorted_path


def grouped(handle):
    for key, group in itertools.groupby(handle):
        yield key, sum(1 for _ in group)


def multiset_difference(left, right):
    """Merge sorted hashes, preserving the multiplicity of every paired record."""
    removed = added = 0
    with left.open("rb") as a, right.open("rb") as b:
        old, new = grouped(a), grouped(b)
        x, y = next(old, None), next(new, None)
        while x is not None or y is not None:
            if y is None or (x is not None and x[0] < y[0]):
                removed += x[1]
                x = next(old, None)
            elif x is None or y[0] < x[0]:
                added += y[1]
                y = next(new, None)
            else:
                removed += max(0, x[1] - y[1])
                added += max(0, y[1] - x[1])
                x, y = next(old, None), next(new, None)
    return removed, added


def compare(left_r1, left_r2, right_r1, right_r2, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix="pair-hashes-", dir=output_dir) as temporary:
        directory = Path(temporary)
        left, old = fingerprint(left_r1, left_r2, directory, "left")
        right, new = fingerprint(right_r1, right_r2, directory, "right")
        removed, added = multiset_difference(old, new)
    result = dict(left=left, right=right,
                  same_ordered_pairs=left["ordered_pairs_sha256"] == right["ordered_pairs_sha256"],
                  same_pair_id_order=left["ordered_pair_ids_sha256"] == right["ordered_pair_ids_sha256"],
                  same_pair_multiset=not (removed or added),
                  pairs_only_in_left=removed, pairs_only_in_right=added)
    (output_dir / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("left-r1", "left-r2", "right-r1", "right-r2", "outdir"):
        parser.add_argument("--" + option, type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.left_r1, args.left_r2, args.right_r1, args.right_r2, args.outdir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
