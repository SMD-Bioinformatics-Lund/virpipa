#!/usr/bin/env python3
"""Create de-identified resistance fixtures and keep the re-identification key private."""
import argparse
import json
import os
import secrets
from pathlib import Path


def within(path, root):
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def rewrite_fasta(source, destination, anonymous_id):
    lines = source.read_text(encoding="utf-8").splitlines()
    destination.write_text(
        ">"
        + anonymous_id
        + "\n"
        + "\n".join(line for line in lines if not line.startswith(">"))
        + "\n",
        encoding="utf-8",
    )


def rewrite_gff(source, destination, anonymous_id):
    lines = []
    for raw in source.read_text(encoding="utf-8").splitlines():
        if raw.startswith("#"):
            lines.append(raw)
            continue
        fields = raw.split("\t")
        if len(fields) == 9:
            fields[0] = anonymous_id
        lines.append("\t".join(fields))
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--key", type=Path, required=True)
    parser.add_argument("--expectations", type=Path, required=True)
    parser.add_argument("--repo", type=Path, action="append", default=[])
    args = parser.parse_args()
    if any(within(args.key, repo) for repo in args.repo):
        raise SystemExit(
            "Refusing to write the re-identification key inside a repository"
        )

    sources = sorted(path for path in args.source.iterdir() if path.is_dir())
    aliases = [f"G2P_CASE_{index:03d}" for index in range(1, len(sources) + 1)]
    if args.key.exists():
        saved = {}
        for line in args.key.read_text(encoding="utf-8").splitlines()[1:]:
            alias, source = line.split("\t")
            saved[alias] = source
        by_name = {path.name: path for path in sources}
        if set(saved) != set(aliases) or set(saved.values()) != set(by_name):
            raise SystemExit("Existing key does not match the supplied source cases")
        sources = [by_name[saved[alias]] for alias in aliases]
    else:
        secrets.SystemRandom().shuffle(sources)
    expected = json.loads(args.expectations.read_text(encoding="utf-8"))
    args.output.mkdir(parents=True, exist_ok=True)
    mapping = []
    anonymous_expected = {}
    for alias, directory in zip(aliases, sources):
        source_id = directory.name
        fasta_candidates = list(directory.glob(f"{source_id}-0.15-iupac.fasta"))
        gff_candidates = list(directory.glob("*.vadr.*.gff"))
        if len(fasta_candidates) != 1 or len(gff_candidates) != 1:
            raise SystemExit(
                f"Expected one IUPAC FASTA and one VADR GFF for {source_id}"
            )
        rewrite_fasta(fasta_candidates[0], args.output / f"{alias}.fasta", alias)
        rewrite_gff(gff_candidates[0], args.output / f"{alias}.gff3", alias)
        mapping.append((alias, source_id))
        anonymous_expected[alias] = expected[source_id]
    (args.output / "expected.json").write_text(
        json.dumps(anonymous_expected, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    args.key.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.chmod(args.key.parent, 0o700)
    args.key.write_text(
        "anonymous_id\tsource_id\n"
        + "".join(f"{alias}\t{source}\n" for alias, source in mapping),
        encoding="utf-8",
    )
    os.chmod(args.key, 0o600)


if __name__ == "__main__":
    main()
