#!/usr/bin/env python3
"""Re-call resistance in existing VirPipa sample archives without rerunning assembly."""
import argparse, shutil, subprocess, sys, tempfile
from datetime import datetime, timezone
from pathlib import Path

SUFFIXES = (
    "_resistance.tsv",
    "_resistance.bed",
    "_resistance.gff",
    "_resistance_by_drug.tsv",
    "_resistance.json",
    "_resistance_sites.gff3",
)


def subtype(report):
    for line in report.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[0] == "reference":
            return fields[1].split("-", 1)[0]
    raise ValueError(f"No reference subtype in {report}")


def inputs(directory):
    sample = directory.name
    sequence = directory / f"{sample}-0.15-iupac.fasta"
    report = directory / f"{sample}-0.15-iupac.report.tsv"
    annotations = sorted(directory.glob(f"{sample}.vadr.*.gff"))
    cram = directory / f"{sample}-0.15-iupac.cram"
    if not sequence.exists() or not report.exists() or len(annotations) != 1:
        return None
    return sample, sequence, report, annotations[0], cram if cram.exists() else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path)
    parser.add_argument("--rules", type=Path, required=True)
    parser.add_argument("--h77-fasta", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--backup-dir", type=Path)
    args = parser.parse_args()
    candidates = [
        (directory, inputs(directory))
        for directory in sorted(args.archives.iterdir())
        if directory.is_dir()
    ]
    candidates = [item for item in candidates if item[1]]
    if not args.apply:
        for directory, item in candidates:
            print(
                f"WOULD RECALL\t{item[0]}\t{directory}\t{'with CRAM evidence' if item[4] else 'baseline only'}"
            )
        print(f"Dry run: {len(candidates)} sample(s); use --apply to write outputs")
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup_root = args.backup_dir or args.archives.with_name(
        f"{args.archives.name}-resistance-backup-{stamp}"
    )
    for directory, (sample, sequence, report, annotation, cram) in candidates:
        with tempfile.TemporaryDirectory(dir=directory) as temporary:
            command = [
                sys.executable,
                str(Path(__file__).with_name("annotate_vcf_resistance.py")),
                "--fasta",
                str(sequence),
                "--gff",
                str(annotation),
                "--subtype",
                subtype(report),
                "--rules",
                str(args.rules.resolve()),
                "--h77-fasta",
                str(args.h77_fasta.resolve()),
                "--sample-name",
                sample,
                "--output-dir",
                temporary,
            ]
            if cram:
                command += ["--cram", str(cram)]
            subprocess.run(command, check=True)
            generated = [Path(temporary) / f"{sample}{suffix}" for suffix in SUFFIXES]
            if not all(path.exists() for path in generated):
                raise RuntimeError(f"Incomplete resistance output for {sample}")
            backup = backup_root / sample
            backup.mkdir(parents=True, exist_ok=True)
            for suffix, source in zip(SUFFIXES, generated):
                destination = directory / f"{sample}{suffix}"
                if destination.exists():
                    shutil.copy2(destination, backup / destination.name)
                source.replace(destination)
        print(f"RECALLED\t{sample}\t{directory}")
    print(
        f"Completed {len(candidates)} sample(s); replaced files were backed up under {backup_root}"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        raise SystemExit(1)
