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


def subtype(blast):
    for line in blast.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("query acc.ver"):
            continue
        fields = line.split("\t")
        if len(fields) >= 2 and fields[1]:
            return fields[1].split("_", 1)[0]
        raise ValueError(f"Malformed BLAST hit in {blast}")
    raise ValueError(f"No subtype hit in {blast}")


def inputs(directory):
    sample = directory.name
    sequence = directory / f"{sample}-0.15-iupac.fasta"
    blast = directory / f"{sample}-0.15-iupac.fasta.blast"
    annotations = sorted(directory.glob(f"{sample}.vadr.*.gff"))
    cram = directory / f"{sample}-0.15-iupac.cram"
    if not sequence.exists() or not blast.exists() or len(annotations) != 1:
        return None
    return sample, sequence, blast, annotations[0], cram if cram.exists() else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path)
    parser.add_argument("--rules", type=Path, required=True)
    parser.add_argument(
        "--h77-fasta",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "refgenomes/resistance/AF011751.1.fasta",
    )
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
                f"WOULD RECALL\t{item[0]}\t{directory}\tsubtype={subtype(item[2])}\t{'with CRAM evidence' if item[4] else 'baseline only'}"
            )
        print(f"Dry run: {len(candidates)} sample(s); use --apply to write outputs")
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup_root = args.backup_dir or args.archives.with_name(
        f"{args.archives.name}-resistance-backup-{stamp}"
    )
    for directory, (sample, sequence, blast, annotation, cram) in candidates:
        with tempfile.TemporaryDirectory(dir=directory) as temporary:
            command = [
                sys.executable,
                str(Path(__file__).with_name("annotate_vcf_resistance.py")),
                "--fasta",
                str(sequence),
                "--gff",
                str(annotation),
                "--subtype",
                subtype(blast),
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
