#!/usr/bin/env python3
"""Preview or transactionally re-call resistance in completed sample archives."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from build_qc_summary import read_resistance
from consensus_metadata import load as load_consensus

SUFFIXES = ("_resistance.tsv", "_resistance.bed", "_resistance.gff",
            "_resistance_by_drug.tsv", "_resistance.json", "_resistance_sites.gff3")
OUTPUT_KEYS = ("resistance_tsv", "resistance_bed", "resistance_gff",
               "resistance_by_drug_tsv", "resistance_json", "resistance_sites_gff3")


def subtype(blast):
    for line in blast.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("query acc.ver"):
            continue
        fields = line.split("\t")
        if len(fields) >= 2 and fields[1]:
            return fields[1].split("_", 1)[0]
        raise ValueError(f"Malformed BLAST hit in {blast}")
    raise ValueError(f"No subtype hit in {blast}")


def inputs(directory, sample=None):
    sample = sample or directory.name
    sequence = directory / f"{sample}-0.15-iupac.fasta"
    blast = directory / f"{sample}-0.15-iupac.fasta.blast"
    annotations = sorted(directory.glob(f"{sample}.vadr.*_mod.gff"))
    if not annotations:
        annotations = sorted(directory.glob(f"{sample}.vadr.*.gff"))
    qc = directory / f"{sample}_qc_summary.json"
    if not sequence.exists() or not blast.exists() or len(annotations) != 1 or not qc.exists():
        raise ValueError("requires consensus, BLAST, QC JSON, and one normalized VADR GFF")
    metadata = load_consensus(sequence)
    if metadata["mode"] != "iupac" or metadata["min_af"] != 0.15 or metadata["minimum_usable_depth"] != 1 or metadata["reference_fill"]:
        raise ValueError("legacy/unverified consensus: regenerate consensus or rerun pipeline first")
    scaffold = directory / f"{sample}.fasta"
    if scaffold.exists():
        scaffold_metadata = load_consensus(scaffold)
        if metadata["coordinate_source_sha256"] != scaffold_metadata["coordinate_source_sha256"] or scaffold_metadata["mode"] != "scaffold":
            raise ValueError("main scaffold and analytical consensus lack a shared coordinate origin")
    cram = directory / f"{sample}-0.15-iupac.cram"
    return sample, sequence, blast, annotations[0], cram if cram.exists() else None


def dump(path, payload):
    fd, name = tempfile.mkstemp(prefix=".resistance-json-", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def replace_file(source, destination):
    """Copy then rename on the destination filesystem; never truncate in place."""
    fd, name = tempfile.mkstemp(prefix=".resistance-recall-", dir=destination.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        shutil.copy2(source, temporary)
        if destination.exists():
            temporary.chmod(destination.stat().st_mode & 0o777)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def release_lock(archives):
    lock = archives / ".resistance-recall-lock"
    if lock.exists():
        (lock / "backup.txt").unlink(missing_ok=True)
        lock.rmdir()


def restore(backup_root, archives):
    manifest_path = backup_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if Path(manifest["archives"]).resolve() != archives.resolve():
        raise ValueError("Backup belongs to a different archive root")
    lock_path = archives / ".resistance-recall-lock" / "backup.txt"
    if lock_path.exists() and Path(lock_path.read_text().strip()).resolve() != backup_root.resolve():
        raise ValueError("A different recall transaction owns the archive lock")
    if manifest["state"] == "restored":
        release_lock(archives)
        return
    if manifest["state"] not in {"prepared", "applying", "complete"}:
        raise ValueError("Backup was not completely prepared")
    targets = []
    for item in manifest["files"]:
        relative = Path(item["path"])
        destination = archives / relative
        source = backup_root / "files" / relative
        if relative.is_absolute() or ".." in relative.parts or not destination.resolve().is_relative_to(archives.resolve()) or not source.resolve().is_relative_to(backup_root.resolve()):
            raise ValueError("Unsafe path in backup manifest")
        if item["existed"] and not source.is_file():
            raise ValueError(f"Incomplete backup: {relative}")
        targets.append((item, source, destination))
    for item, source, destination in targets:
        if item["existed"]:
            replace_file(source, destination)
        else:
            destination.unlink(missing_ok=True)
    manifest["state"] = "restored"
    dump(manifest_path, manifest)
    release_lock(archives)


def apply_transaction(archives, replacements, backup_root):
    """Back up the complete set before writing; restore it on any caught failure."""
    lock = archives / ".resistance-recall-lock"
    lock.mkdir()  # Existing lock requires explicit recovery, not a second writer.
    manifest = {"archives": str(archives.resolve()), "state": "preparing", "files": []}
    prepared = False
    try:
        backup_root.mkdir(parents=True, exist_ok=False)
        (lock / "backup.txt").write_text(str(backup_root.resolve())+"\n")
        for destination in replacements:
            relative = destination.relative_to(archives)
            if not destination.resolve().is_relative_to(archives.resolve()):
                raise ValueError(f"Archive target escapes root: {destination}")
            if destination.is_symlink():
                raise ValueError(f"Refusing to replace symlinked archive output: {destination}")
            item = {"path": str(relative), "existed": destination.exists()}
            manifest["files"].append(item)
            if item["existed"]:
                target = backup_root / "files" / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(destination, target)
        manifest["state"] = "prepared"
        dump(backup_root / "manifest.json", manifest)
        prepared = True
        manifest["state"] = "applying"
        dump(backup_root / "manifest.json", manifest)
        for destination, source in replacements.items():
            replace_file(source, destination)
        manifest["state"] = "complete"
        dump(backup_root / "manifest.json", manifest)
    except BaseException:
        if prepared:
            restore(backup_root, archives)
        raise
    finally:
        if lock.exists() and (not prepared or manifest["state"] == "complete" or json.loads((backup_root / "manifest.json").read_text())["state"] == "restored"):
            release_lock(archives)


def stage_sample(directory, item, args, staged):
    sample, sequence, blast, annotation, cram = item
    staged.mkdir()
    command = [sys.executable, str(Path(__file__).with_name("annotate_vcf_resistance.py")),
               "--fasta", str(sequence), "--gff", str(annotation), "--subtype", subtype(blast),
               "--rules", str(args.rules.resolve()), "--h77-fasta", str(args.h77_fasta.resolve()),
               "--sample-name", sample, "--output-dir", str(staged)]
    if cram:
        command += ["--cram", str(cram)]
    main_fasta = directory / f"{sample}.fasta"
    if main_fasta.exists():
        command += ["--track-reference-fasta", str(main_fasta)]
    subprocess.run(command, check=True)
    generated = [staged / f"{sample}{suffix}" for suffix in SUFFIXES]
    if not all(path.is_file() for path in generated):
        raise RuntimeError(f"Incomplete resistance output for {sample}")
    evidence = json.loads(generated[4].read_text())
    if evidence["schema_version"] != "3.1" or not evidence["drug_outcomes"]:
        raise ValueError(f"Invalid resistance JSON for {sample}")
    qc_path = directory / f"{sample}_qc_summary.json"
    qc = json.loads(qc_path.read_text())
    if qc.get("sample_id") != sample:
        raise ValueError(f"QC sample identity mismatch for {sample}")
    qc["resistance"] = read_resistance(generated[0], generated[3], generated[4])
    for key, path in zip(OUTPUT_KEYS, generated):
        qc.setdefault("outputs", {})[key] = path.name
    qc["resistance_recalled_at_utc"] = datetime.now(timezone.utc).isoformat()
    dump(staged / qc_path.name, qc)
    return generated, qc_path, qc


def stage_aggregates(archives, updates, staged):
    """Refresh existing aggregate formats, preserving records for other samples."""
    current = {}
    for path in archives.glob("**/*_qc_summary.json"):
        record = updates.get(path) or json.loads(path.read_text())
        current[(record.get("run_name"), record.get("sample_id"))] = record
    current.update(((record.get("run_name"), record.get("sample_id")), record) for record in updates.values())
    replacements = {}
    for path in sorted(archives.glob("**/qc_summary.json*")):
        if path.name not in {"qc_summary.json", "qc_summary.jsonl"}:
            continue
        records = ([json.loads(line) for line in path.read_text().splitlines() if line.strip()]
                   if path.suffix == ".jsonl" else json.loads(path.read_text()))
        if not isinstance(records, list):
            raise ValueError(f"Unsupported aggregate layout: {path}")
        records = [current.get((r.get("run_name"), r.get("sample_id")), r) for r in records]
        target = staged / f"aggregate-{len(replacements)}{path.suffix}"
        if path.suffix == ".jsonl":
            target.write_text("".join(json.dumps(r, sort_keys=True)+"\n" for r in records))
        else:
            dump(target, records)
        replacements[path] = target
    return replacements


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path)
    parser.add_argument("--rules", type=Path)
    parser.add_argument("--h77-fasta", type=Path, default=Path(__file__).resolve().parents[1] / "refgenomes/resistance/AF011751.1.fasta")
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--backup-dir", type=Path)
    parser.add_argument("--restore-backup", type=Path)
    args = parser.parse_args()
    args.archives = args.archives.resolve()
    if args.restore_backup:
        restore(args.restore_backup.resolve(), args.archives)
        print("Restored resistance files and QC summaries from backup")
        return
    if not args.rules:
        parser.error("--rules is required unless restoring a backup")
    lock = args.archives / ".resistance-recall-lock"
    if lock.exists():
        raise RuntimeError(f"Interrupted/concurrent recall: inspect {lock}/backup.txt and use --restore-backup")
    candidates = []
    for child in sorted(args.archives.iterdir()):
        if not child.is_dir() or child.name.startswith(".") or child.name == "pipeline_info":
            continue
        directory = child if (child / f"{child.name}_qc_summary.json").exists() else (child / "results" if (child / "results").is_dir() else child)
        try:
            item = inputs(directory, child.name)
            value = subtype(item[2])
        except ValueError as exc:
            print(f"SKIP\t{child.name}\t{exc}")
            continue
        candidates.append((directory, item))
        print(f"WOULD RECALL\t{item[0]}\t{directory}\tsubtype={value}\t{'with CRAM evidence' if item[4] else 'baseline only'}")
    if not args.apply or not candidates:
        print(f"{'Dry run' if not args.apply else 'No eligible inputs'}: {len(candidates)} sample(s)")
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup_root = (args.backup_dir or args.archives.with_name(f"{args.archives.name}-resistance-backup-{stamp}")).resolve()
    if backup_root.is_relative_to(args.archives):
        raise ValueError("Backup directory must be outside the archive root")
    with tempfile.TemporaryDirectory(prefix="virpipa-recall-") as temporary:
        staged_root = Path(temporary)
        replacements, qc_replacements, updates = {}, {}, {}
        for index, (directory, item) in enumerate(candidates):
            staged = staged_root / str(index)
            generated, qc_path, qc = stage_sample(directory, item, args, staged)
            replacements.update((directory / path.name, path) for path in generated)
            qc_replacements[qc_path] = staged / qc_path.name
            updates[qc_path] = qc
        replacements.update(qc_replacements)  # Sample completion markers last.
        replacements.update(stage_aggregates(args.archives, updates, staged_root))
        apply_transaction(args.archives, replacements, backup_root)
    print(f"Completed {len(candidates)} sample(s); backup: {backup_root}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        raise SystemExit(1)
