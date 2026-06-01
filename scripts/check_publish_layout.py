#!/usr/bin/env python3
"""Validate the fixture-backed published output layout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check VirPipa published output layout")
    parser.add_argument("--sample-dir", required=True, help="Published <outdir>/<run>/<sample> directory")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--mode", choices=["routine", "debug"], required=True, help="Expected publish mode")
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def main() -> None:
    args = parse_args()
    sample_dir = Path(args.sample_dir)
    sample_id = args.sample_id
    qc_path = sample_dir / f"{sample_id}_qc_summary.json"

    require(sample_dir.is_dir(), f"Missing sample directory: {sample_dir}")
    require(qc_path.is_file(), f"Missing QC summary: {qc_path}")
    require(not (sample_dir / "results").exists(), f"Unexpected legacy results directory: {sample_dir / 'results'}")
    require(not (sample_dir / "lid").exists(), f"Unexpected LID directory: {sample_dir / 'lid'}")

    data = json.loads(qc_path.read_text(encoding="utf-8"))
    outputs = data.get("outputs") or {}
    require("lid_2limsrs" not in outputs, "QC outputs still includes lid_2limsrs")
    require(outputs.get("display_rug_kde_plot") == f"{sample_id}_display_rug_kde_plot.png", "Unexpected display plot path")

    for key, rel_path in outputs.items():
        if rel_path is None:
            continue
        output_path = sample_dir / rel_path
        require(not Path(rel_path).is_absolute(), f"{key} is absolute: {rel_path}")
        require(output_path.exists(), f"{key} does not resolve: {output_path}")
        require(sample_dir.resolve() in [output_path.resolve(), *output_path.resolve().parents], f"{key} escapes sample dir: {rel_path}")

    for rel_path in [
        f"{sample_id}.fasta.fai",
        f"{sample_id}-0.15-iupac.fasta.fai",
        f"{sample_id}.cram.crai",
        f"{sample_id}-0.15-iupac.cram.crai",
        f"{sample_id}-pilon-m0.15.vcf.gz.csi",
    ]:
        require((sample_dir / rel_path).is_file(), f"Missing required sidecar: {rel_path}")

    if args.mode == "routine":
        for name in ["bam", "vcf", "fastq", "fasta", "spades", "mummer", "pilon"]:
            require(not (sample_dir / name).exists(), f"Unexpected debug directory in routine mode: {name}")
    else:
        require(any((sample_dir / name).exists() for name in ["bam", "vcf", "fastq", "fasta", "spades", "mummer", "pilon"]), "No debug directories found")


if __name__ == "__main__":
    main()
