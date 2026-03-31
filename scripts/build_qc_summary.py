#!/usr/bin/env python3
"""Build a machine-readable QC summary from a sample results directory."""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path


AF_KEYS = ["0.05", "0.1", "0.15", "0.2", "0.3", "0.4"]
FILTERED_STATS_PATTERN = re.compile(r"-m(0\.\d+)\.vcf\.gz\.stats$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a QC summary JSON for one sample")
    parser.add_argument("--results-dir", required=True, help="Sample results directory")
    parser.add_argument("--run-name", required=True, help="Pipeline run name")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--lid", default="", help="Optional LID/sample_name")
    parser.add_argument("--sample-info", default="", help="Optional clarity_sample_info.json path")
    parser.add_argument("--output", required=True, help="Output JSON path")
    return parser.parse_args()


def read_first_blast_hit(blast_path: Path) -> dict:
    if not blast_path.exists():
        return {
            "main_blast_reference": None,
            "main_blast_genotype": None,
            "main_blast_identity": None,
        }

    with blast_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("query acc.ver"):
                continue
            fields = line.split("\t")
            subject = fields[1] if len(fields) > 1 else None
            return {
                "main_blast_reference": subject,
                "main_blast_genotype": subject.split("_")[0] if subject else None,
                "main_blast_identity": float(fields[2]) if len(fields) > 2 else None,
            }

    return {
        "main_blast_reference": None,
        "main_blast_genotype": None,
        "main_blast_identity": None,
    }


def read_report(report_path: Path) -> dict:
    report = {}
    if not report_path.exists():
        return report

    with report_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                report[parts[0]] = parts[1]
    return report


def maybe_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def maybe_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    return int(float(value))


def read_coverage_tsv(coverage_path: Path) -> dict:
    if not coverage_path.exists():
        return {}

    with coverage_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)
        if not row:
            return {}
        return {
            "1x": maybe_float(row.get("1x")),
            "10x": maybe_float(row.get("10x")),
            "100x": maybe_float(row.get("100x")),
            "1000x": maybe_float(row.get("1000x")),
        }


def read_af_counts(results_dir: Path, sample_id: str) -> dict:
    af_counts = {key: None for key in AF_KEYS}
    for stats_path in sorted(results_dir.glob(f"{sample_id}-pilon-m*.vcf.gz.stats")):
        match = FILTERED_STATS_PATTERN.search(stats_path.name)
        if not match:
            continue
        min_frac = match.group(1)
        if min_frac not in af_counts:
            continue
        with stats_path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("AF\t0\t0.000000\t"):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 4:
                        af_counts[min_frac] = maybe_int(fields[3])
                    break
    return af_counts


def read_hostile(hostile_path: Path) -> dict:
    if not hostile_path.exists():
        return {
            "enabled": False,
            "reads_in": None,
            "reads_removed_proportion": None,
        }

    data = json.loads(hostile_path.read_text(encoding="utf-8"))
    return {
        "enabled": True,
        "reads_in": data.get("reads_in"),
        "reads_removed_proportion": data.get("reads_removed_proportion"),
    }


def read_resistance(resistance_path: Path) -> dict:
    if not resistance_path.exists():
        return {
            "has_resistance": False,
            "mutation_count": 0,
            "genes": [],
            "drugs": [],
            "predictions": [],
        }

    rows = []
    with resistance_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows.extend(reader)

    genes = sorted({row["gene"] for row in rows if row.get("gene")})
    predictions = sorted({row["prediction"] for row in rows if row.get("prediction")})
    drugs = sorted(
        {
            drug.strip()
            for row in rows
            for drug in row.get("drugs", "").split(",")
            if drug.strip()
        }
    )
    return {
        "has_resistance": bool(rows),
        "mutation_count": len(rows),
        "genes": genes,
        "drugs": drugs,
        "predictions": predictions,
    }


def read_sample_info(sample_info_path: Path | None, sample_id: str) -> dict:
    result = {
        "ct": None,
        "library_concentration_ng_ul": None,
        "library_fragment_length_bp": None,
    }
    if sample_info_path is None or not sample_info_path.exists():
        return result

    data = json.loads(sample_info_path.read_text(encoding="utf-8"))
    entries = data.values() if isinstance(data, dict) else data
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        if str(entry.get("clarity_sample_id", "")).strip() != sample_id:
            continue
        result["ct"] = maybe_float(entry.get("CT"))
        result["library_concentration_ng_ul"] = maybe_float(entry.get("Library concentration (ng/ul)"))
        result["library_fragment_length_bp"] = maybe_int(entry.get("Library fragment length (bp)"))
        return result

    return result


def build_output_paths(results_dir: Path, sample_id: str) -> dict:
    path_names = {
        "main_fasta": f"{sample_id}.fasta",
        "main_cram": f"{sample_id}.cram",
        "main_blast": f"{sample_id}.fasta.blast",
        "iupac_fasta": f"{sample_id}-0.15-iupac.fasta",
        "iupac_cram": f"{sample_id}-0.15-iupac.cram",
        "iupac_report": f"{sample_id}-0.15-iupac.report.tsv",
        "coverage_tsv": f"{sample_id}-coverage.tsv",
        "vadr_gff": f"{sample_id}.vadr.pass_mod.gff",
        "vadr_bed": f"{sample_id}.vadr.bed",
        "resistance_tsv": f"{sample_id}_resistance.tsv",
        "resistance_bed": f"{sample_id}_resistance.bed",
        "resistance_gff": f"{sample_id}_resistance.gff",
        "resistance_by_drug_tsv": f"{sample_id}_resistance_by_drug.tsv",
        "filtered_vcf_m015": f"{sample_id}-pilon-m0.15.vcf.gz",
        "filtered_vcf_m015_stats": f"{sample_id}-pilon-m0.15.vcf.gz.stats",
        "hostile_json": "hostile.json",
    }
    return {
        key: value if (results_dir / value).exists() else None
        for key, value in path_names.items()
    }


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    sample_id = args.sample_id

    blast_info = read_first_blast_hit(results_dir / f"{sample_id}.fasta.blast")
    report = read_report(results_dir / f"{sample_id}-0.15-iupac.report.tsv")
    coverage_thresholds = read_coverage_tsv(results_dir / f"{sample_id}-coverage.tsv")
    af_counts = read_af_counts(results_dir, sample_id)
    host_filter = read_hostile(results_dir / "hostile.json")
    resistance = read_resistance(results_dir / f"{sample_id}_resistance.tsv")
    sample_info = read_sample_info(Path(args.sample_info) if args.sample_info else None, sample_id)

    summary = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "run_name": args.run_name,
        "sample_id": sample_id,
        "lid": args.lid or None,
        "typing": {
            **blast_info,
            "report_subtype": report.get("subtype"),
            "report_reference": report.get("reference"),
        },
        "qc": {
            "coverage_pct": maybe_float(report.get("coverage")),
            "mean_depth": maybe_float(report.get("meandepth")),
            "coverage_thresholds_pct": coverage_thresholds,
        },
        "variants": {
            "snps": maybe_int(report.get("snps")),
            "multiallelic_snps": maybe_int(report.get("multiallelic_snps")),
            "af_counts": af_counts,
        },
        "host_filter": host_filter,
        "sample_metadata": sample_info,
        "resistance": resistance,
        "outputs": build_output_paths(results_dir, sample_id),
    }

    output_path = Path(args.output)
    output_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
