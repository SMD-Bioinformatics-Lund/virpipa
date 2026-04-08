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
    if isinstance(data, list):
        data = next((item for item in data if isinstance(item, dict)), {})
    elif not isinstance(data, dict):
        data = {}

    return {
        "enabled": True,
        "reads_in": data.get("reads_in"),
        "reads_removed_proportion": data.get("reads_removed_proportion"),
    }

def mutation_label(gene: str | None, ref_aa: str | None, aa_pos: str | int | None, alt_aa: str | None) -> str | None:
    if not gene or not ref_aa or aa_pos in (None, "") or not alt_aa:
        return None
    return f"{gene}:{ref_aa}{aa_pos}{alt_aa}"


def prediction_priority(prediction: str) -> tuple[int, str]:
    value = prediction.strip().lower()
    if "resistant" in value:
        return (0, value)
    if "probable" in value or "intermediate" in value or "reduced" in value:
        return (1, value)
    if "possible" in value:
        return (2, value)
    return (3, value)


def best_prediction(predictions: list[str]) -> str | None:
    cleaned = sorted({item.strip() for item in predictions if item and item.strip()}, key=prediction_priority)
    return cleaned[0] if cleaned else None


def read_resistance(resistance_path: Path, resistance_by_drug_path: Path) -> dict:
    analysis_present = resistance_path.exists() or resistance_by_drug_path.exists()
    if not analysis_present:
        return {
            "analysis_present": False,
            "has_resistance": False,
            "mutation_count": 0,
            "genes": [],
            "drugs": [],
            "predictions": [],
            "by_drug": [],
            "mutations": [],
        }

    rows = []
    if resistance_path.exists():
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

    mutations = []
    for row in rows:
        mutation = {
            "gene": row.get("gene"),
            "aa_pos": maybe_int(row.get("aa_pos")) if row.get("aa_pos") not in (None, "") else None,
            "ref_aa": row.get("ref_aa"),
            "alt_aa": row.get("alt_aa"),
            "mutation_label": mutation_label(row.get("gene"), row.get("ref_aa"), row.get("aa_pos"), row.get("alt_aa")),
            "rule_definition": row.get("rule_definition") or None,
            "prediction": row.get("prediction") or None,
            "drugs": [drug.strip() for drug in row.get("drugs", "").split(",") if drug.strip()],
            "genomic_start": maybe_int(row.get("genomic_start")),
            "genomic_end": maybe_int(row.get("genomic_end")),
            "strand": row.get("strand") or None,
            "reference": row.get("reference") or None,
        }
        mutations.append(mutation)

    derived_by_drug: dict[str, dict[str, object]] = {}
    for mutation in mutations:
        for drug in mutation.get("drugs", []):
            record = derived_by_drug.setdefault(drug, {"drug": drug, "predictions": [], "mutations": []})
            if mutation.get("prediction"):
                record["predictions"].append(mutation["prediction"])
            if mutation.get("mutation_label"):
                record["mutations"].append(mutation["mutation_label"])

    by_drug_records = []
    if resistance_by_drug_path.exists():
        current_drug = None
        current_rows: list[dict] = []
        with resistance_by_drug_path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("## "):
                    if current_drug is not None:
                        by_drug_records.append(
                            {
                                "drug": current_drug,
                                "prediction": best_prediction([item.get("prediction", "") for item in current_rows]),
                                "mutations": [
                                    label
                                    for label in (
                                        mutation_label(item.get("gene"), item.get("ref_aa"), item.get("aa_pos"), item.get("alt_aa"))
                                        for item in current_rows
                                    )
                                    if label
                                ],
                            }
                        )
                    current_drug = line[3:].strip()
                    current_rows = []
                    continue
                if line.startswith("#") or line.startswith("sample\t"):
                    continue
                if current_drug is None:
                    continue
                fields = line.split("\t")
                if len(fields) < 12:
                    continue
                current_rows.append(
                    {
                        "sample": fields[0],
                        "gene": fields[1],
                        "genomic_start": fields[2],
                        "genomic_end": fields[3],
                        "ref_nuc": fields[4],
                        "alt_nuc": fields[5],
                        "aa_pos": fields[6],
                        "ref_aa": fields[7],
                        "alt_aa": fields[8],
                        "rule_definition": fields[9],
                        "prediction": fields[10],
                        "reference": fields[11],
                    }
                )
        if current_drug is not None:
            by_drug_records.append(
                {
                    "drug": current_drug,
                    "prediction": best_prediction([item.get("prediction", "") for item in current_rows]),
                    "mutations": [
                        label
                        for label in (
                            mutation_label(item.get("gene"), item.get("ref_aa"), item.get("aa_pos"), item.get("alt_aa"))
                            for item in current_rows
                        )
                        if label
                    ],
                }
            )
    if not by_drug_records and derived_by_drug:
        by_drug_records = [
            {
                "drug": drug,
                "prediction": best_prediction(record["predictions"]),
                "mutations": sorted(set(record["mutations"])),
            }
            for drug, record in sorted(derived_by_drug.items())
        ]

    return {
        "analysis_present": True,
        "has_resistance": bool(rows),
        "mutation_count": len(rows),
        "genes": genes,
        "drugs": drugs,
        "predictions": predictions,
        "by_drug": by_drug_records,
        "mutations": mutations,
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


def first_existing_relative_path(results_dir: Path, *candidates: str) -> str | None:
    for candidate in candidates:
        if (results_dir / candidate).exists():
            return candidate
    return None


def build_output_paths(results_dir: Path, sample_id: str, lid: str | None) -> dict:
    selected_vadr_gff = first_existing_relative_path(
        results_dir,
        f"{sample_id}.vadr.pass_mod.gff",
        f"{sample_id}.vadr.fail_mod.gff",
    )

    lid_value = lid or ""
    lid_fasta = f"lid/{lid_value}.fasta" if lid_value else None
    lid_iupac_fasta = f"lid/{lid_value}-0.15-iupac.fasta" if lid_value else None
    lid_rug_kde_plot = f"lid/{lid_value}_rug_kde_plot.png" if lid_value else None
    lid_2limsrs = f"lid/{lid_value}-2limsrs.txt" if lid_value else None

    path_names = {
        "main_fasta": f"{sample_id}.fasta",
        "main_fasta_index": f"{sample_id}.fasta.fai",
        "main_cram": f"{sample_id}.cram",
        "main_cram_index": f"{sample_id}.cram.crai",
        "main_blast": f"{sample_id}.fasta.blast",
        "iupac_fasta": f"{sample_id}-0.15-iupac.fasta",
        "iupac_cram": f"{sample_id}-0.15-iupac.cram",
        "iupac_cram_index": f"{sample_id}-0.15-iupac.cram.crai",
        "iupac_report": f"{sample_id}-0.15-iupac.report.tsv",
        "coverage_tsv": f"{sample_id}-coverage.tsv",
        "rug_kde_plot": f"{sample_id}_rug_kde_plot.png",
        "vadr_pass_gff": f"{sample_id}.vadr.pass_mod.gff",
        "vadr_fail_gff": f"{sample_id}.vadr.fail_mod.gff",
        "vadr_bed": f"{sample_id}.vadr.bed",
        "resistance_tsv": f"{sample_id}_resistance.tsv",
        "resistance_bed": f"{sample_id}_resistance.bed",
        "resistance_gff": f"{sample_id}_resistance.gff",
        "resistance_by_drug_tsv": f"{sample_id}_resistance_by_drug.tsv",
        "filtered_vcf_m005": f"{sample_id}-pilon-m0.05.vcf.gz",
        "filtered_vcf_m01": f"{sample_id}-pilon-m0.1.vcf.gz",
        "filtered_vcf_m015": f"{sample_id}-pilon-m0.15.vcf.gz",
        "filtered_vcf_m02": f"{sample_id}-pilon-m0.2.vcf.gz",
        "filtered_vcf_m03": f"{sample_id}-pilon-m0.3.vcf.gz",
        "filtered_vcf_m04": f"{sample_id}-pilon-m0.4.vcf.gz",
        "filtered_vcf_m015_stats": f"{sample_id}-pilon-m0.15.vcf.gz.stats",
        "hostile_json": "hostile.json",
    }
    outputs = {
        key: value if (results_dir / value).exists() else None
        for key, value in path_names.items()
    }
    outputs["lid_fasta"] = lid_fasta if lid_fasta and (results_dir / lid_fasta).exists() else None
    outputs["lid_iupac_fasta"] = lid_iupac_fasta if lid_iupac_fasta and (results_dir / lid_iupac_fasta).exists() else None
    outputs["lid_rug_kde_plot"] = lid_rug_kde_plot if lid_rug_kde_plot and (results_dir / lid_rug_kde_plot).exists() else None
    outputs["lid_2limsrs"] = lid_2limsrs if lid_2limsrs and (results_dir / lid_2limsrs).exists() else None
    outputs["selected_vadr_gff"] = selected_vadr_gff
    outputs["vadr_gff"] = selected_vadr_gff
    outputs["display_rug_kde_plot"] = outputs["lid_rug_kde_plot"] or outputs["rug_kde_plot"]
    outputs["export_fasta"] = outputs["lid_fasta"] or outputs["main_fasta"]
    outputs["export_iupac_fasta"] = outputs["lid_iupac_fasta"] or outputs["iupac_fasta"]
    return outputs


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    sample_id = args.sample_id

    blast_info = read_first_blast_hit(results_dir / f"{sample_id}.fasta.blast")
    report = read_report(results_dir / f"{sample_id}-0.15-iupac.report.tsv")
    coverage_thresholds = read_coverage_tsv(results_dir / f"{sample_id}-coverage.tsv")
    af_counts = read_af_counts(results_dir, sample_id)
    host_filter = read_hostile(results_dir / "hostile.json")
    resistance = read_resistance(
        results_dir / f"{sample_id}_resistance.tsv",
        results_dir / f"{sample_id}_resistance_by_drug.tsv",
    )
    sample_info = read_sample_info(Path(args.sample_info) if args.sample_info else None, sample_id)

    summary = {
        "schema_version": 4,
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "pipeline_name": "virpipa",
        "virus": "HCV",
        "run_name": args.run_name,
        "sample_id": sample_id,
        "sample_run_id": f"{sample_id}_{args.run_name}",
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
        "outputs": build_output_paths(results_dir, sample_id, args.lid or None),
    }

    output_path = Path(args.output)
    output_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
