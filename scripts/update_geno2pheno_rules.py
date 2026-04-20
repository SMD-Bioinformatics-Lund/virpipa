#!/usr/bin/env python3
"""
Download and normalize HCV geno2pheno resistance rules.

This script can either:
1. Download the rules page directly from geno2pheno, or
2. Rebuild the normalized JSON from an existing CSV snapshot.

The normalized JSON is the pipeline-facing artifact.
"""

import argparse
import csv
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_URL = "https://hcv.geno2pheno.org/index.php?page=Rules"
DEFAULT_JSON = "hcv_geno2pheno_rules.json"


def parse_subtype_pattern(pattern):
    """Expand geno2pheno subtype selectors into explicit subtype strings."""
    pattern = (pattern or "").strip()
    if not pattern:
        return []

    expansions = {
        "1": ["", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"],
        "2": ["", "a", "b", "c", "d", "e", "f", "g", "k", "l", "m"],
        "3": ["", "a", "b", "c", "d", "e", "f", "g", "h", "i", "k", "l", "m"],
        "4": ["", "a", "b", "c", "d", "e", "f", "g", "h", "i", "k", "l", "m", "n", "o", "p", "q", "r"],
        "5": ["", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q"],
        "6": ["", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "t", "u", "v", "w", "x"],
    }

    subtypes = set()
    for part in [value.strip() for value in pattern.split(",") if value.strip()]:
        if part in expansions:
            subtypes.update(f"{part}{suffix}" for suffix in expansions[part])
        else:
            subtypes.add(part)

    return sorted(subtypes)


def parse_rule_definition(rule_definition):
    """Split one rule definition into its amino-acid parts."""
    rule_definition = (rule_definition or "").strip()
    if not rule_definition:
        return []

    parsed = []
    for part in re.split(r"\s+and\s+", rule_definition, flags=re.IGNORECASE):
        match = re.match(r"(\d+)\s*(del|[A-Za-z*]+)?", part.strip())
        if not match:
            continue
        parsed.append(
            {
                "position": int(match.group(1)),
                "aa": match.group(2) or "",
                "raw": part.strip(),
            }
        )

    return parsed


def fetch_rules_html(url):
    import requests

    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text


def extract_rules_rows_from_html(html):
    from bs4 import BeautifulSoup

    soup = BeautifulSoup(html, "html.parser")
    table = soup.find("table", {"class": "hbv_result"})
    if table is None:
        raise ValueError("Could not find geno2pheno rules table with class 'hbv_result'")

    rows = []
    header = None
    for row in table.find_all("tr"):
        cells = row.find_all(["td", "th"])
        values = []
        for cell in cells:
            divs = cell.find_all("div")
            text = divs[0].get_text(strip=True) if divs else cell.get_text(strip=True)
            values.append(text)

        if not values:
            continue

        if header is None:
            header = values
            continue

        if len(values) == len(header):
            rows.append(dict(zip(header, values)))

    if header is None or not rows:
        raise ValueError("No geno2pheno rule rows were extracted from the rules table")

    return header, rows


def load_rules_rows_from_csv(csv_path):
    lines = Path(csv_path).read_text(encoding="utf-8").splitlines()
    header = None
    data_lines = []
    for line in lines:
        if line.startswith("#"):
            continue
        if header is None:
            header = next(csv.reader([line]))
            continue
        data_lines.append(line)

    if header is None:
        raise ValueError(f"Could not find a CSV header in {csv_path}")

    return header, list(csv.DictReader(data_lines, fieldnames=header))


def normalize_rows(rows, source_url=None, retrieved_at=None):
    normalized_rules = []

    for row in rows:
        if not row.get("Drugs"):
            continue

        subtype_pattern = (row.get("Rule for Subtype") or "").strip()
        rule_definition = (row.get("Rule Definitions") or "").strip()

        base_rule = {
            "drug": (row.get("Drugs") or "").strip(),
            "region": (row.get("Region") or "").strip(),
            "rule_definition": rule_definition,
            "subtype_pattern": subtype_pattern,
            "drug_licensed_for_genotype": (row.get("Drug licensed for genotype") or "").strip(),
            "prediction": (row.get("Predictions") or "").strip(),
            "reference": (row.get("Reference") or "").strip(),
        }

        normalized_rules.append(base_rule)

    columns = [
        "drug",
        "region",
        "rule_definition",
        "subtype_pattern",
        "drug_licensed_for_genotype",
        "prediction",
        "reference",
    ]

    return {
        "metadata": {
            "source_url": source_url,
            "retrieved_at": retrieved_at,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "total_rules": len(normalized_rules),
            "drugs": sorted({rule["drug"] for rule in normalized_rules}),
            "regions": sorted({rule["region"] for rule in normalized_rules}),
        },
        "columns": columns,
        "rules": [[rule[column] for column in columns] for rule in normalized_rules],
    }


def write_csv(header, rows, output_csv, source_url=None, retrieved_at=None):
    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        handle.write(f"# Downloaded from: {source_url or DEFAULT_URL}\n")
        if retrieved_at:
            handle.write(f"# Last updated: {retrieved_at}\n")
        handle.write("#\n")
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow([(row.get(column) or "").strip() for column in header])


def write_json(data, output_json):
    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_rules_json(json_path):
    data = json.loads(Path(json_path).read_text(encoding="utf-8"))
    columns = data.get("columns")
    if columns and data.get("rules"):
        data["rules"] = [dict(zip(columns, row)) for row in data["rules"]]
    return data


def main():
    parser = argparse.ArgumentParser(description="Download and normalize HCV geno2pheno resistance rules")
    parser.add_argument("--url", default=DEFAULT_URL, help=f"Rules page URL (default: {DEFAULT_URL})")
    parser.add_argument("--input-csv", help="Reuse an existing CSV snapshot instead of downloading")
    parser.add_argument("--output-json", default=DEFAULT_JSON, help=f"Normalized JSON output (default: {DEFAULT_JSON})")
    parser.add_argument("--output-csv", help="Optional CSV snapshot output path")
    args = parser.parse_args()

    try:
        if args.input_csv:
            header, rows = load_rules_rows_from_csv(args.input_csv)
            source_url = DEFAULT_URL
            retrieved_at = None
        else:
            html = fetch_rules_html(args.url)
            header, rows = extract_rules_rows_from_html(html)
            source_url = args.url
            retrieved_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        if args.output_csv:
            write_csv(header, rows, args.output_csv, source_url=source_url, retrieved_at=retrieved_at)

        normalized = normalize_rows(rows, source_url=source_url, retrieved_at=retrieved_at)
        write_json(normalized, args.output_json)

        print(f"Wrote normalized rules JSON to {args.output_json}")
        if args.output_csv:
            print(f"Wrote CSV snapshot to {args.output_csv}")
        print(f"Rules: {normalized['metadata']['total_rules']}")
        return 0
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
