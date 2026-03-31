#!/usr/bin/env python3
"""Aggregate per-sample QC summaries into run-level JSON and JSONL."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate sample QC summary JSON files")
    parser.add_argument("--output-json", required=True, help="Combined JSON output path")
    parser.add_argument("--output-jsonl", required=True, help="Combined JSONL output path")
    parser.add_argument("inputs", nargs="+", help="Sample QC summary JSON files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = [
        json.loads(Path(path).read_text(encoding="utf-8"))
        for path in args.inputs
    ]
    records.sort(key=lambda item: (item.get("run_name") or "", item.get("sample_id") or ""))

    Path(args.output_json).write_text(json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with Path(args.output_jsonl).open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
