"""Generate and validate consensus sidecars; FASTA headers contain identifiers only."""
import argparse
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path


def sequence(path):
    lines = Path(path).read_text().splitlines()
    if not lines or not lines[0].startswith(">") or any(line.startswith(">") for line in lines[1:]):
        raise ValueError("Consensus requires one FASTA record")
    return lines[0][1:].split()[0], "".join(lines[1:]).upper()


def sequence_checksum(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


def load(fasta, path=None, required=True):
    path = Path(path or str(fasta)+".consensus.json")
    if not path.exists() and not required:
        return None
    if not path.exists():
        raise ValueError("Missing consensus sidecar; regenerate consensus: "+str(path))
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError("Consensus sidecar must be a JSON object")
    identifier, seq = sequence(fasta)
    if data.get("schema_version") != "1.0" or data.get("converter_version") != "4.0":
        raise ValueError("Unsupported consensus provenance; regenerate consensus")
    if data.get("contig") != identifier or data.get("length") != len(seq) or data.get("sequence_sha256") != sequence_checksum(seq):
        raise ValueError("Consensus sidecar does not match FASTA identity, length or sequence checksum")
    if not re.fullmatch(r"[0-9a-f]{64}", data.get("coordinate_source_sha256", "")):
        raise ValueError("Invalid coordinate-source checksum")
    depths = data.get("usable_depths")
    if not isinstance(depths, list) or len(depths) != len(seq) or any(type(d) is not int or d < 0 for d in depths):
        raise ValueError("Invalid consensus usable depths")
    reasons = data.get("masking_reasons", {})
    if not isinstance(reasons, dict) or any(not str(p).isdigit() or not 1 <= int(p) <= len(seq) for p in reasons):
        raise ValueError("Invalid consensus masking coordinates")
    if data.get("mode") not in {"scaffold", "iupac", "analytical_majority"}:
        raise ValueError("Invalid consensus mode")
    if data.get("mode") != "scaffold" and any(depth == 0 and base != "N" for depth, base in zip(depths, seq)):
        raise ValueError("Analytical consensus contains a base without usable evidence")
    return data


def generate(vcf, reference, output, identifier, mode, min_af=0.15, reference_metadata=None):
    if mode not in {"scaffold", "iupac", "analytical_majority"} or not 0 < min_af <= 1:
        raise ValueError("Invalid consensus mode or allele cutoff")
    reference_id, ref = sequence(reference)
    parent = load(reference, reference_metadata, required=reference_metadata is not None) if mode != "scaffold" else None
    if parent and parent["mode"] != "scaffold":
        raise ValueError("Analytical consensus requires a scaffold coordinate source")
    minimum = 7 if mode == "scaffold" else 1
    awk_mode = "majority" if mode == "scaffold" else mode
    with tempfile.TemporaryDirectory(prefix="consensus-evidence-") as temporary:
        evidence = Path(temporary)/"depths.tsv"
        stats = Path(temporary)/"stats.tsv"
        result = subprocess.run(["awk", "-v", "MODE="+awk_mode, "-v", "MIN_AF="+str(min_af), "-v", "MIN_DP="+str(minimum), "-v", "EVIDENCE_FILE="+str(evidence), "-v", "STATS_FILE="+str(stats), "-f", str(Path(__file__).with_name("vcf_to_iupac.awk")), str(vcf), str(reference)], check=True, capture_output=True, text=True)
        indels, ties = map(int, stats.read_text().split())
        seq = "".join(result.stdout.splitlines()[1:]).upper()
        depths, reasons = [], {}
        for line in evidence.read_text().splitlines():
            position, depth, reason = line.split("\t")
            if int(position) != len(depths)+1:
                raise ValueError("Consensus evidence coordinates are not contiguous")
            depths.append(int(depth))
            if reason:
                reasons[position] = reason
    if len(seq) != len(ref) or len(depths) != len(ref):
        raise ValueError("Consensus conversion changed coordinate length")
    data = {
        "schema_version": "1.0", "converter_version": "4.0",
        "contig": identifier, "length": len(seq), "sequence_sha256": sequence_checksum(seq),
        "mode": mode, "min_af": min_af if mode == "iupac" else None,
        "minimum_usable_depth": minimum, "low_depth_range": [1, 6],
        "reference_fill": mode == "scaffold", "usable_depths": depths,
        "indel_records": indels, "allele_count_ties": ties,
        "masking_reasons": reasons,
        "input_reference": {"contig": reference_id, "sequence_sha256": sequence_checksum(ref)},
        "coordinate_source_sha256": parent["coordinate_source_sha256"] if parent else sequence_checksum(ref),
    }
    Path(output).write_text(">"+identifier+"\n"+"\n".join(seq[i:i+80] for i in range(0,len(seq),80))+"\n")
    Path(str(output)+".consensus.json").write_text(json.dumps(data, indent=2)+"\n")
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--reference-metadata")
    parser.add_argument("--output", required=True)
    parser.add_argument("--identifier", required=True)
    parser.add_argument("--mode", choices=["scaffold", "iupac", "analytical_majority"], required=True)
    parser.add_argument("--min-af", type=float, default=0.15)
    args = parser.parse_args()
    generate(args.vcf, args.reference, args.output, args.identifier, args.mode, args.min_af, args.reference_metadata)


if __name__ == "__main__":
    main()
