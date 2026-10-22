#!/usr/bin/env python3
"""Make coordinate-preserving display copies for the analytical IUPAC IGV view."""
import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import pysam
import pysam.bcftools

from consensus_metadata import load, sequence


ANNOTATIONS = {
    "resistance_gff": "_resistance.gff",
    "resistance_bed": "_resistance.bed",
    "resistance_sites_gff3": "_resistance_sites.gff3",
    "vadr_bed": ".vadr.bed",
    "coverage_1x_bed": "-coverage-1x.bed",
}
THRESHOLDS = ("001", "005", "01", "015", "02", "03", "04")


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rename_annotation(source, destination, old, new, length):
    lines = []
    bed = source.suffix == ".bed"
    for line in source.read_text().splitlines(keepends=True):
        if line.startswith("##FASTA"):
            raise ValueError("Embedded FASTA is not supported in review annotations")
        if line.startswith("##sequence-region"):
            fields = line.split()
            if len(fields) != 4 or fields[1] != old or not 1 <= int(fields[2]) <= int(fields[3]) <= length:
                raise ValueError(f"Invalid sequence-region in {source}")
            line = f"##sequence-region {new} {fields[2]} {fields[3]}\n"
        elif line.strip() and not line.startswith("#"):
            fields = line.rstrip("\r\n").split("\t")
            minimum = 3 if bed else 9
            if len(fields) < minimum or fields[0] != old:
                raise ValueError(f"Unexpected annotation contig/record in {source}")
            start, end = map(int, fields[1:3] if bed else fields[3:5])
            if not (0 if bed else 1) <= start <= end <= length:
                raise ValueError(f"Out-of-range annotation in {source}")
            line = new + line[len(old):]
        lines.append(line)
    destination.write_text("".join(lines))


def prepare(directory, sample, output_dir=None, overrides=None):
    directory = Path(directory)
    output_dir = Path(output_dir or directory)
    overrides = overrides or {}
    scaffold = directory / f"{sample}.fasta"
    analytical = directory / f"{sample}-0.15-iupac.fasta"
    sm, im = load(scaffold), load(analytical)
    old, source_seq = sequence(scaffold)
    new, display_seq = sequence(analytical)
    if (sm["mode"] != "scaffold" or im["mode"] != "iupac" or im["min_af"] != 0.15
            or im["minimum_usable_depth"] != 1 or im["reference_fill"]
            or sm["coordinate_source_sha256"] != im["coordinate_source_sha256"]
            or len(source_seq) != len(display_seq)):
        raise ValueError("IGV review requires coordinate-compatible scaffold and analytical 15% IUPAC")
    cram = analytical.with_suffix(".cram")
    for path in (Path(str(analytical)+".fai"), cram, Path(str(cram)+".crai")):
        if not path.is_file():
            raise ValueError(f"Missing review reference/read file: {path}")
    with pysam.FastaFile(str(analytical)) as handle:
        if tuple(handle.references) != (new,) or tuple(handle.lengths) != (len(display_seq),) or handle.fetch(new).upper() != display_seq:
            raise ValueError("IUPAC FASTA index does not match review reference")
    with pysam.AlignmentFile(str(cram), "rc", reference_filename=str(analytical)) as handle:
        if handle.references != (new,) or handle.lengths != (len(display_seq),):
            raise ValueError("IUPAC CRAM sequence dictionary does not match review reference")
        sq = handle.header.to_dict()["SQ"][0]
        if sq.get("M5") != hashlib.md5(display_seq.encode()).hexdigest():
            raise ValueError("IUPAC CRAM reference checksum does not match review reference")
        handle.check_index()
    prefix = f"{sample}-0.15-iupac.review."
    outputs = {"igv_iupac_fasta": analytical.name, "igv_iupac_fasta_index": analytical.name+".fai",
               "igv_iupac_cram": cram.name, "igv_iupac_cram_index": cram.name+".crai"}
    tracks = []
    sources = {key: directory / f"{sample}{suffix}" for key, suffix in ANNOTATIONS.items()}
    vadr = [p for p in directory.glob(f"{sample}.vadr.*_mod.gff")]
    if len(vadr) != 1:
        raise ValueError("Review requires exactly one selected normalized VADR GFF")
    sources["selected_vadr_gff"] = vadr[0]
    for threshold in THRESHOLDS:
        sources[f"filtered_vcf_m{threshold}"] = directory / f"{sample}-pilon-m0.{threshold[1:]}.vcf.gz"
    # Build in isolation: invalid input never publishes a partial new manifest.
    # Keep staging on the destination filesystem: Hopper work/results may be
    # NFS while the system temporary directory is node-local.
    with tempfile.TemporaryDirectory(prefix=".igv-review-", dir=output_dir) as temporary:
        staging = Path(temporary)
        mapping = staging / "contigs.tsv"
        mapping.write_text(f"{old}\t{new}\n")
        headers = staging / "headers.txt"
        headers.write_text('##igv_review="Contig-name-only display copy; REF and ALT retain original scaffold-derived calling reference"\n')
        for key, original in sources.items():
            source = Path(overrides.get(key, original))
            if not source.is_file():
                raise ValueError(f"Missing review source: {source}")
            suffix = original.name[len(sample):].lstrip("_-.")
            name = prefix+suffix
            target = staging / name
            is_vcf = source.name.endswith(".vcf.gz")
            if is_vcf:
                with pysam.VariantFile(str(source)) as handle:
                    if tuple(handle.header.contigs) != (old,) or handle.header.contigs[old].length != len(source_seq):
                        raise ValueError(f"Unexpected VCF sequence dictionary: {source}")
                    for record in handle:
                        if record.contig != old or not 0 <= record.start < record.stop <= len(source_seq):
                            raise ValueError(f"Unexpected VCF contig/coordinates: {source}")
                        if any("[" in a or "]" in a for a in record.alts or ()):
                            raise ValueError("Breakend alleles cannot be contig-renamed as review tracks")
                # Use the bcftools bundled with pysam: this also works inside the
                # pipeline's Python container, without launching nested containers.
                pysam.bcftools.annotate("--rename-chrs", str(mapping), "--header-lines", str(headers),
                                       "--no-version", "-Oz", "-o", str(target), str(source), catch_stdout=False)
                pysam.bcftools.index("--csi", str(target), catch_stdout=False)
            else:
                rename_annotation(source, target, old, new, len(source_seq))
            output_key = "igv_iupac_"+key
            outputs[output_key] = name
            if is_vcf:
                outputs[output_key+"_index"] = name+".csi"
            tracks.append({"source_key": key, "output_key": output_key, "source": original.name,
                           "display": name, "source_sha256": checksum(source), "display_sha256": checksum(target),
                           "kind": "vcf" if is_vcf else "annotation",
                           "analytical_source": "iupac" if key.startswith("resistance") else "scaffold-derived"})
        manifest = {"schema_version": "1.0", "view": "analytical_iupac", "sample_id": sample,
                    "coordinate_source_sha256": im["coordinate_source_sha256"],
                    "source_reference": {"path": scaffold.name, "contig": old, "sequence_sha256": sm["sequence_sha256"]},
                    "display_reference": {"path": analytical.name, "contig": new, "sequence_sha256": im["sequence_sha256"], "length": len(display_seq)},
                    "tracks": tracks, "outputs": outputs.copy()}
        manifest_name = prefix+"json"
        (staging / manifest_name).write_text(json.dumps(manifest, indent=2)+"\n")
        for path in sorted(staging.iterdir(), key=lambda p: p.name == manifest_name):
            if path.name.startswith(prefix):
                path.replace(output_dir / path.name)
    outputs["igv_review_manifest"] = manifest_name
    return manifest, outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    args = parser.parse_args()
    prepare(args.results_dir, args.sample_id)


if __name__ == "__main__":
    main()
