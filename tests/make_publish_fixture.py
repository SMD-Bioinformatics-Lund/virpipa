"""Upgrade the legacy publication fixture with genuine synthetic consensus evidence."""
import shutil
import sys
from pathlib import Path

import pysam

from test_igv_review import fixture
from consensus_metadata import generate


def main():
    repository, output = map(Path, sys.argv[1:])
    shutil.copytree(repository / "assets/test_data/qc_summary/results", output)
    fixture(output, "SAMPLE001")
    pysam.faidx(str(output / "SAMPLE001.fasta"))
    generate(output / "calls.vcf", output / "SAMPLE001.fasta", output / "SAMPLE001-majority.fasta",
             "SAMPLE001-majority", "analytical_majority", reference_metadata=str(output / "SAMPLE001.fasta.consensus.json"))
    pysam.faidx(str(output / "SAMPLE001-majority.fasta"))
    with pysam.AlignmentFile(str(output / "SAMPLE001.cram"), "wc", reference_filename=str(output / "SAMPLE001.fasta"),
                             header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "SAMPLE001", "LN": 8}]},
                             format_options=[b"embed_ref=1"]) as handle:
        pass
    pysam.index(str(output / "SAMPLE001.cram"))


if __name__ == "__main__":
    main()
