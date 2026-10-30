"""Display-copy regressions: analytical evidence and source tracks never change."""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from consensus_metadata import generate
from prepare_igv_review import ANNOTATIONS, THRESHOLDS, checksum, prepare
from build_qc_summary import build_output_paths
import recall_resistance as recall


def fixture(root, sample="ANON"):
    ref = root / "input.fa"
    ref.write_text(">input\nACGTACGT\n")
    calls = root / "calls.vcf"
    calls.write_text("input\t1\t.\tA\tG\t.\t.\t.\tAD:DP\t1,1:2\n")
    scaffold = root / f"{sample}.fasta"
    generate(calls, ref, scaffold, sample, "scaffold")
    calls.write_text(calls.read_text().replace("input\t", sample+"\t"))
    analytical = root / f"{sample}-0.15-iupac.fasta"
    generate(calls, scaffold, analytical, sample+"-0.15-iupac", "iupac", reference_metadata=str(scaffold)+".consensus.json")
    pysam.faidx(str(analytical))
    cram = analytical.with_suffix(".cram")
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": sample+"-0.15-iupac", "LN": 8,
              "M5": hashlib.md5(b"RNNNNNNN").hexdigest()}]}
    with pysam.AlignmentFile(str(cram), "wc", header=header, reference_filename=str(analytical), format_options=[b"embed_ref=1"]) as handle:
        pass
    pysam.index(str(cram))
    for key, suffix in ANNOTATIONS.items():
        (root / (sample+suffix)).write_text(sample+"\t0\t1\n" if suffix.endswith("bed") else f"##gff-version 3\n##sequence-region {sample} 1 8\n{sample}\tfixture\tgene\t1\t3\t.\t+\t.\tName=ANON;sample_codon=RNN\n")
    (root / f"{sample}.vadr.pass_mod.gff").write_text(f"{sample}\tfixture\tgene\t1\t8\t.\t+\t.\tName=NS3\n")
    for threshold in THRESHOLDS:
        vcf = root / f"{sample}-pilon-m0.{threshold[1:]}.vcf"
        vcf.write_text(f"##fileformat=VCFv4.2\n##reference={sample}.fasta\n##contig=<ID={sample},length=8>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n{sample}\t1\t.\tA\tG\t42\tPASS\tDP=2\tDP:AD\t2:1,1\n")
        pysam.tabix_compress(str(vcf), str(vcf)+".gz", force=True)
        pysam.tabix_index(str(vcf)+".gz", preset="vcf", csi=True)


class ReviewTests(unittest.TestCase):
    def test_display_copies_preserve_payloads_and_publish_paths(self):
        for run in ("RUN_A", "RUN_B"):
            with self.subTest(run=run), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                fixture(root)
                before = {p.name: checksum(p) for p in root.iterdir()}
                with patch("prepare_igv_review.tempfile.TemporaryDirectory", wraps=tempfile.TemporaryDirectory) as scratch:
                    manifest, outputs = prepare(root, "ANON")
                    self.assertEqual(scratch.call_args.kwargs["dir"], root)
                for name, sha in before.items():
                    self.assertEqual(checksum(root / name), sha)
                for track in manifest["tracks"]:
                    source, display = root / track["source"], root / track["display"]
                    if track["kind"] == "vcf":
                        with pysam.VariantFile(str(source)) as a, pysam.VariantFile(str(display)) as b:
                            original, shown = list(a), list(b)
                            self.assertEqual([str(r).split("\t")[1:] for r in original], [str(r).split("\t")[1:] for r in shown])
                            self.assertEqual(tuple(b.header.contigs), ("ANON-0.15-iupac",))
                            self.assertIn("##reference=ANON.fasta", str(b.header))
                        with pysam.VariantFile(str(display)) as handle:
                            self.assertEqual(len(list(handle.fetch("ANON-0.15-iupac", 0, 2))), 1)
                    else:
                        self.assertEqual(display.read_text().replace("ANON-0.15-iupac", "ANON"), source.read_text())
                self.assertEqual(build_output_paths(root, "ANON", None)["igv_iupac_cram"], "ANON-0.15-iupac.cram")
                qc = root / "ANON_qc_summary.json"
                subprocess.run([sys.executable, str(ROOT / "scripts/build_qc_summary.py"), "--results-dir", str(root), "--sample-id", "ANON", "--run-name", run, "--output", str(qc)], check=True, capture_output=True)
                self.assertEqual(json.loads(qc.read_text())["igv_review"], manifest)
                self.assertTrue(all((root / p).exists() for p in outputs.values()))

    def test_empty_annotations_are_valid(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture(root)
            (root / "ANON_resistance.gff").write_text("##gff-version 3\n")
            (root / "ANON_resistance.bed").write_text("")
            prepare(root, "ANON")

    def test_incompatible_inputs_fail_before_manifest_publication(self):
        for failure in ("origin", "length", "multiple", "cram", "cram_dictionary", "cram_checksum", "track", "vcf", "index"):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                fixture(root)
                meta = root / "ANON-0.15-iupac.fasta.consensus.json"
                data = json.loads(meta.read_text())
                if failure == "origin":
                    data["coordinate_source_sha256"] = "b"*64
                elif failure == "length":
                    data["length"] += 1
                elif failure == "multiple":
                    with (root / "ANON-0.15-iupac.fasta").open("a") as handle:
                        handle.write(">extra\nA\n")
                elif failure == "cram":
                    (root / "ANON-0.15-iupac.cram").unlink()
                elif failure in {"cram_dictionary", "cram_checksum"}:
                    name = "OTHER" if failure == "cram_dictionary" else "ANON-0.15-iupac"
                    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": name, "LN": 8, "M5": "0"*32}]}
                    with pysam.AlignmentFile(str(root / "ANON-0.15-iupac.cram"), "wc", header=header) as handle:
                        pass
                elif failure == "track":
                    (root / "ANON_resistance.bed").write_text("OTHER\t0\t1\n")
                elif failure == "vcf":
                    source = root / "ANON-pilon-m0.05.vcf"
                    source.write_text(source.read_text().replace("ANON\t", "OTHER\t"))
                    pysam.tabix_compress(str(source), str(source)+".gz", force=True)
                elif failure == "index":
                    (root / "ANON-0.15-iupac.cram.crai").unlink()
                meta.write_text(json.dumps(data))
                with self.assertRaises((ValueError, OSError)):
                    prepare(root, "ANON")
                self.assertFalse((root / "ANON-0.15-iupac.review.json").exists())

    def test_review_updates_participate_in_recall_rollback(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "archive"
            archive.mkdir()
            fixture(archive)
            _, outputs = prepare(archive, "ANON")
            target = archive / outputs["igv_iupac_resistance_gff"]
            old = target.read_bytes()
            staged = root / "staged"
            staged.mkdir()
            override = staged / "new.gff"
            override.write_text("ANON\tfixture\tgene\t2\t4\t.\t+\t.\tName=changed\n")
            prepare(archive, "ANON", output_dir=staged, overrides={"resistance_gff": override})
            replacements = {archive / p.name: p for p in staged.iterdir() if p.name.startswith("ANON-")}
            backup = root / "backup"
            recall.apply_transaction(archive, replacements, backup)
            self.assertIn("changed", target.read_text())
            recall.restore(backup, archive)
            self.assertEqual(target.read_bytes(), old)

    def test_recall_stages_review_tracks_and_qc_in_same_transaction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, staged = root / "archive", root / "staged"
            archive.mkdir()
            fixture(archive)
            manifest, outputs = prepare(archive, "ANON")
            qc_path = archive / "ANON_qc_summary.json"
            qc_path.write_text(json.dumps({"sample_id": "ANON", "outputs": outputs, "igv_review": manifest}))
            blast = archive / "ANON.blast"
            blast.write_text("ANON\t1a_HQ850279\t100\t4\t0\t0\t1\t4\t1\t4\t0\t8\n")

            def annotate(command, **kwargs):
                for suffix in recall.SUFFIXES:
                    source = archive / ("ANON"+suffix)
                    target = staged / source.name
                    if source.exists():
                        shutil.copy2(source, target)
                    else:
                        target.write_text("")
                (staged / "ANON_resistance.json").write_text(json.dumps({"schema_version": "3.1", "drug_outcomes": [{"drug": "Dasabuvir", "outcome": "susceptible"}]}))
                (staged / "ANON_resistance.gff").write_text("ANON\tfixture\tgene\t2\t4\t.\t+\t.\tName=recalled\n")

            args = SimpleNamespace(rules=ROOT / "assets/hcv_geno2pheno_rules.csv", h77_fasta=ROOT / "refgenomes/resistance/AF011751.1.fasta")
            item = ("ANON", archive / "ANON-0.15-iupac.fasta", blast, archive / "ANON.vadr.pass_mod.gff", archive / "ANON-0.15-iupac.cram")
            with patch.object(recall.subprocess, "run", side_effect=annotate):
                generated, _, qc = recall.stage_sample(archive, item, args, staged)
            review = staged / qc["outputs"]["igv_iupac_resistance_gff"]
            self.assertIn(review, generated)
            self.assertIn("recalled", review.read_text())
            track = next(t for t in qc["igv_review"]["tracks"] if t["source_key"] == "resistance_gff")
            self.assertEqual(track["source_sha256"], checksum(staged / "ANON_resistance.gff"))
            self.assertNotIn("recalled", (archive / "ANON_resistance.gff").read_text())


if __name__ == "__main__":
    unittest.main()
