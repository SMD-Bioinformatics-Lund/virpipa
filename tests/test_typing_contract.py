"""The authoritative assignment survives uncertainty, publication and re-call."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from rank_blast import summarize, validate
from consensus_metadata import generate
from recall_resistance import typing_evidence


class TypingContractTests(unittest.TestCase):
    def test_unknown_typing_emits_all_drugs_unassessed(self):
        for status in ("unassigned", "ambiguous"):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                source, vcf, fasta = root / "ref.fa", root / "calls.vcf", root / "sample.fa"
                source.write_text(">ref\nACGTACGT\n")
                vcf.write_text("ref\t1\t.\tA\tG\t.\t.\t.\tAD:DP\t1,1:2\n")
                generate(vcf, source, fasta, "sample", "iupac", 0.15)
                search = root / "sample.blast"
                search.write_text("")
                evidence = summarize(search, fasta)
                evidence["status"] = status
                typing_path = root / "typing.json"
                typing_path.write_text(json.dumps(evidence))
                annotation = root / "sample.gff"
                annotation.write_text("##gff-version 3\n")
                subprocess.run([sys.executable, str(ROOT / "scripts/annotate_vcf_resistance.py"),
                                "--fasta", str(fasta), "--gff", str(annotation),
                                "--typing-json", str(typing_path), "--rules", str(ROOT / "assets/hcv_geno2pheno_rules.csv"),
                                "--output-dir", str(root)], check=True, capture_output=True, text=True)
                result = json.loads((root / "sample_resistance.json").read_text())
                self.assertIsNone(result["subtype"])
                self.assertEqual(len(result["drug_outcomes"]), 16)
                self.assertTrue(all(r["outcome"] == "not assessed" and r["licensed"] is None
                                    and not r["mutations"] and r["warnings"] for r in result["drug_outcomes"]))
                self.assertFalse(result["rules"])
                self.assertEqual(result["typing_evidence"]["status"], status)

    def test_qc_publishes_new_evidence_and_preserves_legacy_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "results"
            shutil.copytree(ROOT / "assets/test_data/qc_summary/results", root)
            fasta = root / "SAMPLE001-0.15-iupac.fasta"
            search = root / "iupac.blast"
            search.write_text("")
            evidence = summarize(search, fasta)
            (root / "SAMPLE001_typing.json").write_text(json.dumps(evidence))
            output = root / "summary.json"
            subprocess.run([sys.executable, str(ROOT / "scripts/build_qc_summary.py"),
                            "--results-dir", str(root), "--sample-id", "SAMPLE001",
                            "--run-name", "fixture", "--output", str(output)], check=True, capture_output=True)
            result = json.loads(output.read_text())
            self.assertEqual(result["typing"]["evidence"], evidence)
            self.assertIsNotNone(result["typing"]["main_blast_genotype"])
            self.assertEqual(result["outputs"]["typing_json"], "SAMPLE001_typing.json")

    def test_recall_reuses_provenance_but_rejects_stale_assignment(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta, search = root / "sample.fa", root / "sample.blast"
            fasta.write_text(">sample\nACGT\n")
            search.write_text("sample\t1a_ref\t100\t4\t0\t0\t1\t4\t1\t4\t0\t8\n")
            evidence = summarize(search, fasta)
            evidence["database"] = {"db.nsq": "a" * 64}
            (root / "sample_typing.json").write_text(json.dumps(evidence))
            self.assertEqual(typing_evidence(search, fasta, "sample"), evidence)
            search.write_text("")
            with self.assertRaisesRegex(ValueError, "differs"):
                typing_evidence(search, fasta, "sample")
            fasta.write_text(">sample\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "mismatched"):
                validate(evidence, fasta)

    def test_invalid_declared_assignment_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta, search = root / "sample.fa", root / "sample.blast"
            fasta.write_text(">sample\nACGT\n")
            search.write_text("")
            evidence = summarize(search, fasta)
            evidence["subtype"] = "1a"
            with self.assertRaises(ValueError):
                validate(evidence, fasta)
