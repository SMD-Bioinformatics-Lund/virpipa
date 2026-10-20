"""Final analytical FASTAs preserve nonzero evidence without scaffold fill-in."""
import json
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
from consensus_metadata import generate, load, sequence


class ConsensusPolicyTests(unittest.TestCase):
    def convert(self, depth, counts, alt="G", mode="iupac"):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/"ref.fa").write_text(">ref\nAC\n")
            (root/"calls.vcf").write_text(f"ref\t1\t.\tA\t{alt}\t.\t.\t.\tAD:DP\t{counts}:{depth}\n")
            data = generate(root/"calls.vcf", root/"ref.fa", root/"out.fa", "ANON", mode)
            self.assertEqual((root/"out.fa").read_text().splitlines()[0], ">ANON")
            self.assertEqual(load(root/"out.fa"), data)
            return sequence(root/"out.fa")[1], data

    def test_nonzero_depth_alt_is_not_replaced_by_reference(self):
        for depth in (1, 6, 7):
            for mode in ("iupac", "analytical_majority"):
                with self.subTest(depth=depth, mode=mode):
                    seq, data = self.convert(depth, f"0,{depth}", mode=mode)
                    self.assertEqual(seq, "GN")
                    self.assertEqual(data["usable_depths"], [depth, 0])

    def test_low_depth_disagreement_is_preserved_as_iupac(self):
        self.assertEqual(self.convert(6, "5,1")[0], "RN")
        self.assertEqual(self.convert(2, "1,1")[0], "RN")
        self.assertEqual(self.convert(7, "6,1")[0], "AN")

    def test_zero_depth_and_missing_ad_are_not_reference_calls(self):
        seq, data = self.convert(0, "0,0")
        self.assertEqual(seq, "NN")
        self.assertEqual(data["masking_reasons"]["1"], "missing_coverage")
        seq, data = self.convert(5, ".")
        self.assertEqual(seq, "NN")
        self.assertEqual(data["masking_reasons"]["1"], "invalid_evidence")

    def test_all_qualifying_alleles_and_unqualified_ref(self):
        self.assertEqual(self.convert(100, "5,95")[0], "GN")
        self.assertEqual(self.convert(100, "60,20,20", "G,T")[0], "DN")
        self.assertEqual(self.convert(100, "85,15")[0], "RN")
        seq, data = self.convert(4, "1,1,1,1", "C,G,T")
        self.assertEqual(seq, "NN")
        self.assertEqual(data["usable_depths"][0], 4)
        self.assertEqual(data["masking_reasons"]["1"], "all_four_alleles")

    def test_scaffold_retention_is_separate_and_majority_ties_are_deterministic(self):
        self.assertEqual(self.convert(6, "0,6", mode="scaffold")[0], "AC")
        self.assertEqual(self.convert(7, "0,7", mode="scaffold")[0], "GC")
        self.assertEqual(self.convert(2, "1,1", mode="analytical_majority")[0], "AN")
        self.assertEqual(self.convert(2, "0,1,1", "G,T", "analytical_majority")[0], "GN")

    def test_indels_are_masked_without_changing_coordinates(self):
        seq, data = self.convert(8, "0,8", "AG")
        self.assertEqual(seq, "NN")
        self.assertEqual(data["masking_reasons"]["1"], "indel_ambiguity")
        self.assertEqual(data["indel_records"], 1)

    def test_wrong_vcf_contig_cannot_silently_produce_missing_coverage(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/"ref.fa").write_text(">ref\nAC\n")
            (root/"calls.vcf").write_text("wrong\t1\t.\tA\tG\t.\t.\t.\tAD:DP\t0,1:1\n")
            with self.assertRaises(subprocess.CalledProcessError):
                generate(root/"calls.vcf", root/"ref.fa", root/"out.fa", "ANON", "iupac")

    def test_stale_and_missing_sidecars_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/"ref.fa").write_text(">ref\nAC\n")
            (root/"calls.vcf").write_text("ref\t1\t.\tA\tG\t.\t.\t.\tAD:DP\t0,1:1\n")
            output = root/"out.fa"
            generate(root/"calls.vcf", root/"ref.fa", output, "ANON", "iupac")
            output.write_text(">ANON\nAN\n")
            with self.assertRaisesRegex(ValueError, "checksum"):
                load(output)
            with self.assertRaisesRegex(ValueError, "Missing consensus sidecar"):
                load(root/"ref.fa")
            output.write_text(">ANON\nGN\n")
            sidecar = Path(str(output)+".consensus.json")
            data = json.loads(sidecar.read_text()); data["usable_depths"] = [1]
            sidecar.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "usable depths"):
                load(output)

    def test_qc_json_publishes_all_three_roles_and_compact_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root/"input.fa"
            reference.write_text(">ref\nACGT\n")
            vcf = root/"calls.vcf"
            vcf.write_text("ref\t1\t.\tA\tG\t.\t.\t.\tAD:DP\t0,6:6\n")
            scaffold = root/"ANON.fasta"
            generate(vcf, reference, scaffold, "ANON", "scaffold")
            vcf.write_text(vcf.read_text().replace("ref\t", "ANON\t"))
            generate(vcf, scaffold, root/"ANON-0.15-iupac.fasta", "ANON-0.15-iupac", "iupac")
            generate(vcf, scaffold, root/"ANON-majority.fasta", "ANON-majority", "analytical_majority")
            output = root/"ANON_qc_summary.json"
            subprocess.run([sys.executable, str(ROOT/"scripts/build_qc_summary.py"), "--results-dir", str(root), "--sample-id", "ANON", "--run-name", "R", "--output", str(output)], check=True, capture_output=True)
            qc = json.loads(output.read_text())
            self.assertIsNone(qc["pipeline_version"])
            self.assertEqual(set(qc["consensus"]), {"scaffold", "iupac", "analytical_majority"})
            self.assertEqual(qc["outputs"]["analytical_majority_fasta"], "ANON-majority.fasta")
            for key in ("main_consensus_metadata", "iupac_consensus_metadata", "analytical_majority_consensus_metadata"):
                self.assertTrue((root/qc["outputs"][key]).exists())
            self.assertNotIn("usable_depths", qc["consensus"]["iupac"])
            self.assertEqual(sequence(scaffold)[1], "ACGT")
            self.assertEqual(sequence(root/"ANON-majority.fasta")[1], "GNNN")

    def test_qc_json_records_explicit_pipeline_version(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root/"ANON_qc_summary.json"
            subprocess.run([sys.executable, str(ROOT/"scripts/build_qc_summary.py"),
                            "--results-dir", str(root), "--sample-id", "ANON",
                            "--run-name", "R", "--pipeline-version", "0.3.0",
                            "--output", str(output)], check=True, capture_output=True)
            self.assertEqual(json.loads(output.read_text())["pipeline_version"], "0.3.0")
