import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("compare_fastq_pairs", Path(__file__).resolve().parents[1] / "scripts/compare_fastq_pairs.py")
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


class FastqComparisonTests(unittest.TestCase):
    def compare(self, left, right):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = []
            for label, pairs in (("old", left), ("new", right)):
                for mate in (1, 2):
                    path = root / f"{label}-{mate}.fq.gz"
                    with gzip.open(path, "wb") as handle:
                        for name, sequence, quality in pairs:
                            handle.write(f"@{name}/{mate}\n{sequence}\n+\n{quality}\n".encode())
                    paths.append(path)
            return comparison.compare(*paths, root / "comparison")

    def test_identical_payload_ignores_gzip_metadata(self):
        result = self.compare([("anon", "ACGT", "IIII")], [("anon", "ACGT", "IIII")])
        self.assertTrue(result["same_ordered_pairs"])
        self.assertTrue(result["same_pair_multiset"])

    def test_reordering_is_separate_from_changed_content(self):
        pairs = [("one", "ACGT", "IIII"), ("two", "AAAA", "IIII")]
        result = self.compare(pairs, list(reversed(pairs)))
        self.assertFalse(result["same_ordered_pairs"])
        self.assertFalse(result["same_pair_id_order"])
        self.assertTrue(result["same_pair_multiset"])

    def test_base_quality_and_header_changes_are_detected(self):
        for record in [("anon", "ACGA", "IIII"), ("anon", "ACGT", "II!I"), ("other", "ACGT", "IIII")]:
            result = self.compare([("anon", "ACGT", "IIII")], [record])
            self.assertFalse(result["same_pair_multiset"])
            self.assertEqual(result["pairs_only_in_left"], 1)
            self.assertEqual(result["pairs_only_in_right"], 1)

    def test_duplicate_multiplicity_is_preserved(self):
        pair = ("anon", "ACGT", "IIII")
        result = self.compare([pair, pair], [pair])
        self.assertEqual(result["pairs_only_in_left"], 1)
        self.assertEqual(result["pairs_only_in_right"], 0)

    def test_missing_mate_mismatched_id_and_truncation_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left, right = root / "left.fq", root / "right.fq"
            left.write_text("@one/1\nACGT\n+\nIIII\n")
            for text in ("", "@two/2\nACGT\n+\nIIII\n", "@one/2\nACGT\n+\n"):
                right.write_text(text)
                with self.assertRaises(ValueError):
                    comparison.fingerprint(left, right, root, "invalid")

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "keep.txt"
            marker.write_text("keep")
            with self.assertRaises(FileExistsError):
                comparison.compare("unused", "unused", "unused", "unused", root)
            self.assertEqual(marker.read_text(), "keep")
