import importlib.util
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location("rank_blast", Path(__file__).resolve().parents[1] / "scripts/rank_blast.py")
blast = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(blast)


def row(subject, qs, qe, ss, se, score, identity=100):
    length = abs(qe - qs) + 1
    return f"query\t{subject}\t{identity}\t{length}\t0\t0\t{qs}\t{qe}\t{ss}\t{se}\t0\t{score}\n"


class RankingTests(unittest.TestCase):
    def rank(self, rows, sequence="A" * 100):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "query.blast"
            path.write_text("query acc.ver\tsubject acc.ver\n" + rows)
            return blast.rank(path, sequence)

    def test_split_segments_beat_single_first_hit(self):
        result = self.rank(row("1a_A", 1, 100, 1, 100, 100)
                           + row("2b_B", 1, 40, 1, 40, 60)
                           + row("2b_B", 61, 100, 61, 100, 60), "A" * 40 + "N" * 20 + "A" * 40)
        self.assertEqual(result["subtype"], "2b")
        best = result["groups"][0]["references"][0]
        self.assertEqual(best["query_coverage_pct"], 80)
        self.assertEqual(best["observed_coverage_pct"], 100)
        self.assertEqual(best["bit_score"], 120)

    def test_overlaps_are_not_double_counted(self):
        result = self.rank(row("1a_A", 1, 70, 1, 70, 80)
                           + row("1a_A", 31, 100, 31, 100, 90))
        self.assertEqual(result["groups"][0]["bit_score"], 90)

    def test_subject_overlap_order_and_orientation(self):
        for start, end in [(10, 49), (1, 40), (100, 61)]:
            result = self.rank(row("1a_A", 1, 40, 41, 80, 50)
                               + row("1a_A", 61, 100, start, end, 60))
            self.assertEqual(result["groups"][0]["bit_score"], 60)

    def test_reverse_chain(self):
        result = self.rank(row("1a_A", 1, 40, 100, 61, 50)
                           + row("1a_A", 61, 100, 40, 1, 60))
        self.assertEqual(result["groups"][0]["bit_score"], 110)

    def test_no_cross_reference_aggregation(self):
        result = self.rank(row("1a_A", 1, 40, 1, 40, 60)
                           + row("1a_B", 61, 100, 61, 100, 60)
                           + row("2b_C", 1, 100, 1, 100, 100))
        self.assertEqual(result["subtype"], "2b")

    def test_tie_between_subtypes_is_ambiguous(self):
        result = self.rank(row("1a_A", 1, 100, 1, 100, 100)
                           + row("2b_B", 1, 100, 1, 100, 100))
        self.assertEqual(result["status"], "ambiguous")
        self.assertIsNone(result["subtype"])
        self.assertEqual(result["score_margin"], 0)

    def test_duplicates_and_row_order_do_not_change_result(self):
        rows = [row("1a_A", 1, 40, 1, 40, 60)] * 2 + [row("1a_B", 1, 40, 1, 40, 60)]
        self.assertEqual(self.rank("".join(rows)), self.rank("".join(reversed(rows))))
        self.assertEqual(self.rank("".join(rows))["groups"][0]["bit_score"], 60)

    def test_empty_and_unknown_top_do_not_assign(self):
        self.assertEqual(self.rank("")["status"], "unassigned")
        result = self.rank(row("unknown_A", 1, 100, 1, 100, 110)
                           + row("1a_A", 1, 100, 1, 100, 100))
        self.assertIsNone(result["subtype"])

    def test_genotype_only_is_preserved(self):
        self.assertEqual(self.rank(row("1_A", 1, 100, 1, 100, 100))["status"], "genotype_only")

    def test_invalid_rows_fail(self):
        for rows in ["broken\n", row("1a_A", 0, 99, 1, 100, 100),
                     row("1a_A", 1, 101, 1, 101, 100),
                     row("1a_A", 1, 100, 1, 100, "nan")]:
            with self.assertRaises(ValueError):
                self.rank(rows)

    def test_weighted_identity_and_explicit_nident(self):
        result = self.rank(row("1a_A", 1, 20, 1, 20, 30, 50)
                           + row("1a_A", 61, 100, 61, 100, 60, 100))
        self.assertAlmostEqual(result["groups"][0]["references"][0]["identity"], 100 * 50 / 60)
        text = row("1a_A", 1, 100, 1, 100, 100, 99.999).rstrip() + "\t100\n"
        self.assertEqual(self.rank(text)["groups"][0]["references"][0]["identical_bases"], 100)

    def test_provenance_and_single_record_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, search, database = [root / name for name in ("sample.fa", "sample.blast", "db.nsq")]
            fasta.write_text(">sample\n" + "A" * 100 + "\n")
            search.write_text(row("1a_A", 1, 100, 1, 100, 100))
            database.write_bytes(b"reference")
            result = blast.summarize(search, fasta, [database])
            self.assertEqual(result["input"]["role"], "analytical_15pct_iupac")
            self.assertEqual(len(result["database"]["db.nsq"]), 64)
            fasta.write_text(">one\nAAA\n>two\nAAA\n")
            with self.assertRaises(ValueError):
                blast.summarize(search, fasta)
