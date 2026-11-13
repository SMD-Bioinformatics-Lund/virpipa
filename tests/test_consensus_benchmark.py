"""Known-molecule acceptance tests independent of consensus output."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from benchmark_consensus import fixtures, codon_counts
from benchmark_distances import compare_cells, nearest
from benchmark_campaign import validate, fingerprint
from benchmark_native_score import qualify_read_sites


class BenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scratch = tempfile.TemporaryDirectory()
        cls.cases = fixtures(Path(cls.scratch.name))

    @classmethod
    def tearDownClass(cls):
        cls.scratch.cleanup()

    def test_defined_molecules_match_production_extraction(self):
        for case in self.cases:
            with self.subTest(case=case['case']):
                actual = codon_counts(case['bam'], case['reference'], [60, 61, 62])
                for key in ('codons', 'amino_acids', 'depth', 'qualifying_amino_acids'):
                    self.assertEqual(actual[key], case[key])

    def test_synonymous_codons_sum_before_threshold(self):
        case = next(c for c in self.cases if c['case'] == 'synonymous_sum')
        self.assertEqual(case['codons']['GCT'], 8)
        self.assertEqual(case['codons']['GCC'], 8)
        self.assertEqual(case['amino_acids']['A'], 16)
        self.assertIn('A', case['qualifying_amino_acids'])

    def test_missing_evidence_never_becomes_nearest_identity(self):
        cells = {('a', 'b'): dict(substitutions=0, indel_events=0, total_events=0, compared_bases=0),
                 ('a', 'c'): dict(substitutions=1, indel_events=0, total_events=1, compared_bases=100)}
        self.assertEqual(nearest(cells)['a'], ['c'])

    def test_distance_changes_include_denominator(self):
        old = {('a', 'b'): dict(substitutions=0, indel_events=0, total_events=0, compared_bases=90)}
        new = {('a', 'b'): dict(substitutions=1, indel_events=0, total_events=1, compared_bases=100)}
        result = compare_cells(old, new)['pairs'][0]
        self.assertEqual(result['delta']['compared_bases'], 10)
        self.assertEqual(result['delta']['substitutions'], 1)

    def test_neighbours_use_the_same_available_sample_pairs(self):
        cell = dict(substitutions=1, indel_events=0, total_events=1, compared_bases=100)
        old = {('a', 'b'): cell, ('a', 'c'): cell}
        new = {('a', 'b'): cell}
        result = compare_cells(old, new)
        self.assertEqual(result['nearest_neighbour_changes'], {})
        self.assertEqual(result['missing_pairs']['baseline_only'], 1)

    def test_read_score_does_not_invent_codons_or_override_integrity_exclusions(self):
        site = dict(amino_acids=['A', 'F', 'S', 'V'], assessed=True,
                    assessment_reason=None, genomic_start=61,
                    codon_evidence={'codon_counts': {'TTT': 84, 'GCT': 8, 'GCC': 8}})
        sites = {0: dict(site), 1: dict(site, assessed=False, assessment_reason='ambiguous_mapping'),
                 2: dict(site, deletion_anchor=True)}
        qualify_read_sites(sites)
        self.assertEqual(sites[0]['amino_acids'], ['A', 'F'])
        self.assertFalse(sites[1]['assessed'])
        self.assertEqual(sites[2]['assessment_reason'], 'deletion_frequency_not_quantified')


class CampaignTests(unittest.TestCase):
    def test_missing_fastq_and_repeated_sample_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sheet = root/'input.csv'
            read = root/'read.fastq'
            read.write_text('@a\nAC\n+\nII\n')
            header = 'clarity_sample_id,sequencing_run,read1,read2\n'
            sheet.write_text(header+f'ANON,RUN,{read},{read}\n')
            self.assertEqual(len(validate(sheet, 1)), 1)
            with self.assertRaises(ValueError):
                validate(sheet, 240)
            sheet.write_text(sheet.read_text()+f'ANON,OTHER,{read},{read}\n')
            with self.assertRaisesRegex(ValueError, 'unique'):
                validate(sheet, 2)
            sheet.write_text(header+f'ANON,RUN,{read},{root / "missing"}\n')
            with self.assertRaisesRegex(ValueError, 'Missing FASTQ'):
                validate(sheet, 1)

    def test_fingerprints_capture_payload_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'input'
            path.write_bytes(b'ACTG')
            first = fingerprint(path)
            path.write_bytes(b'ACCG')
            self.assertNotEqual(first['sha256'], fingerprint(path)['sha256'])


if __name__ == '__main__':
    unittest.main()
