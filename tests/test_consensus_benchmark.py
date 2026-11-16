"""Known-molecule acceptance tests independent of consensus output."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from benchmark_consensus import fixtures, codon_counts
from benchmark_distances import compare_cells, nearest, run_alignment
from benchmark_campaign import validate, fingerprint
from benchmark_native_score import qualify_read_sites
from annotate_vcf_resistance import read_codon
import pysam
from collections import Counter
from prepare_benchmark_distances import depth_mask
from summarize_benchmark_native import compare


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

    def test_linkage_truth_distinguishes_equal_marginal_mutation_frequencies(self):
        for name in ('linked_mutations', 'unlinked_mutations'):
            case = next(c for c in self.cases if c['case'] == name)
            observed = Counter()
            with pysam.AlignmentFile(case['bam']) as reads:
                for read in reads.fetch():
                    codons = [read_codon(read, positions, 13)
                              for positions in case['linked_codon_positions_zero_based']]
                    observed[':'.join(codons)] += 1
            self.assertEqual(dict(observed), case['haplotypes'])
            self.assertEqual(sum(v for k, v in observed.items() if k.startswith('GTT:')), 16)
            self.assertEqual(sum(v for k, v in observed.items() if k.endswith(':GAA')), 16)
            self.assertEqual(observed['GTT:GAA'], 16 if name == 'linked_mutations' else 0)

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

    def test_partial_or_different_alignment_inputs_are_not_reused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            previous = root/'previous'
            previous.mkdir()
            source = root/'one.input.fasta'
            source.write_text('>a\nACGT\n')
            (previous/source.name).write_text(source.read_text())
            cached = previous/'one.aligned.fasta'
            cached.write_text('>a\nACG\n')
            result = run_alignment(source, root/'first.fasta', 'unused-single-sequence', cached)
            self.assertIsNone(result['reused_from'])
            self.assertEqual(result['sequences'], ['ACGT'])
            cached.write_text('>a\nACGT\n')
            (previous/source.name).write_text('>a\nACGA\n')
            result = run_alignment(source, root/'second.fasta', 'unused-single-sequence', cached)
            self.assertIsNone(result['reused_from'])

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
                 2: dict(site, deletion_anchor=True),
                 3: dict(site, assessed=False, assessment_reason='internal_stop_codon',
                         codon_evidence={'codon_counts': {}})}
        qualify_read_sites(sites)
        self.assertEqual(sites[0]['amino_acids'], ['A', 'F'])
        self.assertFalse(sites[1]['assessed'])
        self.assertEqual(sites[2]['assessment_reason'], 'deletion_frequency_not_quantified')
        self.assertEqual(sites[3]['assessment_reason'], 'internal_stop_codon')

    def test_native_coverage_requires_complete_matching_reference_coordinates(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'depth.tsv'
            path.write_text('edited\t1\t3\nedited\t2\t0\nedited\t3\t1\n')
            self.assertEqual(depth_mask(path, 'edited', 3), [True, False, True])
            with self.assertRaisesRegex(ValueError, 'Incomplete'):
                depth_mask(path, 'edited', 4)
            with self.assertRaisesRegex(ValueError, 'Invalid'):
                depth_mask(path, 'scaffold', 3)
            path.write_text('edited\t1\t3\nedited\t3\t1\n')
            with self.assertRaisesRegex(ValueError, 'Invalid'):
                depth_mask(path, 'edited', 3)

    def test_native_report_separates_lost_assessment_and_changed_amino_acids(self):
        site = dict(gene='NS3', h77_position=155, assessed=True,
                    assessment_reason=None, amino_acids=['R'], read_amino_acid_counts={'R': 100})
        before = dict(sites=[site], drug_outcomes=[])
        lost = dict(sites=[dict(site, assessed=False, assessment_reason='insufficient_complete_codons')],
                    drug_outcomes=[])
        self.assertEqual(compare(before, lost)['lost_assessments'], 1)
        self.assertEqual(compare(before, lost)['changed_assessed_amino_acids'], 0)
        changed = dict(sites=[dict(site, amino_acids=['R', 'K'])], drug_outcomes=[])
        self.assertEqual(compare(before, changed)['changed_assessed_amino_acids'], 1)


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
