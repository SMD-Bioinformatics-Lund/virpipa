import hashlib
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import annotate_vcf_resistance as resistance


class ResistanceReferenceTests(unittest.TestCase):
    def test_panel_sequences_and_gene_coordinates_match_linking_alignments(self):
        directory = ROOT / "refgenomes/resistance"
        panel = json.loads((directory / "panel.json").read_text())
        alignments = {
            gene: resistance.fasta(directory / filename)
            for gene, filename in panel["region_alignments"].items()
        }
        subtypes = [item["subtype"] for item in panel["references"] if item["subtype"]]
        self.assertEqual(len(subtypes), len(set(subtypes)))
        for gene, records in alignments.items():
            self.assertEqual(len(set(map(len, records.values()))), 1, gene)
        for item in panel["references"]:
            records = resistance.fasta(directory / item["genome_fasta"])
            self.assertEqual(list(records), [item["accession"]])
            sequence = records[item["accession"]]
            self.assertEqual(len(sequence), item["length"])
            self.assertEqual(hashlib.sha256(sequence.encode()).hexdigest(), item["sequence_sha256"])
            for gene, feature in item["regions"].items():
                coding = sequence[feature["start"] - 1:feature["end"]]
                self.assertEqual(len(coding) % 3, 0, item["accession"])
                self.assertEqual(alignments[gene][item["accession"]].replace("-", ""), coding)
        incomplete = [
            (item["accession"], gene)
            for item in panel["references"] for gene, feature in item["regions"].items()
            if not feature["complete"]
        ]
        self.assertEqual(incomplete, [("DQ155561.1", "NS5B")])

    def test_subtype_reference_residues_at_scored_positions(self):
        alignment = resistance.fasta(ROOT / "refgenomes/resistance/alignments/NS5A.fasta")
        columns = [i for i, base in enumerate(alignment["AF011751.1"]) if base != "-"]
        for accession, position, amino_acid in (
            ("FJ462433.1", 31, "M"),
            ("FJ462439.1", 28, "I"),
            ("AJ851228.1", 93, "H"),
            ("DQ418786.1", 30, "R"),
            ("FJ462436.1", 30, "R"),
        ):
            codon = "".join(alignment[accession][columns[i]] for i in range((position - 1) * 3, position * 3))
            self.assertEqual(resistance.aas(codon), {amino_acid})

    def test_anonymous_cases_cover_each_additional_subtype_once(self):
        directory = ROOT / "assets/test_data/resistance/rare_cases"
        cases = json.loads((directory / "expected.json").read_text())["cases"]
        labels = []
        for alias, item in cases.items():
            records = resistance.fasta(directory / f"{alias}.fasta")
            self.assertEqual(list(records), [alias])
            self.assertRegex(alias, r"^G2P_RARE_\d{3}$")
            self.assertEqual(hashlib.sha256(records[alias].encode()).hexdigest(), item["sequence_sha256"])
            self.assertEqual(len(item["drugs"]), 16)
            self.assertEqual(len({region["service_subtype"] for region in item["regions"].values()}), 1)
            labels.append(item["regions"]["NS5A"]["service_subtype"])
            if labels[-1] == "1*2":
                self.assertEqual(item["regions"]["NS5A"]["genotype"], "1")
                self.assertIsNone(item["regions"]["NS5A"]["subtype"])
        self.assertCountEqual(labels, ["2c", "2i", "2k", "3h", "4c", "4d", "4k", "4m", "4r", "6e", "1*2"])


if __name__ == "__main__":
    unittest.main()
