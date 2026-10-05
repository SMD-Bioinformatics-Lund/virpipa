import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "resistance", ROOT / "scripts/annotate_vcf_resistance.py"
)
RESISTANCE = importlib.util.module_from_spec(SPEC)
sys.path.insert(0, str(ROOT / "scripts"))
SPEC.loader.exec_module(RESISTANCE)


class ResistanceUnitTests(unittest.TestCase):
    def test_iupac_codon_reports_every_possible_amino_acid(self):
        self.assertEqual(RESISTANCE.aas("CAY"), {"H"})
        self.assertEqual(RESISTANCE.aas("TAY"), {"Y"})
        self.assertEqual(RESISTANCE.aas("MAY"), {"H", "N"})

    def test_alignment_maps_positions_across_an_insertion_and_deletion(self):
        mapping = RESISTANCE.position_map("ABCDE", "ABXCDE")
        self.assertEqual([mapping[index] for index in range(1, 6)], [1, 2, 4, 5, 6])
        mapping = RESISTANCE.position_map("ABCDE", "ABDE")
        self.assertIsNone(mapping[3])

    def test_compound_rule_is_partial_until_all_parts_match(self):
        rule = {
            "drug": "Drug",
            "region": "NS5A",
            "rule_definition": "30R and 93H",
            "subtype_pattern": "1a",
            "drug_licensed_for_genotype": "1",
            "prediction": "resistant",
            "reference": "ref",
        }
        sites = {
            ("NS5A", 30): {"assessed": True, "amino_acids": ["R"], "h77_aa": "Q"},
            ("NS5A", 93): {"assessed": True, "amino_acids": ["Y"], "h77_aa": "Y"},
        }
        self.assertEqual(
            RESISTANCE.evaluate([rule], sites, "1a")[0]["match_state"], "partial"
        )
        sites[("NS5A", 93)]["amino_acids"] = ["H"]
        self.assertEqual(
            RESISTANCE.evaluate([rule], sites, "1a")[0]["match_state"], "full"
        )


class OnlineParityCases(unittest.TestCase):
    def test_anonymized_consensuses_recover_online_mutations(self):
        fixture_dir = ROOT / "assets/test_data/resistance/cases"
        expected = json.loads((fixture_dir / "expected.json").read_text())
        with tempfile.TemporaryDirectory() as temporary:
            for case, expectation in expected.items():
                output = Path(temporary) / case
                subprocess.run(
                    [
                        sys.executable,
                        str(ROOT / "scripts/annotate_vcf_resistance.py"),
                        "--fasta",
                        str(fixture_dir / f"{case}.fasta"),
                        "--gff",
                        str(fixture_dir / f"{case}.gff3"),
                        "--subtype",
                        expectation["subtype"],
                        "--rules",
                        str(ROOT / "assets/hcv_geno2pheno_rules.csv"),
                        "--h77-fasta",
                        str(ROOT / "refgenomes/1a-AF009606.fa"),
                        "--sample-name",
                        case,
                        "--output-dir",
                        str(output),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                payload = json.loads((output / f"{case}_resistance.json").read_text())
                observed = {
                    f"{site['gene']}:{site['h77_position']}{aa}"
                    for site in payload["sites"]
                    for aa in site["amino_acids"]
                    if aa != site["h77_aa"]
                }
                self.assertTrue(set(expectation["mutations"]).issubset(observed), case)
                resistance_rows = (output / f"{case}_resistance.tsv").read_text()
                self.assertNotIn("\t156\tA\tT\t", resistance_rows, case)


if __name__ == "__main__":
    unittest.main()
