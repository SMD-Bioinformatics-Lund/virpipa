import importlib.util
import json
import shutil
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
    def test_geno2pheno_reference_preserves_scored_residues(self):
        reference = ROOT / "refgenomes/resistance/AF011751.1.fasta"
        records = RESISTANCE.fasta(reference)
        self.assertEqual(list(records), ["AF011751.1"])
        self.assertEqual(len(records["AF011751.1"]), 9599)
        genes = RESISTANCE.h77_genes(reference)
        previous = RESISTANCE.h77_genes(ROOT / "refgenomes/1a-AF009606.fa")
        self.assertEqual(genes["NS3"]["protein"][175], "G")
        for rule in RESISTANCE.rules(ROOT / "assets/hcv_geno2pheno_rules.csv"):
            for part in RESISTANCE.parts(rule["rule_definition"]):
                gene, position = rule["region"], part["position"] - 1
                self.assertEqual(
                    genes[gene]["protein"][position],
                    previous[gene]["protein"][position],
                    (gene, position + 1),
                )

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
                self.assertEqual(payload["h77_reference"], "AF011751.1")
                observed = {
                    f"{site['gene']}:{site['h77_position']}{aa}"
                    for site in payload["sites"]
                    for aa in site["amino_acids"]
                    if aa != site["h77_aa"]
                }
                self.assertTrue(set(expectation["mutations"]).issubset(observed), case)
                resistance_rows = (output / f"{case}_resistance.tsv").read_text()
                self.assertNotIn("\t156\tA\tT\t", resistance_rows, case)


class HistoricalRecallTests(unittest.TestCase):
    def test_recall_uses_blast_subtype_and_geno2pheno_reference(self):
        with tempfile.TemporaryDirectory() as temporary:
            archives = Path(temporary) / "archives"
            name = "G2P_CASE_008"
            sample = archives / name
            sample.mkdir(parents=True)
            fixture = ROOT / "assets/test_data/resistance/cases"
            shutil.copy2(fixture / f"{name}.fasta", sample / f"{name}-0.15-iupac.fasta")
            shutil.copy2(fixture / f"{name}.gff3", sample / f"{name}.vadr.pass_mod.gff")
            (sample / f"{name}-0.15-iupac.fasta.blast").write_text(
                "query acc.ver\tsubject acc.ver\n"
                f"{name}\t3a_D17763\t99\t1000\t0\t0\t1\t1000\t1\t1000\t0\t1000\n"
            )
            (sample / f"{name}-0.15-iupac.report.tsv").write_text(
                "reference\t1a-AF009606\n"
            )
            previous = sample / f"{name}_resistance.tsv"
            previous.write_text("previous result\n")
            backup = Path(temporary) / "backup"
            command = [
                sys.executable,
                str(ROOT / "scripts/recall_resistance.py"),
                str(archives),
                "--rules", str(ROOT / "assets/hcv_geno2pheno_rules.csv"),
            ]
            preview = subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertIn("subtype=3a", preview.stdout)
            self.assertEqual(previous.read_text(), "previous result\n")
            subprocess.run(
                command + ["--apply", "--backup-dir", str(backup)],
                check=True, capture_output=True, text=True,
            )
            payload = json.loads((sample / f"{name}_resistance.json").read_text())
            self.assertEqual(payload["subtype"], "3a")
            self.assertEqual(payload["h77_reference"], "AF011751.1")
            self.assertIn("\t93\tY\tH\t", previous.read_text())
            self.assertEqual(
                (backup / name / f"{name}_resistance.tsv").read_text(),
                "previous result\n",
            )


if __name__ == "__main__":
    unittest.main()
