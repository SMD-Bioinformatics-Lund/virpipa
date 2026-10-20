"""Regression tests for unsafe consensus, mapping and rule evaluation."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import annotate_vcf_resistance as caller
import recall_resistance as recall
import consensus_metadata as consensus


def sidecar(path, mode="iupac", origin="a"*64, depths=None):
    identifier, seq = consensus.sequence(path)
    data = {"schema_version": "1.0", "converter_version": "4.0", "contig": identifier,
            "length": len(seq), "sequence_sha256": consensus.sequence_checksum(seq),
            "coordinate_source_sha256": origin, "mode": mode, "min_af": 0.15 if mode == "iupac" else None,
            "minimum_usable_depth": 1 if mode == "iupac" else 7, "reference_fill": mode == "scaffold",
            "usable_depths": depths or [10]*len(seq), "masking_reasons": {}}
    Path(str(path)+".consensus.json").write_text(json.dumps(data))
    return data


def rule(definition):
    return {"drug": "Drug", "region": "NS5A", "rule_definition": definition,
            "subtype_pattern": "1a", "drug_licensed_for_genotype": "1",
            "prediction": "resistant", "reference": "test"}


class SafetyTests(unittest.TestCase):
    def test_every_committed_rule_is_parsed_completely_including_or(self):
        for row in caller.rules(ROOT / "assets/hcv_geno2pheno_rules.csv"):
            self.assertTrue(caller.parts(row["rule_definition"]))
        sites = {("NS5A", 168): {"assessed": True, "amino_acids": ["H"]}}
        self.assertEqual(caller.evaluate([rule("168A or 168H or 168E")], sites, "1a")[0]["match_state"], "full")
        sites[("NS5A", 168)]["h77_aa"] = "D"
        self.assertEqual(caller.evaluate([rule("168any")], sites, "1a")[0]["match_state"], "full")
        sites[("NS5A", 168)]["amino_acids"] = ["D"]
        self.assertEqual(caller.evaluate([rule("168any")], sites, "1a")[0]["match_state"], "none")
        with self.assertRaises(ValueError):
            caller.parts("30R and unsupported")

    def test_unresolved_alternative_prevents_susceptible(self):
        sites = {("NS5A", 30): {"assessed": True, "amino_acids": ["Q"], "h77_aa": "Q"}}
        rows = caller.evaluate([rule("30R"), rule("93H")], sites, "1a")
        self.assertEqual(caller.outcome(rows, sites, True, rows), ("not assessed", None))

    def test_compound_disproved_even_with_missing_other_position(self):
        sites = {("NS5A", 30): {"assessed": True, "amino_acids": ["Q"], "h77_aa": "Q"}}
        rows = caller.evaluate([rule("30R and 93H")], sites, "1a")
        self.assertEqual(rows[0]["match_state"], "none")
        self.assertEqual(caller.outcome(rows, sites, True, rows), ("susceptible", None))

    def test_warning_only_for_sole_unresolved_compound(self):
        sites = {("NS5A", 30): {"assessed": True, "amino_acids": ["R"], "h77_aa": "Q"}}
        compound = caller.evaluate([rule("30R and 93H")], sites, "1a")
        self.assertTrue(caller.compound_warning(compound))
        self.assertFalse(caller.compound_warning(compound + caller.evaluate([rule("93H")], sites, "1a")))
        sites[("NS5A", 93)] = {"assessed": True, "amino_acids": ["Y"], "h77_aa": "Y"}
        self.assertFalse(caller.compound_warning(caller.evaluate([rule("30R and 93H")], sites, "1a")))

    def test_truncation_is_not_deletion(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"][:93]
        sites = caller.make_sites(sequence, {"NS5A": {"start": 1, "end": len(sequence), "strand": "+"}}, h77, [rule("32del")])
        self.assertFalse(sites[("NS5A", 32)]["assessed"])
        self.assertNotIn("del", sites[("NS5A", 32)]["amino_acids"])

    def test_frameshift_does_not_create_downstream_resistance(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        sequence = sequence[:60] + "A" + sequence[60:]
        sites = caller.make_sites(sequence, {"NS5A": {"start": 1, "end": len(sequence), "strand": "+"}}, h77, [rule("30H"), rule("20R")])
        self.assertFalse(sites[("NS5A", 30)]["assessed"])

    def test_read_codon_rejects_internal_insertion(self):
        import pysam
        read = pysam.AlignedSegment()
        read.query_sequence = "AGGCC"
        read.query_qualities = [40] * 5
        read.reference_start = 0
        read.cigartuples = [(0, 1), (1, 2), (0, 2)]
        self.assertIsNone(caller.read_codon(read, [0, 1, 2], 13))

    def test_n_is_unknown_even_when_its_translations_match_a_rule(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        sequence = sequence[:276] + "NNN" + sequence[279:]
        sites = caller.make_sites(sequence, {"NS5A": {"start": 1, "end": len(sequence), "strand": "+"}}, h77, [rule("93H")])
        self.assertFalse(sites[("NS5A", 93)]["assessed"])

    def test_coding_uncertainty_does_not_generate_compound_coverage_warning(self):
        sites = {("NS5A", 30): {"assessed": True, "amino_acids": ["R"]},
                 ("NS5A", 93): {"assessed": False, "assessment_reason": "coding_integrity"}}
        self.assertIsNone(caller.compound_warning(caller.evaluate([rule("30R and 93H")], sites, "1a")))

    def test_full_rule_wins_over_unresolved_alternative(self):
        sites = {("NS5A", 30): {"assessed": True, "amino_acids": ["R"]}}
        rows = caller.evaluate([rule("30R"), rule("93H")], sites, "1a")
        self.assertEqual(caller.outcome(rows, sites, True, rows), ("resistant", "resistant"))
        self.assertIsNone(caller.compound_warning(rows))

    def test_confirmed_internal_deletion_is_assessed_and_anchored(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        sequence = sequence[:93] + sequence[96:]
        site = caller.make_sites(sequence, {"NS5A": {"start": 1, "end": len(sequence), "strand": "+"}}, h77, [rule("32del")])[("NS5A", 32)]
        self.assertTrue(site["assessed"])
        self.assertEqual(site["amino_acids"], ["del"])
        self.assertEqual(site["genomic_start"], site["genomic_end"])

    def test_upstream_sites_remain_usable_after_downstream_frameshift(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        sequence = sequence[:600] + "A" + sequence[600:]
        sites = caller.make_sites(sequence, {"NS5A": {"start": 1, "end": len(sequence), "strand": "+"}}, h77, [rule("93H"), rule("250H")])
        self.assertTrue(sites[("NS5A", 93)]["assessed"])
        self.assertFalse(sites[("NS5A", 250)]["assessed"])

    def test_near_start_frameshifts_are_not_hidden_as_terminal_gaps(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        for position in (1, 2, 3):
            seq = sequence[:position]+"A"+sequence[position:]
            sites = caller.make_sites(seq, {"NS5A": {"start": 1, "end": len(seq), "strand": "+"}}, h77, [rule("93H")])
            self.assertFalse(sites[("NS5A", 93)]["assessed"], position)

    def test_ambiguous_deletion_placement_is_not_a_confirmed_deletion(self):
        reference = "ATG"+"GCT"*5+"TGG"
        h77 = {"NS5A": {"nuc": reference}}
        query = "ATG"+"GCT"*4+"TGG"
        site = caller.make_sites(query, {"NS5A": {"start": 1, "end": len(query), "strand": "+"}}, h77, [rule("4del")])[("NS5A", 4)]
        self.assertFalse(site["assessed"])

    def test_definite_stop_blocks_downstream_but_possible_stop_is_local(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        for codon, assessed in (("TAA", False), ("TRG", True)):
            seq = sequence[:60]+codon+sequence[63:]
            alignment = {}
            sites = caller.make_sites(seq, {"NS5A": {"start": 1, "end": len(seq), "strand": "+"}}, h77, [rule("93H")], alignments=alignment)
            self.assertEqual(sites[("NS5A", 93)]["assessed"], assessed)
            self.assertTrue(any("stop" in warning for warning in alignment["NS5A"]["warnings"]))


class ConsensusTests(unittest.TestCase):
    def convert(self, ad, alt="G", mode="iupac", depth=100, af="0.15"):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "ref.fa").write_text(">ref\nAC\n")
            (directory / "calls.vcf").write_text(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                f"ref\t1\t.\tA\t{alt}\t.\t.\t.\tAD:DP\t{ad}:{depth}\n")
            return subprocess.check_output(["awk", "-v", f"MODE={mode}", "-v", f"MIN_AF={af}", "-v", "MIN_DP=7", "-f", str(ROOT / "scripts/vcf_to_iupac.awk"), str(directory / "calls.vcf"), str(directory / "ref.fa")], text=True).splitlines()[1]

    def test_all_alt_excludes_reference(self):
        self.assertEqual(self.convert("0,100"), "GN")

    def test_multiallelic_and_exact_cutoff(self):
        self.assertEqual(self.convert("60,20,20", "G,T"), "DN")
        self.assertEqual(self.convert("85,15"), "RN")
        self.assertEqual(self.convert("86,14"), "AN")

    def test_missing_depth_masks_not_reference(self):
        self.assertEqual(self.convert("0,0", depth=0), "NN")

    def test_majority_scaffold_retains_uncovered_reference(self):
        self.assertEqual(self.convert("0,100", mode="majority"), "GC")
        self.assertEqual(self.convert("50,50", mode="majority"), "AC")

    def test_invariant_record_with_single_ad_is_observed(self):
        self.assertEqual(self.convert("100", alt="."), "AN")

    def test_indels_do_not_change_length_or_become_nucleotide_calls(self):
        self.assertEqual(self.convert("0,100", alt="AG"), "NN")
        self.assertEqual(self.convert("0,100", alt="AG", mode="majority"), "AC")

    def test_missing_ad_is_not_reference_evidence(self):
        self.assertEqual(self.convert("."), "NN")

    def test_analytical_consensus_inherits_sidecar_origin_and_headers_are_clean(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "ref.fa").write_text(">ref\nAC\n")
            sidecar(directory / "ref.fa", "scaffold")
            (directory / "calls.vcf").write_text("ref\t1\t.\tA\tG\t.\t.\t.\tAD:DP\t0,100:100\n")
            data = consensus.generate(directory/"calls.vcf", directory/"ref.fa", directory/"out.fa", "ANON", "iupac")
            self.assertEqual(data["coordinate_source_sha256"], "a"*64)
            self.assertEqual((directory/"out.fa").read_text().splitlines()[0], ">ANON")
            data = consensus.generate(directory/"calls.vcf", directory/"ref.fa", directory/"scaffold.fa", "ANON", "scaffold")
            self.assertEqual(data["coordinate_source_sha256"], consensus.sequence_checksum("AC"))


class RecallTransactionTests(unittest.TestCase):
    def test_failed_rollback_keeps_lock_for_explicit_recovery(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archives, backup = root / "archives", root / "backup"
            archives.mkdir()
            target, source = archives / "a", root / "new"
            target.write_text("old"); source.write_text("new")
            with patch.object(recall, "replace_file", side_effect=OSError("write unavailable")):
                with self.assertRaises(OSError):
                    recall.apply_transaction(archives, {target: source}, backup)
            self.assertTrue((archives / ".resistance-recall-lock").exists())
            recall.restore(backup, archives)
            self.assertEqual(target.read_text(), "old")
            self.assertFalse((archives / ".resistance-recall-lock").exists())

    def test_restored_journal_can_release_lock_after_interrupted_cleanup(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archives, backup = root / "archives", root / "backup"
            archives.mkdir(); backup.mkdir()
            lock = archives / ".resistance-recall-lock"
            lock.mkdir(); (lock / "backup.txt").write_text(str(backup))
            (backup / "manifest.json").write_text(json.dumps({"archives": str(archives), "state": "restored", "files": []}))
            recall.restore(backup, archives)
            self.assertFalse(lock.exists())

    def test_failed_replacement_restores_all_old_files_and_removes_new_ones(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archives, staged, backup = root / "archives", root / "staged", root / "backup"
            archives.mkdir(); staged.mkdir()
            for name in ("a", "b"):
                (archives / name).write_text("old-"+name)
                (staged / name).write_text("new-"+name)
            (staged / "c").write_text("new-c")
            replacements = {archives / name: staged / name for name in ("c", "a", "b")}
            original = recall.replace_file
            calls = 0
            def fail_once(source, destination):
                nonlocal calls
                calls += 1
                if calls == 3:
                    raise OSError("injected replacement failure")
                return original(source, destination)
            with patch.object(recall, "replace_file", side_effect=fail_once):
                with self.assertRaisesRegex(OSError, "injected"):
                    recall.apply_transaction(archives, replacements, backup)
            self.assertEqual((archives / "a").read_text(), "old-a")
            self.assertEqual((archives / "b").read_text(), "old-b")
            self.assertFalse((archives / "c").exists())
            self.assertFalse((archives / ".resistance-recall-lock").exists())
            self.assertEqual(json.loads((backup / "manifest.json").read_text())["state"], "restored")

    def test_completed_backup_is_non_overwriting_and_explicitly_restorable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archives, backup = root / "archives", root / "backup"
            archives.mkdir()
            destination, source = archives / "result", root / "new"
            destination.write_text("old"); source.write_text("new")
            recall.apply_transaction(archives, {destination: source}, backup)
            self.assertEqual(destination.read_text(), "new")
            with self.assertRaises(FileExistsError):
                recall.apply_transaction(archives, {destination: source}, backup)
            self.assertEqual((backup / "files/result").read_text(), "old")
            recall.restore(backup, archives)
            self.assertEqual(destination.read_text(), "old")

    def test_existing_aggregates_receive_fresh_qc_but_keep_other_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample, staged = root / "S", root / "staged"
            sample.mkdir(); staged.mkdir()
            path = sample / "S_qc_summary.json"
            old = {"sample_id": "S", "run_name": "R", "resistance": "old"}
            other = {"sample_id": "OTHER", "run_name": "R", "resistance": "untouched"}
            path.write_text(json.dumps(old))
            (root / "qc_summary.json").write_text(json.dumps([old, other]))
            (root / "qc_summary.jsonl").write_text("\n".join(map(json.dumps, [old, other])))
            fresh = dict(old, resistance="new")
            replacements = recall.stage_aggregates(root, {path: fresh}, staged)
            self.assertEqual(json.loads(replacements[root / "qc_summary.json"].read_text()), [fresh, other])
            self.assertEqual([json.loads(line) for line in replacements[root / "qc_summary.jsonl"].read_text().splitlines()], [fresh, other])


class AnnotationOutputTests(unittest.TestCase):
    def test_low_depth_qualifies_resistant_and_susceptible_calls_without_changing_outcome(self):
        reference = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")["NS5A"]["nuc"]
        for mutation in (False, True):
            sequence = reference[:276]+("CAT" if mutation else reference[276:279])+reference[279:]
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                depths = [10]*len(sequence)
                depths[276:279] = [1, 6, 7]
                low = self.run_caller(directory, sequence, depths=depths)
                high = self.run_caller(directory, sequence)
                self.assertEqual([d["outcome"] for d in low["drug_outcomes"]], [d["outcome"] for d in high["drug_outcomes"]])
                site = next(s for s in low["sites"] if s["gene"] == "NS5A" and s["h77_position"] == 93)
                self.assertTrue(site["low_depth"])
                self.assertEqual(site["nucleotide_depths"], [1, 6, 7])
                self.assertTrue(any("Low-depth support" in w for d in low["drug_outcomes"] for w in d["warnings"]))
                self.assertFalse(any("Low-depth support" in w for d in high["drug_outcomes"] for w in d["warnings"]))

    def test_four_allele_n_is_unassessed_not_missing_coverage_or_a_resistance_hit(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        sequence = sequence[:276]+"NAT"+sequence[279:]
        sites = caller.make_sites(sequence, {"NS5A": {"start": 1, "end": len(sequence), "strand": "+"}}, h77, [rule("93H")])
        caller.qualify_sites(sites, {"usable_depths": [4]*len(sequence), "length": len(sequence), "masking_reasons": {"277": "all_four_alleles"}})
        site = sites[("NS5A", 93)]
        self.assertFalse(site["assessed"])
        self.assertEqual(site["assessment_reason"], "all_four_alleles")
        evaluation = caller.evaluate([rule("93H")], sites, "1a")[0]
        self.assertEqual(evaluation["match_state"], "not_assessed")

    def test_non_resistant_mixture_allele_is_not_written_as_resistant(self):
        sequence = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")["NS5A"]["nuc"]
        sequence = sequence[:276]+"CAW"+sequence[279:]
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = self.run_caller(directory, sequence)
            site = next(s for s in payload["sites"] if s["gene"] == "NS5A" and s["h77_position"] == 93)
            self.assertEqual(site["amino_acids"], ["H", "Q"])
            import csv
            with (directory / "ANON_resistance.tsv").open() as handle:
                row = next(row for row in csv.DictReader(handle, delimiter="\t") if row["aa_pos"] == "93")
            self.assertEqual(row["alt_aa"], "H")
            self.assertIn("observed_aa=H%2CQ", (directory / "ANON_resistance_sites.gff3").read_text())
            self.assertIn("resistance_aa=H", (directory / "ANON_resistance_sites.gff3").read_text())

    def test_vadr_normalization_preserves_identifier_partial_flags_and_inference(self):
        raw = "ANON_WITH_UNDERSCORES\tVADR\tmat_peptide\t10\t101\t.\t+\t.\tID:mat_peptide1;product=NS3 protease/helicase;trunc3=1;\n"
        raw += "ANON_WITH_UNDERSCORES\tVADR\tmisc_feature\t200\t300\t.\t+\t.\tID:misc_feature1;note=similar to NS5A protein;\n"
        result = subprocess.check_output(["awk", "-f", str(ROOT / "scripts/normalize_vadr_gff.awk")], input=raw, text=True)
        self.assertIn("ANON_WITH_UNDERSCORES\tVADR\tgene", result)
        self.assertIn("trunc3=1;", result)
        self.assertIn("annotation_inferred=true;", result)
        self.assertIn("\t+\t.\t", result)

    def test_vadr_runtime_failure_cannot_emit_fake_bed(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            sequence = directory / "input.fasta"
            sequence.write_text(">ANON\nACGT\n")
            result = subprocess.run(["bash", str(ROOT / "scripts/vadr_annotate.sh"), str(sequence), str(directory), "ANON"], env=dict(os.environ, VADR_CONTAINER_RUNTIME="false"), capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((directory / "vadr/ANON.vadr.bed").exists())

    def run_caller(self, directory, sequence, track_reference=False, scaffold_origin="a"*64, depths=None, scaffold_sequence=None):
        sequence_path, annotation = directory / "sequence.fasta", directory / "genes.gff"
        sequence_path.write_text(">ANON-0.15-iupac\n"+sequence+"\n")
        sidecar(sequence_path, depths=depths)
        annotation.write_text(f"ANON\tfixture\tgene\t1\t{len(sequence)}\t.\t+\t.\tgene=NS5A;\n")
        command = [sys.executable, str(ROOT / "scripts/annotate_vcf_resistance.py"), "--fasta", str(sequence_path), "--gff", str(annotation), "--sample-name", "ANON", "--subtype", "1a", "--rules", str(ROOT / "assets/hcv_geno2pheno_rules.csv"), "--output-dir", str(directory)]
        if track_reference:
            scaffold = directory / "scaffold.fasta"
            scaffold.write_text(">ANON\n"+(scaffold_sequence or sequence)+"\n")
            sidecar(scaffold, "scaffold", origin=scaffold_origin)
            command += ["--track-reference-fasta", str(scaffold)]
        subprocess.run(command, check=True, capture_output=True, text=True)
        return json.loads((directory / "ANON_resistance.json").read_text())

    def test_low_depth_alt_can_differ_from_mapping_scaffold_without_changing_coordinates(self):
        scaffold = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")["NS5A"]["nuc"]
        analytical = scaffold[:276]+"CAT"+scaffold[279:]
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            payload = self.run_caller(directory, analytical, True, depths=[1]*len(analytical), scaffold_sequence=scaffold)
            self.assertEqual(payload["coordinate_reference"]["contig"], "ANON")
            self.assertTrue(any(d["outcome"] == "resistant" for d in payload["drug_outcomes"]))

    def test_track_rebasing_rejects_a_different_coordinate_origin(self):
        sequence = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")["NS5A"]["nuc"]
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(subprocess.CalledProcessError) as failure:
                self.run_caller(Path(temporary), sequence, True, scaffold_origin="b"*64)
            self.assertIn("shared coordinate-source checksum", failure.exception.stderr)

    def test_tracks_use_actual_reference_and_concise_names(self):
        h77 = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")
        sequence = h77["NS5A"]["nuc"]
        sequence = sequence[:276]+"CAT"+sequence[279:]
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for main in (False, True):
                payload = self.run_caller(directory, sequence, main)
                contig = "ANON" if main else "ANON-0.15-iupac"
                self.assertEqual(payload["coordinate_reference"]["contig"], contig)
                rows = [line.split("\t") for line in (directory / "ANON_resistance_sites.gff3").read_text().splitlines() if not line.startswith("#")]
                self.assertTrue(all(row[0] == contig for row in rows))
                self.assertTrue(any("Name=NS5A:Y93H" in row[8] for row in rows))
                self.assertFalse(any("Name=site:" in row[8] for row in rows))

    def test_truncated_and_deleted_genes_write_complete_outputs(self):
        sequence = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")["NS5A"]["nuc"]
        for seq, deletion in ((sequence[:93], False), (sequence[:93]+sequence[96:], True)):
            with self.subTest(deletion=deletion), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                payload = self.run_caller(directory, seq)
                self.assertTrue(all((directory / ("ANON"+suffix)).exists() for suffix in recall.SUFFIXES))
                site = next(s for s in payload["sites"] if s["gene"] == "NS5A" and s["h77_position"] == 32)
                self.assertEqual(site["assessed"], deletion)
                if deletion:
                    self.assertEqual(site["amino_acids"], ["del"])

    def test_cram_indel_disagreement_is_flagged_without_repair_and_wrong_reference_rejected(self):
        import pysam
        sequence = caller.h77_genes(ROOT / "refgenomes/resistance/AF011751.1.fasta")["NS5A"]["nuc"]
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, cram = directory / "ref.fasta", directory / "reads.cram"
            reference.write_text(">ANON\n"+sequence+"\n")
            pysam.faidx(str(reference))
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "ANON", "LN": len(sequence)}]}
            with pysam.AlignmentFile(str(cram), "wc", header=header, reference_filename=str(reference), format_options=[b"embed_ref=1"]) as handle:
                for index in range(10):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read-{index}"
                    read.query_sequence = sequence[:60]+"AC"+sequence[60:]
                    read.query_qualities = [40]*(len(sequence)+2)
                    read.flag = 0; read.reference_id = 0; read.reference_start = 0
                    read.mapping_quality = 60
                    read.cigartuples = [(0,60),(1,2),(0,len(sequence)-60)]
                    handle.write(read)
            pysam.index(str(cram))
            alignments = {"NS5A": {"genomic_start": 1, "genomic_end": len(sequence), "warnings": []}}
            caller.indel_evidence(str(cram), str(reference), alignments, 7, 13)
            events = alignments["NS5A"]["read_indel_evidence"]
            self.assertEqual(len(events), 1)
            self.assertEqual((events[0]["start"], events[0]["end"]), (61,90))
            self.assertEqual(events[0]["non_triplet_fraction"], 1)
            self.assertTrue(alignments["NS5A"]["warnings"])
            self.assertEqual(next(iter(caller.fasta(reference).values())), sequence)
            wrong = directory / "wrong.fasta"
            wrong.write_text(">ANON\n"+"A"+sequence[1:]+"\n")
            with self.assertRaisesRegex(ValueError, "checksum"):
                caller.evidence({}, str(cram), str(wrong), 13, 7)


if __name__ == "__main__":
    unittest.main()
