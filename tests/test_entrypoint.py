"""Strict entry-point parsing and isolated completion-handler integration checks."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@unittest.skipUnless(NEXTFLOW, "Nextflow must be on PATH for entry-point checks")
class EntrypointTests(unittest.TestCase):
    def run_nextflow(self, directory, *arguments):
        return subprocess.run(
            [NEXTFLOW, *arguments], cwd=directory, capture_output=True, text=True,
            env={**os.environ, "NXF_OFFLINE": "true", "NXF_SYNTAX_PARSER": "v2"},
            timeout=120,
        )

    def test_real_entrypoints_pass_strict_lint_without_execution(self):
        with tempfile.TemporaryDirectory(prefix="virpipa-entrypoint-") as temporary:
            result = self.run_nextflow(temporary, "lint", str(ROOT / "main.nf"),
                                       str(ROOT / "test_module.nf"))
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_workflow_tuple_appends_preserve_order_and_nested_artifacts(self):
        source = (ROOT / "workflows/hcvpipe.nf").read_text()
        expressions = []
        for prefix in ("key, input, subtype, metadata, scaffold_metadata",
                       "key, input, artifacts"):
            match = re.search(re.escape(prefix) + r" -> (.*?) \}", source)
            self.assertIsNotNone(match)
            expressions.append(match.group(1))
        with tempfile.TemporaryDirectory(prefix="virpipa-tuples-") as temporary:
            script = Path(temporary) / "tuples.nf"
            script.write_text(
                "workflow {\n"
                "    def input = tuple('RUN', 'SAMPLE', 'GFF', 'FASTA', 'CRAM', 'CRAI', 'SCAFFOLD')\n"
                "    def metadata = 'IUPAC_METADATA'\n"
                "    def scaffold_metadata = 'SCAFFOLD_METADATA'\n"
                "    def resistance = " + expressions[0] + "\n"
                "    assert resistance == ['RUN', 'SAMPLE', 'GFF', 'FASTA', 'CRAM', 'CRAI', 'SCAFFOLD', 'IUPAC_METADATA', 'SCAFFOLD_METADATA']\n"
                "    assert resistance.getClass() == input.getClass()\n"
                "    assert input.size() == 7\n"
                "    input = tuple((0..<37).collect { index -> index })\n"
                "    def artifacts = [file('majority.fasta'), file('metadata.json')]\n"
                "    def final_input = " + expressions[1] + "\n"
                "    assert final_input.size() == 38\n"
                "    assert final_input[0..36] == input\n"
                "    assert final_input[37] == artifacts\n"
                "    assert final_input.getClass() == input.getClass()\n"
                "    assert input.size() == 37\n"
                "}\n"
            )
            result = self.run_nextflow(temporary, "run", str(script), "-ansi-log", "false")
            log = Path(temporary) / ".nextflow.log"
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr +
                             (log.read_text()[-6000:] if result.returncode and log.exists() else ""))

    def check_handler(self, *, fail=False, enabled=True, run_name="", column="run_name"):
        with tempfile.TemporaryDirectory(prefix="virpipa-completion-") as temporary:
            directory = Path(temporary)
            # Exercise the real entry point/handler, replacing only the analysis workflow.
            (directory / "main.nf").write_text((ROOT / "main.nf").read_text())
            (directory / "workflows").mkdir()
            (directory / "workflows/hcvpipe.nf").write_text(
                "process COMPLETION_TEST_FAILURE {\n"
                "    script:\n"
                "    '''\n"
                "    echo 'intentional completion-handler test failure' >&2\n"
                "    exit 13\n"
                "    '''\n"
                "}\n"
                "workflow HCVPIPE {\n"
                "    if (params.test_fail) {\n"
                "        COMPLETION_TEST_FAILURE()\n"
                "    }\n"
                "}\n"
            )
            (directory / "nextflow.config").write_text(
                "params {\n    run_name = ''\n    input = ''\n    csv = ''\n"
                "    completion_log_dir = ''\n    test_fail = false\n}\n"
                "process.executor = 'local'\n"
            )
            csv = directory / "samples.csv"
            csv.write_text(f"sample_id,{column}\nANON,csv_run\n")
            logs = directory / "completion"
            arguments = ["run", "main.nf", "-ansi-log", "false",
                         "--csv", str(csv)]
            if run_name:
                arguments += ["--run_name", run_name]
            if enabled:
                arguments += ["--completion_log_dir", str(logs)]
            if fail:
                arguments += ["--test_fail", "true"]
            result = self.run_nextflow(directory, *arguments)
            output = result.stdout + result.stderr
            self.assertEqual(result.returncode == 0, not fail, output)
            self.assertNotIn("Failed to invoke", output)
            self.assertEqual(len(list(directory.rglob(".command.sh"))), int(fail), output)
            if not enabled:
                self.assertFalse(logs.exists())
                self.assertFalse(list(directory.rglob("*-HCV.complete")))
                return
            expected = logs / f"{run_name or 'csv_run'}-HCV.complete"
            self.assertEqual(list(logs.iterdir()), [expected])
            text = expected.read_text()
            self.assertIn(f"Success     : {str(not fail).lower()}", text)
            self.assertIn(f"csv         : {csv}", text)
            if fail:
                self.assertIn("intentional completion-handler test failure", text)

    def test_success_uses_explicit_run_name(self):
        self.check_handler(run_name="explicit_run")

    def test_success_uses_samplesheet_run_name(self):
        self.check_handler()

    def test_failure_uses_samplesheet_sequencing_run(self):
        self.check_handler(fail=True, column="sequencing_run")

    def test_disabled_logging_creates_no_completion_marker(self):
        self.check_handler(enabled=False)
