"""Resolve real configuration and publish markers without running analyses."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
KEYS = ("outdir", "completion_log_dir", "virtitta_import_dir", "virtitta_run_dir")
PRODUCTION = ("/fs1/results/virpipa/hcv", "/fs1/results/cron/logs", "/fs1/results/cron/virtitta", "/access/virpipa/hcv")
DEV = ("/fs1/results_dev/virpipa/hcv", "/fs1/results_dev/cron/logs",
       "/fs1/results_dev/cron/virtitta", "/fs1/results_dev/virpipa/hcv")


@unittest.skipUnless(NEXTFLOW, "Nextflow must be on PATH for routing checks")
class DevRoutingTests(unittest.TestCase):
    def check_routing(self, arguments, expected, *, marker=False, parameters=None):
        with tempfile.TemporaryDirectory(prefix="virpipa-dev-") as temporary:
            directory = Path(temporary)
            shutil.copy(ROOT / "nextflow.config", directory)
            shutil.copytree(ROOT / "conf", directory / "conf")
            (directory / "qc.json").write_text("{}\n")
            (directory / "qc.jsonl").write_text("{}\n")
            (directory / "routing.nf").write_text(
                "include { WRITE_VIRTITTA_IMPORT_MARKER } from '"
                + str(ROOT / "modules/local/virtitta_import_marker/main.nf") + "'\n"
                "workflow {\n"
                "    println('ROUTING=' + groovy.json.JsonOutput.toJson(["
                + ", ".join(f"'{key}': params.{key}.toString()" for key in KEYS)
                + "]))\n"
                + ("    WRITE_VIRTITTA_IMPORT_MARKER(Channel.of(tuple('fixture_run', "
                   "file('qc.json'), file('qc.jsonl'))))\n" if marker else "")
                + "}\n"
            )
            arguments = [argument.replace("{tmp}", temporary) for argument in arguments]
            if parameters is not None:
                (directory / "params.json").write_text(json.dumps(parameters))
                arguments += ["-params-file", "params.json"]
            result = subprocess.run(
                [NEXTFLOW, "run", "routing.nf", "-ansi-log", "false", *arguments],
                cwd=directory, capture_output=True, text=True, timeout=120,
                env={**os.environ, "NXF_OFFLINE": "true", "NXF_SYNTAX_PARSER": "v2"},
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            resolved = json.loads(next(line.removeprefix("ROUTING=")
                                       for line in result.stdout.splitlines()
                                       if line.startswith("ROUTING=")))
            expected = tuple(value.replace("{tmp}", temporary) for value in expected)
            self.assertEqual(resolved, dict(zip(KEYS, expected)))
            if marker:
                markers = list(directory.rglob("*.sqlimport"))
                if not expected[2].strip():
                    self.assertEqual(markers, [])
                    self.assertFalse(list(directory.rglob(".command.sh")))
                else:
                    published = Path(expected[2]) / "fixture_run.sqlimport"
                    self.assertTrue(published.is_file())
                    self.assertEqual(published.read_text(),
                                     f"--run-dir {expected[3].rstrip('/')}/fixture_run/\n")

    def test_defaults(self):
        for arguments, expected in (([], PRODUCTION), (["--dev", "false"], PRODUCTION),
                                    (["--dev"], DEV),
                                    (["-profile", "tiny"], ("{tmp}/results_tiny", *PRODUCTION[1:])),
                                    (["-profile", "test"], ("{tmp}/results_test", *PRODUCTION[1:])),
                                    (["-profile", "slurm,hpc,apptainer"], PRODUCTION),
                                    (["-profile", "slurm,hpc,apptainer", "--dev"], DEV)):
            with self.subTest(arguments=arguments):
                self.check_routing(arguments, expected)

    def test_local_profiles_keep_local_destinations(self):
        for profile in ("local", "local_containers"):
            for flags in ([], ["--dev"]):
                for suffix, outdir in (("", "results"), (",tiny", "{tmp}/results_tiny"),
                                       (",test", "{tmp}/results_test")):
                    with self.subTest(profile=profile + suffix, flags=flags):
                        self.check_routing(["-profile", profile + suffix, *flags],
                                           (outdir, "", "", outdir), marker=True)

    def test_local_explicit_overrides_publish_markers(self):
        self.check_routing(
            ["-profile", "local", "--dev", "--outdir", "{tmp}/results",
             "--completion_log_dir", "{tmp}/logs", "--virtitta_import_dir", "{tmp}/markers"],
            ("{tmp}/results", "{tmp}/logs", "{tmp}/markers", "{tmp}/results"), marker=True)
        self.check_routing(
            ["-profile", "local_containers", "--virtitta_import_dir", "{tmp}/markers",
             "--virtitta_run_dir", "/custom/local"],
            ("results", "", "{tmp}/markers", "/custom/local"), marker=True)

    def test_profile_results_determine_dev_marker_target(self):
        for profile, outdir in (("tiny", "{tmp}/results_tiny"), ("test", "{tmp}/results_test")):
            with self.subTest(profile=profile):
                self.check_routing(["-profile", profile, "--dev"],
                                   (outdir, DEV[1], DEV[2], outdir))

    def test_custom_results_and_overrides_publish_real_markers(self):
        for target in ("{tmp}/results", "relative_results"):
            with self.subTest(target=target):
                self.check_routing(
                    ["--dev", "--outdir", target, "--completion_log_dir", "{tmp}/logs",
                     "--virtitta_import_dir", "{tmp}/markers"],
                    (target, "{tmp}/logs", "{tmp}/markers", target), marker=True)
        self.check_routing(
            ["--dev", "--outdir", "{tmp}/results", "--completion_log_dir=",
             "--virtitta_import_dir", "{tmp}/markers", "--virtitta_run_dir", "/custom/access///"],
            ("{tmp}/results", "", "{tmp}/markers", "/custom/access///"), marker=True)

    def test_disabled_markers_remain_disabled(self):
        self.check_routing(["--dev", "--virtitta_import_dir="],
                           (DEV[0], DEV[1], "", DEV[3]), marker=True)
        for disabled in ("", "   "):
            with self.subTest(disabled=disabled):
                self.check_routing(["--dev"], (DEV[0], DEV[1], disabled, DEV[3]),
                                   marker=True, parameters={"virtitta_import_dir": disabled})
