# Maintenance / laptop-reboot handoff, 2026-11-13

The evaluation is unfinished. No production converter or resistance behavior has
been changed. The branch is `consensus-benchmark`; its baseline analysis code is
pinned separately from the experimental native-tool workflow.

## Durable state

- Local artifacts formerly in `/tmp` were copied to
  `test_output_consensus_benchmark/session_20261113/` under this repository.
  This directory is ignored by Git and includes private cohort metadata. Do not
  force-add it. Its `CHECKPOINT.json` records the experiment state.
- The validated archived inventory covers 240 analyses from 40 runs. Its source
  hashes were verified after the complete codon scan; the archived results were
  unchanged. Main results: `consensus-benchmark-cohort-01/`.
- Completed 36-case containerized probes: `consensus-benchmark-sif-01/`.
  Pinned raw-versus-normalized comparisons: `consensus-benchmark-all-02/`.
- Completed synthetic distance comparisons and heatmaps:
  `consensus-pinned-probe-distances-plots/`.
- Historical command/provenance records retain their original `/tmp` paths.
  The files themselves are saved with the same directory structure. Portable
  manifests, where supplied, are separate from the historical provenance.

## Running Hopper baseline

- Job: **6405694**. It reached host filtering, mapping and assembly successfully.
- Frozen analysis source:
  `/fs1/jonas/virpipa/consensus_benchmark_36bcbda_20261113/source`
- Output/work/cache root:
  `/fs1/jonas/virpipa/consensus_benchmark_36bcbda_20261113/baseline`
- Log:
  `/fs1/jonas/virpipa/consensus_benchmark_36bcbda_20261113/baseline-6405694.log`
- Input: `/fs1/jonas/hcv/bjorn-rerun/240samples.csv`; all 240 identities match the
  archived cohort, IDs are unique and all 480 FASTQ paths passed preflight.
  `campaign.json` contains the full raw-input hashes.
- Completion sentinel: `baseline/BASELINE_COMPLETE`. It is created only after
  Nextflow exits successfully and verifies the frozen CSV.
- Import-marker and completion-log publication are explicitly disabled.
- Earlier job **6405693** stopped before analysis. Preserve it as failed-startup
  evidence; do not resume that attempt. Hopper Git is 1.8.3.1, lacks worktrees,
  and is unavailable on the compute node. Slurm also spools scripts. The current
  wrapper uses the submission directory and reads a detached `.git/HEAD` instead.

Slurm runs independently of the laptop. A laptop reboot may interrupt monitoring
or local scripts, but does not terminate the submitted baseline. A Hopper service
window may terminate the job; work/cache files must be preserved for recovery.
Do not submit more jobs until maintenance scheduling is known.

## Recovery after the service window

1. Check `squeue -j 6405694`, `sacct` and the log. Do not start a second coordinator
   while the original job is still active. Check for `BASELINE_COMPLETE` first.
2. Preserve `baseline/work` and `baseline/launch/.nextflow*`; these are required
   for Nextflow resume. Do not clean or move them.
3. Use the updated `runme_consensus_benchmark.sh` from a published checkpoint.
   Submit it while **cd'd into the original detached source above**. This keeps
   the analysis source, project paths, work directory and cache fixed:

   ```bash
   cd /fs1/jonas/virpipa/consensus_benchmark_36bcbda_20261113/source
   sbatch --job-name=virpipa-consensus-resume \
     --output=/fs1/jonas/virpipa/consensus_benchmark_36bcbda_20261113/resume-%j.log \
     /PATH/TO/UPDATED/CHECKPOINT/scripts/runme_consensus_benchmark.sh \
     /fs1/jonas/virpipa/consensus_benchmark_36bcbda_20261113/baseline/input.csv \
     /fs1/jonas/virpipa/consensus_benchmark_36bcbda_20261113/baseline baseline resume
   ```

   The updated wrapper is not in the original `36bcbda` checkout. Obtain it from
   the saved branch checkpoint; do not change the original analysis checkout.
   It checks the recorded source commit, frozen CSV and all recorded raw-input
   hashes, then runs `-resume` in the original launch directory. It does not overwrite the original launch
   wrapper or baseline metadata. The resume branch has syntax checks but has
   not yet been exercised on a real interrupted job.
4. Verify archives/layout, consensus sidecars and reference/alignment identities
   after recovery. Cached work can be reused only when retained and valid.
   Interrupted publication may need focused republishing from retained work;
   do not assume existing published files are complete or delete them blindly.
5. Compare software provenance if maintenance changed shared resources. The
   resume wrapper verifies raw-input content but does not automatically verify
   all shared containers/models again.

## Native-tool experiment status

`benchmark_native.nf`, `benchmark_native_sample.sh`, `benchmark_native_score.py`
and `benchmark_native_manifest.py` are an experimental prototype, not a release.
They reconstruct native FASTAs and remap a fixed eligible primary-read pool,
without rerunning UMI consensus. Both native callers use BAQ disabled and
anomalous pairs included; AWK controls on the same calls expose these differences
from production defaults. Coverage is measured in each native reference frame.

The scorecard reuses current numbering, read-codon extraction and catalogue
evaluation, holding typing fixed. It aggregates synonymous codons before 15%,
retains coding-integrity exclusions and requires seven complete codons for an
assessed read-based result. It leaves deletion frequencies explicitly
unquantified. Projected archived VADR gene boundaries are an evaluation adapter,
not fresh annotations or a validated replacement resistance caller.

Strict Nextflow syntax, shell syntax and unit checks passed. The local helper
smoke test reached mapping, where the laptop Sentieon image lacked a working
license endpoint. The **Hopper** Sentieon image has its own license configuration
and baseline mapping worked. No native Hopper jobs were submitted before the
maintenance warning. Run a tiny synthetic helper test there first, then the two
reported samples, before expanding the native cohort.

iVar image verified locally and on Hopper:
`assets/containers/ivar_1.4.4.sif` and
`/fs1/resources/containers/ivar_1.4.4.sif`, SHA-256
`e6cab8972e28c0aa269f9f178716a5fe43dde70e8f3adb096e726ffc8ebce2f0`.
Never commit images to Git.

Remaining: execute native pilots/cohort, complete missing truth controls,
prepare real-cohort trimming/coverage manifests, compare distance matrices,
review all changed resistance calls and deliver the recommendation. No IQ-TREE
topology comparison or production deployment is part of this checkpoint.
