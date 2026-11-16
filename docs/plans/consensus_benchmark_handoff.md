# Maintenance / laptop-reboot handoff, 2026-11-13

The evaluation is complete as of the November 16 follow-up below; see
`consensus_benchmark_results.md`. No production converter or resistance behavior
has been changed. The branch is `consensus-benchmark`; baseline analysis code is
pinned separately from the experimental native-tool workflow. Historical recovery
instructions are retained below.

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

## Completed Hopper baseline (checked 2026-11-16)

- Job: **6405694**, COMPLETED, exit 0, elapsed 01:35:53.
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
The November 16 scheduler check showed no reservations and available normal nodes.
The laptop's `/tmp` was cleared; the durable checkpoint above survived.

All 240 published identities, layouts, consensus sidecars and FASTA/CRAM reference
identities passed validation; all 9,400 traced tasks completed successfully.
`scripts/validate_benchmark_baseline.py` records this comparison. Against the
November 11 archives, 195 analyses have at least one changed consensus sequence,
none has a changed subtype and one has changed drug outcomes (two drugs).
These differences precede the native-tool experiment and require separate
investigation. Native comparisons must use a fixed fresh baseline/read pool.
The private validation report is `baseline/validation_20261116.json`.

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
smoke test reached mapping on November 13, where the local license endpoint was
unavailable. Licensed baseline mapping worked on Hopper. No native Hopper jobs
were submitted before the maintenance warning. The later November 16 validation
and local license route are recorded below.

The first November 16 smoke submission (6428711) stopped before Nextflow because
Hopper's older Bash rejects an empty array expansion with `set -u`. The wrapper
now keeps its Nextflow argument array nonempty for both fresh and resumed runs.

iVar image verified locally and on Hopper:
`assets/containers/ivar_1.4.4.sif` and
`/fs1/resources/containers/ivar_1.4.4.sif`, SHA-256
`e6cab8972e28c0aa269f9f178716a5fe43dde70e8f3adb096e726ffc8ebce2f0`.
Never commit images to Git.

At the November 13 checkpoint, native pilots/cohort, additional truth controls,
real-cohort coverage preparation, distance comparison and final review remained.
They are completed in the November 16 follow-up below. No IQ-TREE topology
comparison or production deployment is part of this evaluation.

## November 16 follow-up

The 240-analysis native campaign completed at source checkpoint `43da431`.
All 2,160 native reference/alignment identities, frozen input hashes and image
manifests passed local validation. Local mapper controls at frozen source
`2681e6c` completed all 41 cases across nine arms: 369/369 known-codon and
template-retention checks pass. The pinned local/Hopper images match byte for byte.

The user supplied a local Sentieon route. Its license port/tunnel must be opened
for the current session before local mapping. Tiny fixtures and routine analysis
run locally; large mapping/full pipeline runs belong on Hopper. No further
remote analysis is needed for this checkpoint.

The supplied November 11 work directory confirmed absence of `--reorder` in all
240 archived Hostile scripts. An inspected sample retained the same filtered
pair multiset but replaced 522,139 of one million polishing pairs after ordering
changed. Treat old-versus-fresh sequence differences as confounded by this change.

Durable local outputs are under `test_output_consensus_benchmark/session_20261116/`.
`CHECKPOINT.json` records current paths. `cohort-scorecard-qualified/` preserves
upstream integrity reasons while verifying that the correction leaves all
amino-acid calls and assessment states unchanged. `depth-sensitivity.json`
separates one- versus seven-complete-codon gates. Original scorecards remain intact.

The local distance job completed in `cohort-distances-03/`, with two independent
single-thread MAFFT workers and reuse of verified completed alignments. It produced
426 matrices, 246,088 pair comparisons and 294 difference figures. Earlier
`cohort-distances/` and `cohort-distances-02/` retain interrupted sequential/partial
calculations; do not present them as final. The current script source was copied
into the current output directory. If interrupted, use a new output directory
and `--reuse-alignments` against the latest complete/partial output root; only
complete outputs with matching source FASTA bytes/settings/version are reused.
Keep the prepared manifest/masks and previous aligner.json files.

`distance-review.json` records all 446 changed neighbour records with before/after
counts and denominators; `distance-extremes.json` decomposes representative largest
changes. `distance-overview.png`/`.pdf` use pooled alignments to keep the cohort
constant despite incompatible native iVar majority interfaces. The completed
results and recommendation are in `consensus_benchmark_results.md`. No benchmark
compute jobs remain running. Production remains unchanged.
