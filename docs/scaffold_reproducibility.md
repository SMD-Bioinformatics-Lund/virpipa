# Scaffold reproducibility investigation

Investigation started 2026-11-06 after comparing the retained resistance validation
runs at 9364731 and 84e4283. No pipeline behaviour, output schema or importer policy
is changed by these diagnostic tools. Case identities and detailed comparisons are
stored in a private diagnostic directory outside Git; tests use synthetic reads.

## Retained-data findings

Two differing samples (case A and case B) and one unchanged control were inspected
through their original Nextflow task work directories. Published debug FASTQs alone
are unsafe for this comparison: reference-selection and polishing subsampling both
publish the same FASTQ filenames into the shared per-sample `fastq/` directory.
The retained work directories distinguish their actual outputs.

In all three cases, old/new Hostile outputs contain identical paired-record multisets
but differ in order. The comparison retains full headers, sequences, qualities and
pairing, and preserves duplicate multiplicity; gzip metadata is not compared.

| Case | Hostile pairs, each run | Assembly/reference subsample | Polishing subsample |
| --- | ---: | --- | --- |
| A | 3,457,629 | 250,000 pairs; 208,059 replaced | 1,000,000 pairs; 638,288 replaced |
| B | 1,323,459 | 250,000 pairs; 173,541 replaced | 1,000,000 pairs; 209,653 replaced |
| Control | 194,069 | all pairs retained, reordered | all pairs retained, reordered |

"Replaced" counts pairs unique to each subsample, not twice that number. Commands
and CPU/memory requests for filtering, subsampling, assembly and polishing match
between the two runs in cases A and B. Current image paths are known; absent
historical image checksums must not be described as proof of identical historical
container bytes.

The first observed divergence is read order after filtering. Subsampling uses
`seqtk sample -s100`: the fixed seed does not make selection independent of input
ordering. Different retained pairs enter SPAdes and polishing in cases A and B.
Their assemblies and hybrid genomes already differ before Pilon. The control's
assembly, hybrid genome and final consensus are unchanged despite reordered reads.

## Pilon repeatability and cycling are separate

The fixed-input Hostile experiment also completed: two executions against unchanged
raw FASTQs retained the same 3,457,629 paired records but produced different ordered
payload and pair-ID hashes. Input-file SHA-256 values verified unchanged before/after
the experiment. Thus output reordering is reproducible without changing the input
reads. The installed Hostile 1.1.0 code defaults `reorder` to false and supports
passing native Bowtie2 `--reorder` on Linux. That candidate correction has not yet
been tested or applied to the pipeline.

Three fixed-input Pilon repetitions were run for each differing case on Hopper,
using one frozen hybrid genome plus the retained round-one paired/unpaired BAMs
and indexes. No remapping occurred between repetitions. All three output sequences
and change files agree within each case; sequences also match the original round-one
output. Frozen-input checksums are verified after the runs. This establishes
repeatability for the two tested inputs/settings, not universal determinism.

The iterative pipeline remaps after changing the genome, so its inputs change on
each round. In the retained case B run, rounds 6/8/10 repeat one sequence and 7/9
repeat another, ending at the ten-round limit. The earlier run also revisits sequence
states. Case A converges within each run but starts with different hybrid genomes.

The current legacy-derived stopping rule checks equality of consecutive change
counts or exhaustion of ten rounds, not sequence equality or repeated-state cycles.
Identical counts alone do not establish convergence. No stopping policy is changed
as part of this investigation.

## Downstream effects observed

For the three investigated cases, resistance-site assessment states, codons, amino
acids, low-depth flags and drug mutation lists agree between runs. All drug outcomes
agree across the six-sample validation batch. This does not establish that all
future datasets are unaffected.

Other outputs are affected: case A changes VADR fail/pass classification; case B's
target-gene coordinates shift by two bases. Some scaffold lengths and analytical
sequences differ. Compare aligned sequences, not naive same-index bases, when
interpreting distance or coordinate effects. Generating-depth sidecars and local
CRAM windows show zero/low support at several differing sequence ends. Case A also
has internal differences near bases 1522–1543 at approximately 731–1816 usable
observations, with substantial local clipping/indel signals. The discrepancy cannot
be dismissed as low coverage alone. Comparison blocks are approximate sequence
windows, not assertions of unique biological indel placement. These mappings use
different references and selected reads, so their counts cannot establish which
alternative is correct. Read-support adjudication remains necessary before choosing
a polishing policy; unchanged resistance calls alone do not resolve the issue.

## Diagnostic tools

`scripts/compare_fastq_pairs.py` streams paired FASTQs and externally sorts SHA-256
pair fingerprints with a 32-MB sort buffer and one sorting thread. It reports ordered
content, pair-ID order and order-independent pair multiplicity, rejects broken
pairing/malformed records, and requires a new output directory. Run it with left/right
R1/R2 paths and `--outdir` outside the source/results trees. It does not output raw
read identifiers or sequences in its JSON.

The two `scripts/replay_*_hopper.sh` wrappers are isolated Slurm diagnostic jobs,
not pipeline launchers. Submit from a private directory with external Slurm logs:

```bash
sbatch --output "$DIAGNOSTIC_ROOT/pilon.%j.out" \
  --error "$DIAGNOSTIC_ROOT/pilon.%j.err" \
  /path/to/virpipa/scripts/replay_pilon_hopper.sh \
  /path/to/hybrid.fasta /path/to/paired-1.bam /path/to/unpaired-1.bam \
  "$DIAGNOSTIC_ROOT/new_pilon_replay"

sbatch --output "$DIAGNOSTIC_ROOT/hostile.%j.out" \
  --error "$DIAGNOSTIC_ROOT/hostile.%j.err" \
  /path/to/virpipa/scripts/replay_hostile_hopper.sh \
  /path/to/raw_R1.fastq.gz /path/to/raw_R2.fastq.gz \
  /path/to/virpipa/scripts/compare_fastq_pairs.py \
  "$DIAGNOSTIC_ROOT/new_hostile_replay"
```

Pilon freezes its small reference/BAM inputs once, runs three repetitions and checks
the frozen inputs afterwards. Hostile runs twice against the same raw FASTQs and
checks source-read hashes before/after; it then compares the outputs. Both record
current container hashes and wrapper copies. A `COMPLETE` marker is written only
after the diagnostic finishes successfully. Existing output directories are rejected.

## Remaining decisions

All three replay jobs completed successfully, and the private report/job records are
retained outside Git. The comparison utility passed six synthetic regression tests;
both replay wrappers passed Bash syntax checks and ShellCheck. Next, test the native
Hostile ordering option and propose a separate correction for order-sensitive read
selection and cycle detection/termination, with explicit acceptance criteria. Do not silently sort
production inputs, alter sampling policy, mask observed alleles or select an arbitrary
cycle member to force agreement. Explain or correct material sequence differences
before collection-wide reprocessing. Fixing colliding debug publication names is a
separate small repair; it does not explain differences in the actual workdir inputs.
