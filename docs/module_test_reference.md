# Module Test Reference

This file collects the detailed local test matrix, fixture inventory, and parity caveats that were previously embedded in `AGENTS.md`.

## Canonical local test commands

Before module tests, check the actual pipeline entry point with `nextflow lint main.nf test_module.nf`.
Lint does not execute workflows or write completion markers; do not add `-format`.
With Nextflow on PATH, `python -m unittest discover -s tests -p test_entrypoint.py -v`
also runs isolated success/failure/disabled completion-handler tests using an
analysis stub: success submits no tasks; failure runs one trivial local shell task.
These tests copy the real entry point, not a rewritten handler, and never run analysis.
Without Nextflow, that integration-test class explicitly skips; the lint gate is
still required before a Hopper run. Use production's Nextflow version for the Hopper gate.

Add `--publish_mode debug` to module tests when you need their intermediate files published under `--outdir`.
The resistance completion sequence and representative-data review gates are in
[the resistance validation checklist](resistance_validation.md).

```bash
nextflow run test_module.nf -profile local,tiny --module subsample --subsample_reads 25
nextflow run test_module.nf -profile local,tiny --module subsample_publish --publish_mode debug
nextflow run test_module.nf -profile local,tiny --module subsample_publish --publish_mode routine --outdir /tmp/virpipa_subsample_routine
nextflow run test_module.nf -profile local_containers,tiny --module hostile
nextflow run test_module.nf -profile local --module bam2fasta
nextflow run test_module.nf -profile local_containers --module bam2fasta
nextflow run test_module.nf -profile local --module bestref
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module mapping
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module mapping_noopt
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module polish
nextflow run test_module.nf -profile local --module consensus
nextflow run test_module.nf -profile local_containers --module consensus_contract --publish_mode debug
nextflow run test_module.nf -profile local_containers --module variantcall
nextflow run test_module.nf -profile local_containers --module filter_vcf
nextflow run test_module.nf -profile local_containers --module cram
nextflow run test_module.nf -profile local --module coverage
nextflow run test_module.nf -profile local_containers --module subtype
nextflow run test_module.nf -profile local --module report
nextflow run test_module.nf -profile local_containers --module vadr
nextflow run test_module.nf -profile local_containers --module resistance
nextflow run test_module.nf -profile local_containers --module finalize_results --outdir /tmp/virpipa_finalize_routine
```

## Tested modules and current parity status

- `hostile`: retained pair membership, filenames and QC interface are preserved;
  Nextflow now uses `--reorder` intentionally to stabilize downstream seeded
  subsampling. Compressed checksums and ordering need not match legacy outputs.
  Compare paired-record multisets and, between corrected repetitions, ordered payloads.
- `subsample`: read payloads and emitted filenames are unchanged. Debug publication
  separates requested limits under `fastq/subsample_<requested_pairs>/`.
- `subsample_publish`: runs both aliases on the same tiny sample with limits 25 and
  50. On Nextflow 26.04.3, check both published R1/R2 pairs under
  `tiny_run/SAMPLE001/fastq/subsample_25/` and `subsample_50/`, including record counts
  and mate pairing. Repeat with `--publish_mode routine` and a fresh output directory;
  neither subsampled set should be published. This exercises publication, not assembly.
- `bam2fasta`: consensus generation intentionally differs from the legacy converter: majority mode selects the most abundant nucleotide rather than REF+ALT ambiguity; uncovered mapping-scaffold bases remain reference. VCF/stats can still differ only in embedded metadata.
- `bestref`: verified identical selected FASTA for the fixture.
- `mapping`: SAM payload verified; stats differ only in embedded command-line/path metadata.
- `mapping_noopt`: SAM payload verified; stats differ only in embedded command-line/path metadata.
- `polish`: final FASTAs and BAM SAM payloads verified; `samtools stats` differ only in embedded command-line/path metadata.
- `consensus`: publishes masked analytical IUPAC and majority FASTAs, indexes and matching JSON sidecars. All qualifying alleles are retained at any nonzero usable depth; zero coverage and uninterpretable evidence are masked without changing coordinate length. Headers contain only identifiers.
- `variantcall`: VCF body verified; header differs only in embedded reference path metadata.
- `filter_vcf`: filtered VCF bodies verified; stats differ only in embedded filename/path metadata.
- `cram`: SAM payload verified; CRAM/header/index differ only in embedded reference path metadata.
- `coverage`: aggregate TSV verified identical; the additional BED records merged depth >=1 intervals.
- `subtype`: now requests all database subjects, appends exact identical-base counts,
  and runs `RANK_SUBTYPE`. Verify the selected call and metrics; raw 12-column byte
  parity with the historical search is no longer the acceptance criterion.
- `resistance_typing`: connected BLAST → ranking → resistance fixture. This uses the
  existing legacy comparison FASTA/CRAM; its consensus is explicitly unverified and
  is not a release biological-validation sample.
- `report`: verified identical.
- `vadr`: verified identical.
- `resistance`: debug publication must contain all six resistance files, 16 drug outcomes and codon evidence. The sparse standard fixture is not a complete negative control: insufficient coverage/mapping can yield not assessed. Anonymized cases separately cover positive calls including NS5A 93H.
- `finalize_results`: flat Virtitta publish layout is fixture-backed; validate with `scripts/check_publish_layout.py`. A tiny synthetic versioned scaffold/IUPAC/CRAM fixture is generated at runtime for review-bundle preparation; unrelated historical outputs remain publication fixtures, not biological controls.

## Fixture inventory

Repo-local fixtures live under `assets/test_data/`:

- `tiny/`
- `bam2fasta/`
- `bestref/`
- `mapping/`
- `mapping_noopt/`
- `polish/`
- `consensus/`
- `variantcall/`
- `filter_vcf/`
- `cram/`
- `coverage/`
- `subtype/`
- `report/`
- `vadr/`
- `resistance/`
- `qc_summary/`
- `finalize/`

## Important parity caveats

### General

- Some output formats embed paths, timestamps, `UR:` references, or command lines. Compare payloads, not just raw file hashes.
- The mounted test sample on the laptop is `SAMPLE001`; `TEST001` references are often really about `sample_name` rather than the FASTQ basename.

### BAM / CRAM / stats

- For BAM/CRAM outputs, compare `samtools view` payload or headers rather than raw file bytes.
- `samtools stats` often differs only in embedded invocation metadata.

### Mapping and polishing

- Use the bash-original subsampled FASTQs for mapping and polishing parity checks, not the original full FASTQs.
- `MAP_READS` emits both legacy outputs and `*_with_ref` outputs; use the metadata-preserving channels when downstream logic depends on the chosen reference.
- The `0.15-iupac` no-opt mapping branch is fixture-backed and should stay aligned with the bash `umimapnoopt()` behavior.

### Consensus / report / resistance

- `bam2fasta` must emit the replacement `SAMPLE001.fasta` with the sample-name header for downstream steps.
- The `0.15-iupac` report branch must key on `sample_id`, not `run_name`.
- `FILTER_VCF.out.stats` emits a list; downstream code must pick the `m0.15` member where required.
- The resistance branch uses the analytical 15% IUPAC FASTA, not the filtered VCF,
  and receives the ranked IUPAC typing result, including explicit uncertainty.

### VADR

- The bash pipeline runs VADR on `results/${id}.fasta` and publishes `${id}.vadr.pass_mod.gff` plus `${id}.vadr.bed`.
- Local VADR tests should use the replacement `SAMPLE001.fasta`, not `SAMPLE001-0.15-iupac.fasta`.

## Suggested verification approach

`consensus_contract` is a connected scaffold → analytical consensus → resistance fixture run, with the same sample ID in two different runs. It verifies `(run_name, sample_id)` sidecar joins and clean headers, publishes both analytical FASTAs and six resistance artifacts per run, and deliberately omits CRAM evidence to test sidecar-based baseline depth qualifications. It does not run assembly, polishing or VADR.

For a changed module:

1. Run the smallest local fixture that exercises the changed behavior.
2. Compare the payload that matters, not just raw bytes.
3. Check published filenames and channel metadata, not only process success.
4. Only move to a larger run after the fixture-backed parity check is understood.
