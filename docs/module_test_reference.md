# Module Test Reference

This file collects the detailed local test matrix, fixture inventory, and parity caveats that were previously embedded in `AGENTS.md`.

## Canonical local test commands

Add `--publish_mode debug` to module tests when you need their intermediate files published under `--outdir`.

```bash
nextflow run test_module.nf -profile local,tiny --module subsample --subsample_reads 25
nextflow run test_module.nf -profile local_containers,tiny --module hostile
nextflow run test_module.nf -profile local --module bam2fasta
nextflow run test_module.nf -profile local_containers --module bam2fasta
nextflow run test_module.nf -profile local --module bestref
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module mapping
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module mapping_noopt
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module polish
nextflow run test_module.nf -profile local --module consensus
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

- `hostile`: verified identical.
- `subsample`: verified identical.
- `bam2fasta`: FASTA identical; VCF/stats differ only in embedded header path/date metadata.
- `bestref`: verified identical selected FASTA for the fixture.
- `mapping`: SAM payload verified; stats differ only in embedded command-line/path metadata.
- `mapping_noopt`: SAM payload verified; stats differ only in embedded command-line/path metadata.
- `polish`: final FASTAs and BAM SAM payloads verified; `samtools stats` differ only in embedded command-line/path metadata.
- `consensus`: verified identical.
- `variantcall`: VCF body verified; header differs only in embedded reference path metadata.
- `filter_vcf`: filtered VCF bodies verified; stats differ only in embedded filename/path metadata.
- `cram`: SAM payload verified; CRAM/header/index differ only in embedded reference path metadata.
- `coverage`: verified identical.
- `subtype`: verified identical.
- `report`: verified identical.
- `vadr`: verified identical.
- `resistance`: positive subtype-3a `NS5A 93H` fixture is covered.
- `finalize_results`: flat Virtitta publish layout is fixture-backed; validate with `scripts/check_publish_layout.py`.

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
- The resistance branch uses the `m0.15` filtered VCF and derives subtype from the first data hit in the `0.15-iupac` BLAST output.

### VADR

- The bash pipeline runs VADR on `results/${id}.fasta` and publishes `${id}.vadr.pass_mod.gff` plus `${id}.vadr.bed`.
- Local VADR tests should use the replacement `SAMPLE001.fasta`, not `SAMPLE001-0.15-iupac.fasta`.

## Suggested verification approach

For a changed module:

1. Run the smallest local fixture that exercises the changed behavior.
2. Compare the payload that matters, not just raw bytes.
3. Check published filenames and channel metadata, not only process success.
4. Only move to a larger run after the fixture-backed parity check is understood.
