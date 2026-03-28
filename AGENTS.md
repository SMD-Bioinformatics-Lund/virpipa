# AGENTS.md - VirPipa Development Guidelines

## Overview

VirPipa is an HCV (Hepatitis C Virus) probe-capture assembly pipeline built with Nextflow DSL2. It refactors a monolithic bash script (`scripts/hcvpipe.sh`) into modular Nextflow processes. The goal is full 1:1 output parity with the bash pipeline.

**Language:** Nextflow DSL2, Bash scripts  
**Testing:** Manual via `test_module.nf` (no nf-test framework)  
**HPC Environment:** Hopper (SLURM) with apptainer containers

---

## Environment and handling

there is a HPC (Hopper) where the pipeline will run but this agent runs on a local laptop
Never add singularity/apptainer containers to git
You have previously analysed the bash pipeline and added the analysis to docs/pipeline_workflow.md

### The laptop

Have too little RAM/CPU to run the pipeline using a regular isolate with full fastq files, but running test with either stubs or a super tiny data set should be fine
There is a mamba environmant on called 'skrotis' which you should use when running things locally. That contains among other things apptainer, nextflow, samtools, bcftools, mafft, snp-sites etc. For conveniance you probably would like to keep this loaded at all times 'mamba activate skrotis'
The local dir for container images is in assets/containers/
Data on /fs1 on the hpc can be reached on /mnt/fs1 here on this computer.
When you have made commits you must also push them to the bare repo on hopper with the command 'git push --set-upstream hopper' before running the script for starting the pipeline. IMPORTANT: The runme_test.sh script automatically pulls changes before running, so you never need to manually pull on HPC - just push and run the test script.

### The HPC

have no internet
have a volume /fs1
the dir for containers is locates in /fs1/resources/containers/
You can run the pipeline yourself with: ssh rs-fe1 'cd /fs1/jonas/src/virpipa && bash -l runme.sh resume'
if you explicitly need to run without the -resume flag use ssh rs-fe1 'cd /fs1/jonas/src/virpipa && bash -l runme.sh'
There are special environment setting, so you will not be able to run the pipeline in any other way that with the runme.sh script.

### The original bash pipeline

It is located in the scripts/ directory
The output from the original bash pipeline is at /fs1/jonas/hcv/results/test_run_bash_original/ on the HPC
You can use something different than datamash, this was just because I needed to do the sorting in bash
The important output corresponds to what is found in the results/ folder on the bash pipeline. A lot of the other files and folders in the bash pipeline are intermediary data, which in nextflow mostly ends up in the workdir for nextflow
the fasta files that are relevent are in the results/ dir and called .fasta and the -0.15-iupac.fasta

### The nextflow pipeline

Prefere to use bash in script blocks instead of grooby when possible
The approach is to implement all the functionality of the bash pipeline one step at at time. Verify the output against the known bash version after every new step/module/function. Do not proceed unless the output files are identical. Make sure to make a solid plan first on all steps taken in the bash version. Then at every step where you introduce a new module, check all parameters first before you start coding, implement the code/module and then verify the output, it should be identical to the bash version both in which files that are emitted and the content of them.
Use coding conventions that fit better with the nextflow paradigms, such as using the container directive and not copying the bash code blindly.
The bash code contains a lot of copying, moving and renaming of files, which is not in the style of nextflow, so the goal is to create the same files in the output, but in several cases the bash script is not the same way as you would do in nextflow.
Checksums might not always work to compare the files from the bash pipeline, since some file types do have paths to working directories embedded in them, so the files will not have the same md5

## Build / Lint / Test Commands

### Validate Pipeline

```bash
# Validate Nextflow syntax and schema
nextflow run . -validate

# Validate against schema
nextflow schema validate nextflow.config

# Check configuration without running
nextflow config -profile hpc,apptainer
```

### Run Tests

```bash
# Single module test (via HPC wrapper)
bash runme_test.sh <module_name>

# Direct nextflow test (on Hopper)
cd /fs1/jonas/src/virpipa
nextflow run test_module.nf -profile hpc,apptainer --module <module_name>

# Run full test profile locally
nextflow run . -profile test,standard,apptainer

# Run tiny local module tests on the laptop
nextflow run test_module.nf -profile local,tiny --module subsample --subsample_reads 25
nextflow run test_module.nf -profile local_containers,tiny --module hostile

# Run fixture-backed bam2fasta test locally
nextflow run test_module.nf -profile local --module bam2fasta
nextflow run test_module.nf -profile local_containers --module bam2fasta

# Run fixture-backed bestref test locally
nextflow run test_module.nf -profile local --module bestref

# Run fixture-backed one-reference sentieon mapping test locally
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module mapping

# Run fixture-backed one-reference sentieon mapping test without bwa tuning
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module mapping_noopt

# Run fixture-backed sentieon + pilon polishing loop locally
SENTIEON_LICENSE=localhost:8990 nextflow run test_module.nf -profile local_containers --module polish

# Run fixture-backed consensus test locally
nextflow run test_module.nf -profile local --module consensus

# Run fixture-backed variant calling and filtered VCF tests locally
nextflow run test_module.nf -profile local_containers --module variantcall
nextflow run test_module.nf -profile local_containers --module filter_vcf

# Run fixture-backed CRAM test locally
nextflow run test_module.nf -profile local_containers --module cram

# Run fixture-backed coverage test locally
nextflow run test_module.nf -profile local --module coverage

# Run fixture-backed subtype BLAST test locally
nextflow run test_module.nf -profile local_containers --module subtype

# Run fixture-backed report test locally
nextflow run test_module.nf -profile local --module report

# Run fixture-backed VADR test locally
nextflow run test_module.nf -profile local_containers --module vadr

# Run fixture-backed resistance annotation test locally
nextflow run test_module.nf -profile local_containers --module resistance
```

### Available Test Modules

Currently tested modules:
- `hostile` - Remove human reads with hostile (verified identical)
- `subsample` - Subsample FASTQ reads with seqtk (verified identical)
- `bam2fasta` - Build consensus FASTA/VCF from a fixture BAM (FASTA identical; VCF/stats differ only in header path/date metadata when compared to bash-original)
- `bestref` - Select `3a-D17763.fa` from fixture mapping stats (verified identical chosen FASTA in workdir)
- `mapping` - Build bash-style `r11b2L25.bwa.umi.sort/filter.sort` BAMs for `SAMPLE001` vs `3a-D17763` (SAM payload verified; stats differ only in embedded command-line path metadata)
- `mapping_noopt` - Build bash-style `r11b2L25.bwa.umi.sort/filter.sort` BAMs for `SAMPLE001` vs `0.15-iupac` without the tuned bwa arguments (SAM payload verified; stats differ only in embedded command-line path metadata)
- `polish` - Run the pilon polishing loop from fixture reads + hybrid FASTA (final FASTAs and BAM SAM payloads verified; `samtools stats` differ only in embedded command-line path metadata)
- `consensus` - Build `0.15-iupac` FASTA from fixture VCF + replacement FASTA (verified identical)
- `variantcall` - Build `SAMPLE001-pilon.vcf.gz` from fixture BAM + replacement FASTA (VCF body verified; header differs only in embedded reference path metadata)
- `filter_vcf` - Build `SAMPLE001-pilon-m*.vcf.gz` from fixture pilon VCF (all filtered VCF bodies verified; stats differ only in embedded filename/path metadata)
- `cram` - Build `SAMPLE001.cram` from fixture BAM + replacement FASTA (SAM payload verified; CRAM/header/index differ only in embedded reference path metadata)
- `coverage` - Build `SAMPLE001-coverage.tsv` from fixture CRAM (verified identical)
- `subtype` - Build `SAMPLE001-0.15-iupac.fasta.blast` from fixture consensus FASTA (verified identical)
- `report` - Build `SAMPLE001-0.15-iupac.report.tsv` and `.fastanucfreq.tsv` from fixture stats/CRAM/FASTA (verified identical)
- `vadr` - Build `SAMPLE001.vadr.pass_mod.gff` and `SAMPLE001.vadr.bed` from fixture `SAMPLE001.fasta` (verified identical)
- `resistance` - Build `SAMPLE001_resistance.tsv`, `.bed`, and `_by_drug.tsv` from a synthetic subtype-3a `NS5A 93H` fixture VCF

### Current Laptop Notes

- The mounted laptop copy of `/fs1` is available under `/mnt/fs1`.
- A faster local copy of the bash-original comparison run is available under `~/hcv/SAMPLE001-nextflow-nfcore-scaffold/`; prefer `~/hcv/SAMPLE001-nextflow-nfcore-scaffold/results/` over `/mnt/fs1/...` when comparing the final `results/` tree on the laptop.
- A local copy of the full test FASTQs is being added under `~/hcv/test_data/`; prefer that over `/mnt/fs1/...` when you need the full-sample inputs on the laptop.
- A faster local hostile cache is available under `~/resources/hostile`; prefer that over `/mnt/fs1/resources/ref/micro/hostile` for laptop runs.
- The repo-local tiny fixture is in `assets/test_data/tiny/` with samplesheet `assets/tiny_samplesheet.csv`.
- The repo-local `bam2fasta` fixture is in `assets/test_data/bam2fasta/` and uses the `SAMPLE001-nextflow-nfcore-scaffold` bash-original outputs so sample naming matches the main test data.
- The repo-local `bestref` fixture is in `assets/test_data/bestref/` and contains only the initial reference-mapping `*.stats` files plus the expected `3a-D17763.fa` reference.
- The repo-local `mapping` fixture is in `assets/test_data/mapping/` and contains the expected bash-original `SAMPLE001-3a-D17763.r11b2L25.bwa.umi.sort/filter.sort` BAMs plus the reference FASTA; the actual test inputs are the bash-original subsampled FASTQs under `/mnt/fs1/jonas/hcv/results/test_run_bash_original/.../fastq/*.sub.fastq.gz`.
- The repo-local `consensus` fixture is in `assets/test_data/consensus/` and uses `SAMPLE001-pilon.vcf.gz` plus the replacement `SAMPLE001.fasta` reference from the same bash-original run.
- The repo-local `polish` fixture is in `assets/test_data/polish/` and contains the bash-original `SAMPLE001` subsampled read pair, the `SAMPLE001.hybrid.fasta` polishing input, the expected final `SAMPLE001`/`SAMPLE001-iupac` FASTAs, the expected final `SAMPLE001-pilon.r11b2L25.bwa.umi.sort/filter.sort` BAMs, and early-round spot-check artifacts (`SAMPLE001-pilon-1.*`, `SAMPLE001-pilon-3.changes`).
- The repo-local `mapping_noopt` fixture is in `assets/test_data/mapping_noopt/` and contains the bash-original `SAMPLE001-0.15-iupac.fasta` plus the expected bash-original no-opt `SAMPLE001-0.15-iupac.r11b2L25.bwa.umi.sort/filter.sort` BAMs.
- The repo-local `variantcall` fixture is in `assets/test_data/variantcall/` and uses the `SAMPLE001` pilon BAM plus replacement `SAMPLE001.fasta`.
- The repo-local `filter_vcf` fixture is in `assets/test_data/filter_vcf/` and contains the bash-original `SAMPLE001-pilon.vcf.gz` plus the expected `SAMPLE001-pilon-m*.vcf.gz` outputs.
- The repo-local `cram` fixture is in `assets/test_data/cram/` and uses the `SAMPLE001` pilon BAM plus replacement `SAMPLE001.fasta`, with expected `SAMPLE001.cram` outputs from the bash-original run.
- The repo-local `coverage` fixture is in `assets/test_data/coverage/` and uses the bash-original `SAMPLE001.cram` plus replacement `SAMPLE001.fasta`, with expected `SAMPLE001-coverage.tsv`.
- The repo-local `subtype` fixture is in `assets/test_data/subtype/` and uses `SAMPLE001-0.15-iupac.fasta` plus the expected bash-original `SAMPLE001-0.15-iupac.fasta.blast`.
- The repo-local `report` fixture is in `assets/test_data/report/` and uses the bash-original `SAMPLE001-0.15-iupac.vcf.gz.stats`, `SAMPLE001-0.15-iupac.cram`, `SAMPLE001-0.15-iupac.fasta`, plus the expected report and nucleotide-frequency outputs.
- The repo-local `vadr` fixture is in `assets/test_data/vadr/` and uses `SAMPLE001.fasta` plus the expected bash-original `SAMPLE001.vadr.pass_mod.gff` and `SAMPLE001.vadr.bed`.
- The repo-local `resistance` fixture is in `assets/test_data/resistance/` and uses a synthetic `SAMPLE001-positive-93H.vcf` together with the VADR/report fixtures to exercise a positive subtype-3a resistance call.
- The `tiny` profile points `params.input` at that tiny samplesheet and writes to `results_tiny` unless `--outdir` overrides it.
- Use `-profile local,tiny` for non-container local tests that rely on tools from the `skrotis` environment.
- Use `-profile local_containers,tiny` for Apptainer-backed local tests; inside Codex this may still require running the command outside the sandbox.
- Use `-profile local --module bam2fasta` for fast local logic checks, and `-profile local_containers --module bam2fasta` when you want the pinned bcftools/samtools versions from the container.
- Use `-profile local --module bestref` for the reference-selection fixture. This module is internal and does not publish to `params.outdir`, so inspect the workdir output or the process log when verifying it manually.
- Use `SENTIEON_LICENSE=localhost:8990 -profile local_containers --module mapping` for the sentieon mapping fixture. The license works locally on this laptop, but the variable must be present in the shell that launches Nextflow.
- Use `SENTIEON_LICENSE=localhost:8990 -profile local_containers --module mapping_noopt` for the no-opt sentieon mapping fixture. This matches the bash `umimapnoopt()` behavior used for the `0.15-iupac` branch.
- Use `SENTIEON_LICENSE=localhost:8990 -profile local_containers --module polish` for the sentieon + pilon fixture. This module is heavy enough that it may need a few minutes locally, but the repo-local fixture is small enough to be practical on the laptop.
- Use `-profile local_containers --module variantcall` and `-profile local_containers --module filter_vcf` for parity checks because they are bcftools-version sensitive.
- Use `-profile local_containers --module cram` for parity checks because CRAM headers embed reference paths and the module relies on samtools behavior matching the pipeline container.
- Use `-profile local --module coverage` for the coverage fixture because `samtools` is available in `skrotis`.
- Use `-profile local --module report` for the report fixture because the module only needs `samtools coverage` plus text processing, and the local `skrotis` toolchain reproduces the bash-original output exactly.
- Use `-profile local_containers --module subtype` for parity checks because `blastn` is not installed in `skrotis`, but the pinned `blast_2.16.0.sif` container reproduces the bash-original output exactly.
- Use `-profile local_containers --module vadr` for parity checks. On this laptop, point `params.vadr_model_dir` at `/home/jonas/resources/vadr/vadr-models-flavi`, not `/mnt/fs1/...`, because VADR startup on the WSL-mounted path is extremely slow.
- Use `-profile local_containers --module resistance` for the resistance fixture. The module accepts either the committed `assets/hbv_result_rules.csv` table or a normalized JSON rules artifact from `scripts/update_geno2pheno_rules.py`.
- Use `~/resources/hostile` as the default laptop hostile cache path. `/mnt/fs1/resources/ref/micro/hostile` remains a fallback if the local copy is unavailable.
- Laptop profiles currently downscale `REMOVE_HOSTILE` to 4 CPUs / 10 GB / 1h, `SUBSAMPLE_READS` to 2 CPUs / 2 GB / 30m, `BAM2FASTA` to 2 CPUs / 4 GB / 30m, `ASSEMBLE_SPADES` to 8 CPUs / 15 GB / 8h, `ASSEMBLE_HYBRID` to 8 CPUs / 15 GB / 4h, `MAP_READS` to 4 CPUs / 8 GB / 1h, `POLISH_PILON_LOOP` to 4 CPUs / 10 GB / 2h, `CREATE_CONSENSUS` to 2 CPUs / 2 GB / 30m, `VARIANT_CALLING` to 2 CPUs / 4 GB / 30m, `FILTER_VCF` to 2 CPUs / 2 GB / 30m, `CREATE_CRAM` to 2 CPUs / 2 GB / 30m, `SUBTYPE_BLAST` to 2 CPUs / 2 GB / 30m, `CREATE_REPORT` to 2 CPUs / 2 GB / 30m, and `ANNOTATE_VADR` to 2 CPUs / 8 GB / 1h so they fit local resources.
- The mounted test sample currently available on the laptop is `SAMPLE001`; if someone refers to `TEST001` they likely mean the sample_name column rather than the FASTQ basename.
- For `bam2fasta`, the generated FASTA is byte-identical to the bash-original fixture. The generated VCF and stats still differ in embedded path/date header metadata, which is expected and acceptable for parity checks.
- For `bestref`, only feed the initial reference-mapping stats into the module test. Do not include later `pilon*` or `0.15-iupac` stats when checking the subtype-selection behavior from the early mapping stage.
- For `mapping`, test against the bash-original subsampled FASTQs, not the original full FASTQs from `/mnt/fs1/jonas/hcv/test_data/`. With the correct `.sub.fastq.gz` inputs, both the raw and filtered BAM SAM payloads match the bash fixture; the remaining `samtools stats` diff is only the embedded command line/path header.
- `MAP_READS` now also emits metadata-preserving channels carrying `ref_name` (`raw_bams_with_ref`, `bams_with_ref`, `stats_with_ref`). Keep the older emits for compatibility, but use the `*_with_ref` channels when the workflow needs to recover the chosen best-reference BAM/CRAM/report path.
- For `mapping_noopt`, use the same repo-local subsampled read pair from `assets/test_data/polish/` together with `assets/test_data/mapping_noopt/SAMPLE001-0.15-iupac.fasta`. The raw and filtered BAM SAM payloads match the bash fixture; the remaining `samtools stats` diff is only the embedded command line/path header.
- For `polish`, final-convergence parity is now verified locally from the repo-local fixture: `SAMPLE001.fasta`, `SAMPLE001-iupac.fasta`, the final raw and filtered pilon BAM SAM payloads, `SAMPLE001-pilon-1.fasta`, `SAMPLE001-pilon-1.changes`, and `SAMPLE001-pilon-3.changes` all match the bash-original fixture. The remaining `*.bam.stats` diff is only the embedded `samtools stats` command-line/path header.
- `bam2fasta` also needs to emit a replacement `SAMPLE001.fasta` with the sample-name header; downstream `variantcall` and `consensus` should use that replacement FASTA, not the published `SAMPLE001-1.0-iupac.fasta`, because the BAM/VCF contig names are `SAMPLE001`.
- After `bam2fasta`, downstream CRAM/coverage steps should use the regenerated replacement `SAMPLE001.fasta` (`1.0-iupac` majority-call content with the `SAMPLE001` header), because that is what the bash pipeline places in `pilon/${id}.fasta` before creating `results/${id}.cram`.
- For `variantcall` and `filter_vcf`, compare VCF content after stripping `##bcftools...` command/version lines or other embedded path metadata. The variant bodies match; remaining diffs are in generated header metadata.
- For `cram`, compare `samtools view` payload or header content rather than raw file bytes. The CRAM binary and CRAI differ because the embedded `UR:` reference path changes with the work directory, but the alignment payload is the same.
- For `subtype`, the bash pipeline writes a literal header row before `blastn -outfmt 6`. That header contains a typo (`s. star ... t`) from the original bash line continuation, and the Nextflow module intentionally reproduces it so the `.blast` file is byte-identical to the bash-original result.
- For `report`, derive the output prefix from the `*.vcf.gz.stats` filename, not from the FASTA name. The bash report appends `samtools coverage` output transposed into key/value rows and omits the `rname` row, so the Nextflow module should do the same.
- `FILTER_VCF.out.stats` emits a list of seven stats files, not a single `Path`, so the workflow must pick the `m0.15` member before building the `0.15-iupac` report.
- The resistance branch should also pick the `m0.15` filtered VCF, not the unfiltered pilon VCF, and derive subtype from the first data hit in `SAMPLE001-0.15-iupac.fasta.blast`.
- The `0.15-iupac` report branch must key on `sample_id`, not `run_name`; otherwise the join never emits and the published `SAMPLE001-0.15-iupac.report.tsv`/`fastanucfreq.tsv` files are missing.
- `CREATE_REPORT` now takes an explicit `report_id` so the workflow can emit bash-style names like `SAMPLE001-0.15-iupac.report.tsv` even when the stats input comes from `SAMPLE001-pilon-m0.15.vcf.gz.stats`.
- For `annotate_vadr`, the bash pipeline runs VADR on `results/${id}.fasta` and publishes `SAMPLE001.vadr.pass_mod.gff` plus `SAMPLE001.vadr.bed`. The laptop test should use the replacement `SAMPLE001.fasta` fixture, not `SAMPLE001-0.15-iupac.fasta`.
- `scripts/update_geno2pheno_rules.py` is the single offline refresh entrypoint for geno2pheno rules. It can download the table or rebuild normalized output from an existing CSV snapshot.
- `POLISH_PILON_LOOP` is now fixture-backed and locally verified against the repo-local `assets/test_data/polish/` set with final-convergence parity plus early-round spot checks.
- The workflow now has a fixture-verified `0.15-iupac` no-opt mapping branch and uses it to build the `0.15-iupac` CRAM/report track after consensus.
- The workflow now also wires the `${sample}-${subtype}` best-reference branch from the early `MAP_READS` outputs: it uses `MAP_READS.out.bams_with_ref` joined to `ch_best_ref_with_name`, runs `BAM2FASTA_BESTREF` to recover `${sample}-${subtype}.fasta` and stats, builds `${sample}-${subtype}.cram` against the original chosen reference FASTA, and reports with `report_id = ${sample}-${subtype}`.

### Clean Run

```bash
# Fresh run with cache clear
nextflow run . -profile hpc,apptainer -resume -clear

# Clean only work directory
rm -rf work/ .nextflow/
```

---

## Code Style Guidelines

### Process Definitions

```nextflow
process PROCESS_NAME {
    tag { "${sample_id}" }           // Required for traceability
    label 'process_medium'           // Resource label from base.config
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    // Publish only specific outputs, not all files
    publishDir "${params.outdir}/${run_name}/${sample_id}/subdir", 
               mode: 'copy', 
               pattern: '*.fastq.gz'
    
    input:
        // Use tuple for multi-inputs; name each element
        tuple val(run_name), val(sample_id), path(read1), path(read2)
        val(optional_config)
    
    output:
        // Always emit named channels for clarity
        tuple val(run_name), val(sample_id), path("*.bam"), emit: bams
        path "*.stats", emit: stats
        path "*.json", emit: json
    
    script:
    // Use Groovy string interpolation for params
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    """
    # Always use bash with strict mode
    set -euo pipefail
    
    # Container execution
    apptainer exec -B ${bind_paths} ${container_dir}/tool.sif command --args
    """
}
```

### Module Organization

```
modules/local/<module_name>/
└── main.nf          // Single process or sub-workflow
```

**Include modules using DSL2:**
```nextflow
include { PROCESS_NAME } from './modules/local/module_name/main'
```

### Naming Conventions

| Element | Convention | Example |
|---------|------------|---------|
| Processes | UpperCamelCase | `REMOVE_HOSTILE`, `BAM2FASTA` |
| Workflows | UpperCamelCase | `HCVPIPE`, `SUBSAMPLE` |
| Channels | LowerCamelCase | `ch_samples`, `ch_reads` |
| Parameters | snake_case | `container_dir`, `bind_paths` |
| Files | kebab-case | `sample-id.fastq.gz` |

### Input/Output Channels

**Always use named emits:**
```nextflow
output:
    tuple val(run_name), val(sample_id), path("*.bam"), emit: bams
    path "*.stats", emit: stats
```

**Prefer tuples over bare paths:**
```nextflow
// Good: carries metadata with data
tuple val(run_name), val(sample_id), path(read1), path(read2)

// Avoid: loses context
path read1
```

### Container Directives

```nextflow
// Prefer container directive over inline apptainer
process Foo {
    container 'image:tag'  // If available in registry
    
    // Or use apptainer exec in script for custom containers
    script:
    def cmd = "apptainer exec -B ${bind_paths} ${container_dir}/tool.sif tool"
}
```

### Shell Script Handling

**Always use strict bash mode:**
```bash
set -euo pipefail

# Don't use 'sh -c "..."' wrapper
# Use direct bash in process script block
```

**Avoid command substitution in process definition:**
```nextflow
// Bad - evaluated at parse time
def ref = $(echo ${params.ref_dir}/genome.fa)

// Good - evaluated at execution time
script:
def ref = "${params.ref_dir}/genome.fa"
```

### Error Handling

**Use Nextflow's built-in error strategy:**
```nextflow
process {
    errorStrategy = 'terminate'  // Default - fail fast
    maxRetries = 1
    
    // For robust processes
    withLabel: robust {
        errorStrategy = 'retry'
        maxRetries = 3
    }
}
```

**Validate inputs at workflow level:**
```nextflow
workflow {
    if (!params.input) {
        error "Missing required parameter: --input <samplesheet.csv>"
    }
    
    Channel.fromPath(params.input, checkIfExists: true)
        .ifEmpty { error "Input file not found: ${params.input}" }
}
```

---

## Project-Specific Patterns

### HPC Environment (Hopper)

**Paths:**
- Test data: `/fs1/jonas/hcv/test_data/`
- Reference genomes: `/fs1/jonas/hcv/refgenomes/`
- Containers: `/fs1/resources/containers/`
- Hostile cache: `/fs1/resources/ref/micro/hostile`

**Modules to load:**
```bash
source /etc/profile
module load Java/23.0.2 nextflow/25.10.0 apptainer
```

### Samplesheet Format

```csv
clarity_sample_id,read1,read2,sample_name,run_name
SAMPLE001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,LID001,test
```

**Expected columns:**
- `clarity_sample_id` - Sample ID (required)
- `read1` / `fastq_1` - Read 1 path (required)
- `read2` / `fastq_2` - Read 2 path (optional)
- `sample_name` / `lid` - Sample name/LID (optional)
- `run_name` - Run name for output directory (optional)

### Configuration

**Global config:** `nextflow.config`  
**Process config:** `conf/base.config`  
**Container dir:** Set via `--container_dir` param

**Key params:**
```nextflow
params {
    input = null              // Samplesheet
    outdir = "results"        // Output directory
    container_dir = '/fs1/resources/containers'
    bind_paths = '/fs1,/fs2,/local'
    cpus = 16
    memory = '32 GB'
    time = '24h'
}
```

---

## Testing Approach

### Module Testing

1. Add module to `test_module.nf`:
```nextflow
include { MODULE_NAME } from './modules/local/module_name/main'

workflow {
    Channel.fromList([...]).set { ch_input }
    MODULE_NAME(ch_input, ...)
}
```

2. Run on HPC:
```bash
bash runme_test.sh <module_name>
```

3. Compare outputs:
```bash
# Compare file checksums
md5sum bash_output/* nextflow_output/*

# Compare statistics
diff <(jq '.reads_in' bash.json) <(jq '.reads_in' nextflow.json)
```

### Output Parity Verification

For each module:
1. Run bash version and capture output
2. Run Nextflow version
3. Compare file contents (not just checksums)
4. Verify statistics in JSON/metadata files

**Known acceptable differences:**
- Read ordering in parallel-processed FASTQ (data identical, order non-deterministic)
- Embedded paths in VCF headers (different container paths)

---

## Branching Strategy

- **Main development:** `refactor/modular-take2`
- **Push to:** `hopper` bare repo on HPC
- **Command:** `git push --set-upstream hopper refactor/modular-take2`

---

## File Locations

| Purpose | Path |
|---------|------|
| Main entry | `main.nf` |
| Workflow | `workflows/hcvpipe.nf` |
| Modules | `modules/local/*/main.nf` |
| Config | `nextflow.config`, `conf/base.config` |
| Test workflow | `test_module.nf` |
| Test runner | `runme_test.sh` |
| Pipeline docs | `docs/pipeline_workflow.md` |
