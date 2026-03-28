# virpipa

VirPipa is a modular Nextflow DSL2 pipeline for HCV probe-capture assembly and reporting. It reproduces the bash-original `results/` contract while adding a maintainable workflow graph, fixture-backed module tests, VADR annotation, and geno2pheno resistance annotation.

## Overview

The pipeline takes paired-end FASTQ input, removes host contamination, selects the best HCV reference, builds a hybrid assembly, polishes it with a pilon loop, creates final consensus and CRAM outputs, annotates with VADR, and reports subtype and resistance calls.

Main outputs per sample include:

- polished consensus FASTA and FAI
- `0.15-iupac` consensus FASTA
- best-reference FASTA, CRAM, and report branch
- filtered `pilon-m*` VCF set
- CRAM / CRAI files
- subtype BLAST outputs
- VADR GFF and BED
- coverage, report, and nucleotide-frequency TSVs
- geno2pheno resistance TSV / BED / by-drug TSV

## Workflow

```text
FASTQ
  -> hostile host removal
  -> subsample for best-reference selection
  -> map to all HCV references
  -> select best reference
  -> SPAdes assembly
  -> hybrid reference build
  -> pilon polishing loop
  -> pilon BAM -> majority-call replacement FASTA
  -> pilon VCF -> filtered pilon-m* VCFs
  -> 0.15-iupac consensus
  -> no-opt remap to 0.15-iupac
  -> final CRAM / report / BLAST branches
  -> VADR annotation
  -> geno2pheno resistance annotation
  -> finalize bash-style results/
```

## Quick Start

### Hopper

Push your current branch to the Hopper bare repo before starting a workflow there:

```bash
git push --set-upstream hopper refactor/modular-take3
```

Run the standard Hopper wrapper from the cluster login node:

```bash
ssh rs-fe1 'cd /fs1/jonas/src/virpipa && bash -l runme.sh resume'
```

For a test-module run on Hopper:

```bash
bash runme_test.sh resistance
```

### Local Laptop

Use the `skrotis` environment for local development:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate skrotis
```

Run a local module test:

```bash
nextflow run test_module.nf -profile local --module resistance
```

Run the full workflow locally with containers:

```bash
export SENTIEON_LICENSE=localhost:8990
NXF_OFFLINE=true nextflow run . \
  -profile local_containers \
  --input /tmp/virpipa_contract_local_samplesheet.csv \
  --outdir /tmp/virpipa_contract_results_local
```

## Inputs

Samplesheet columns:

- `clarity_sample_id` or equivalent sample id column
- `read1`
- `read2`
- optional `sample_name` or `lid`
- optional `run_name` or `sequencing_run`

Example:

```csv
clarity_sample_id,read1,read2,sample_name,run_name
SAMPLE001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,TEST001,test_run
```

## Outputs

Sample outputs are published under:

- `results/<run_name>/<sample_id>/results/`

Run metadata is written under:

- `results/pipeline_info/timeline.html`
- `results/pipeline_info/report.html`
- `results/pipeline_info/trace.txt`
- `results/pipeline_info/dag.svg`

See [`docs/output.md`](/home/jonas/git/virpipa/docs/output.md) for the current output contract and [`docs/pipeline_workflow.md`](/home/jonas/git/virpipa/docs/pipeline_workflow.md) for workflow details.

## Resistance Rules

geno2pheno HCV rules are refreshed manually outside Hopper. The offline refresh entrypoint is:

```bash
python scripts/update_geno2pheno_rules.py --output-csv assets/hcv_geno2pheno_rules.csv
```

The pipeline currently consumes the committed CSV rules table by default via `params.resistance_rules` and the repo default is [`assets/hcv_geno2pheno_rules.csv`](/home/jonas/git/virpipa/assets/hcv_geno2pheno_rules.csv).
