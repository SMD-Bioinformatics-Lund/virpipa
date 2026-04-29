# Agent Environment Reference

This file holds detailed environment and operations notes that are useful when working in `virpipa`, but are too bulky for always-on `AGENTS.md` instructions.

## Execution environments

### Laptop

- The agent runs locally on a laptop.
- Use the `skrotis` mamba environment for local Nextflow and tool runs.
- The laptop does not have enough CPU/RAM for routine full-sample runs.
- Prefer tiny data or fixture-backed module tests locally.
- Local container images are stored in `assets/containers/`.
- Hopper `/fs1` is mounted locally at `/mnt/fs1`.

Helpful local paths:

- Faster hostile cache: `~/resources/hostile`
- Faster VADR models: `~/resources/vadr/vadr-models-flavi`
- Faster bash-original comparison copy: `~/hcv/SAMPLE001-nextflow-nfcore-scaffold/`
- Additional local test data: `~/hcv/test_data/`

### Hopper

- Hopper has no internet.
- Shared storage root is `/fs1`.
- Production containers are under `/fs1/resources/containers/`.
- Hostile cache is under `/fs1/resources/ref/micro/hostile`.
- VADR models are under `/fs1/resources/ref/micro/vadr/vadr-models-flavi`.
- The site wrappers configure the environment; do not assume generic `nextflow run` on Hopper is equivalent.

## Running on Hopper

Preferred wrapper invocations:

```bash
ssh rs-fe1 'cd /fs1/jonas/src/virpipa && bash -l runme.sh resume'
ssh rs-fe1 'cd /fs1/jonas/src/virpipa && bash -l runme.sh'
```

Module test wrapper:

```bash
bash runme_test.sh <module_name>
```

Important workflow for code updates:

- Commit locally.
- Push to the Hopper bare repo.
- Let the Hopper wrapper pull/update as part of its normal flow.
- Do not manually pull in the Hopper working clone unless there is a specific recovery reason.

## Current branch / repo expectations

- `master` is the active base branch after the first release.
- For isolated work, create a new branch from `master` rather than continuing old refactor branches.

## Launcher and CLI notes

Current useful aliases/overrides:

- `--csv` is an alias for `--input`.
- `--queue` is an alias for `--partition`.
- `--sample_info_json` overrides the nearby-samplesheet lookup for `clarity_sample_info.json`.
- `--container_runtime` can be set to `apptainer` or `singularity`; if unset, module scripts prefer `apptainer` and otherwise fall back to `singularity`.

## Built-in pipeline metadata

VirPipa's built-in `timeline/report/trace/dag` output is opt-in via:

- `--pipeline_info`
- `--pipeline_info_dag`

Many site launchers already manage metadata output separately. In that case, leave these disabled.

## Relevant files

- Main entry: `main.nf`
- Workflow: `workflows/hcvpipe.nf`
- Module tests: `test_module.nf`
- Wrapper scripts: `runme.sh`, `runme_test.sh`
- Legacy bash pipeline: `scripts/hcvpipe.sh`
