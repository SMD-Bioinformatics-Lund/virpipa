# Usage

## Input

Provide a CSV samplesheet with at least these columns:

- `clarity_sample_id`
- `read1`
- `read2` (optional; if omitted, `hcvpipe.sh` will infer R2 from R1 naming)
- `sample_name` (optional; passed to `hcvpipe.sh` as LID)
- `run_name` (optional; folder level under `outdir`)
- `sequencing_run` (optional; alias for `run_name`)

Run naming precedence is:

- `run_name` column in the CSV
- `sequencing_run` column in the CSV
- `--run_name`
- fallback to `test`

## Run

```bash
nextflow run . \
  --input samplesheet.csv \
  # or: --csv samplesheet.csv \
  --outdir results \
  --run_name 260224_A00681_1214_BHHG2YDRX7 \
  --partition grace-normal \
  --container_dir /fs1/resources/containers \
  --bind_paths '/fs1,/fs2,/local' \
  -profile slurm,apptainer
```

## Notes

- `--csv` is accepted as an alias for `--input`.
- `--queue` is accepted as an alias for `--partition` for site launchers that use queue terminology for SLURM partitions.
- `--container_runtime` can be set to `apptainer` or `singularity`; if unset, module scripts prefer `apptainer` when available and otherwise fall back to `singularity`.
- `--sample_info_json` overrides the default nearby-CSV lookup for `clarity_sample_info.json` and is useful when the samplesheet has been restored or relocated.
- If `--sample_info_json` is not provided, the pipeline looks for `clarity_sample_info.json` in the samplesheet directory, then its parent, then its grandparent.
- This scaffold intentionally runs the current `scripts/hcvpipe.sh` as a single process.
- Next migration step is to split that process into modular DSL2 components.
