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
  --publish_mode routine \
  --partition grace-normal \
  --container_dir /fs1/resources/containers \
  --bind_paths '/fs1,/fs2,/local' \
  -profile slurm,apptainer
```

## Notes

- `--csv` is accepted as an alias for `--input`.
- `--queue` is accepted as an alias for `--partition` for site launchers that use queue terminology for SLURM partitions.
- `--container_runtime` can be set to `apptainer` or `singularity`; if unset, module scripts prefer `apptainer` when available and otherwise fall back to `singularity`.
- `--publish_mode routine` is the default and publishes flat final sample archives under `<outdir>/<run_name>/<sample_id>/`.
- `--publish_mode debug` keeps the same flat final sample archive and additionally publishes intermediate subdirectories for inspection.
- `--sample_info_json` overrides the default nearby-CSV lookup for `clarity_sample_info.json` and is useful when the samplesheet has been restored or relocated.
- If `--sample_info_json` is not provided, the pipeline looks for `clarity_sample_info.json` in the samplesheet directory, then its parent, then its grandparent.
- `--completion_log_dir` optionally writes a legacy `<run_name>-HCV.complete` workflow summary for external control systems.
- `--virtitta_import_dir` writes a successful-run `<run_name>.sqlimport` marker for the Virtitta importer. It defaults to `/fs1/results/cron/virtitta`; set it to an empty string to disable this marker.
