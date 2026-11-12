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
- `--completion_log_dir` writes a legacy `<run_name>-HCV.complete` workflow summary for external control systems.
- `--virtitta_import_dir` writes a successful-run `<run_name>.sqlimport` marker for the Virtitta importer. It defaults to `/fs1/results/cron/virtitta`; set it to an empty string to disable this marker. The marker points to `<virtitta_run_dir>/<run_name>/`, where `--virtitta_run_dir` defaults to `/access/virpipa/hcv`.

## Output routing

Without `--dev` (or with `--dev false`), shared output destinations use the
production results root `/fs1/results`. `--dev` selects `/fs1/results_dev`:

| Parameter | Production default | Development default |
| --- | --- | --- |
| `outdir` | `/fs1/results/virpipa/hcv` | `/fs1/results_dev/virpipa/hcv` |
| `completion_log_dir` | `/fs1/results/cron/logs` | `/fs1/results_dev/cron/logs` |
| `virtitta_import_dir` | `/fs1/results/cron/virtitta` | `/fs1/results_dev/cron/virtitta` |
| `virtitta_run_dir` | `/access/virpipa/hcv` | The resolved `outdir` |

Production markers use `/access/virpipa/hcv`, the importer-visible results path.
Production launchers no longer need directory overrides after deploying this
version. Completion logging is enabled by default for shared runs; an explicitly
empty `completion_log_dir` disables it.

The laptop profiles `local` and `local_containers` use relative `results` and
empty completion-log and marker-publication directories, with or without `--dev`.
Their marker target defaults to the resolved `outdir` if publication is explicitly
enabled. Profiles `tiny` and `test` retain their result-directory overrides;
combine them with a local profile for laptop tests. Profile order follows
Nextflow precedence, so use `local,tiny` or `local_containers,test`.

Explicit directory overrides take precedence. In particular, `--outdir` also
sets the default dev marker target; `--virtitta_run_dir` overrides that target.
An explicitly empty `virtitta_import_dir` in a config or params file still
disables markers. For command-line disabling, use `--virtitta_import_dir=` so
the empty value cannot be interpreted as a boolean flag.

Selecting an assay named `HCV-dev` in the launcher does not itself activate
pipeline dev mode; the pipeline command must include `--dev`.

See [launcher integration](deployment.md#launcher-integration) for the replacement
launcher entry and [completed-run recovery](dev_output_recovery.md) for existing results.
