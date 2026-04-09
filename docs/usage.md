# Usage

## Input

Provide a CSV samplesheet with at least these columns:

- `clarity_sample_id`
- `read1`
- `read2` (optional; if omitted, `hcvpipe.sh` will infer R2 from R1 naming)
- `sample_name` (optional; passed to `hcvpipe.sh` as LID)
- `run_name` (optional; folder level under `outdir`)

If `run_name` is omitted, the pipeline infers it from `read1` using:

- `.../<run_name>/Data/Intensities/BaseCalls/...`

## Run

```bash
nextflow run . \
  --input samplesheet.csv \
  # or: --csv samplesheet.csv \
  --outdir results \
  --run_name 260224_A00681_1214_BHHG2YDRX7 \
  --container_dir /fs1/resources/containers \
  --bind_paths '/fs1,/fs2,/local' \
  -profile slurm,apptainer
```

## Notes

- `--csv` is accepted as an alias for `--input`.
- This scaffold intentionally runs the current `scripts/hcvpipe.sh` as a single process.
- Next migration step is to split that process into modular DSL2 components.
