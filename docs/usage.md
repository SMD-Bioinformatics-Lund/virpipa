# Usage

## Input

Provide a CSV samplesheet with at least these columns:

- `clarity_sample_id`
- `read1`
- `read2` (optional; if omitted, `hcvpipe.sh` will infer R2 from R1 naming)
- `sample_name` (optional; passed to `hcvpipe.sh` as LID)

## Run

```bash
nextflow run . \
  --input samplesheet.csv \
  --outdir results \
  --container_dir /fs1/resources/containers \
  --bind_paths '/fs1,/fs2,/local' \
  -profile slurm,apptainer
```

## Notes

- This scaffold intentionally runs the current `scripts/hcvpipe.sh` as a single process.
- Next migration step is to split that process into modular DSL2 components.
