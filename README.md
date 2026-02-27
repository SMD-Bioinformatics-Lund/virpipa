# virpipa
![Virpipare](https://cdn.pixabay.com/photo/2012/05/02/22/25/bird-46531_960_720.png)

Virus assembly pipeline. First for HCV.

## Current layout

- `scripts/`: existing Bash pipeline (`hcvpipe.sh`) and helper scripts.
- `main.nf`, `workflows/`, `modules/local/`, `conf/`: nf-core style Nextflow scaffold.

## Run Nextflow scaffold

```bash
nextflow run . \
  --input samplesheet.csv \
  --outdir results \
  --run_name 260224_A00681_1214_BHHG2YDRX7 \
  --container_dir /fs1/resources/containers \
  --bind_paths '/fs1,/fs2,/local' \
  -profile slurm,apptainer
```

The scaffold currently wraps `scripts/hcvpipe.sh` as one process to preserve behavior while preparing modular migration.
