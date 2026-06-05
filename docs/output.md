# Output

For each sample, the pipeline publishes a flat Virtitta-compatible final archive under:

- `results/<run_name>/<sample_id>/`

Example:

- `results/240101_A00000_0001_XXXXXX/SAMPLE001/`

`--publish_mode routine` is the default. `--publish_mode debug` keeps the same final files at the
sample root and also publishes intermediate subdirectories such as `bam/`, `vcf/`, `fastq/`,
`fasta/`, `spades/`, `mummer/`, and `pilon/`.

Current final sample outputs include:

- `SAMPLE001.fasta` and `SAMPLE001.fasta.fai`
- `SAMPLE001.cram` and `SAMPLE001.cram.crai`
- `SAMPLE001.fasta.blast`
- `SAMPLE001-0.15-iupac.fasta`
- `SAMPLE001-0.15-iupac.fasta.fai`
- `SAMPLE001-0.15-iupac.cram` and `SAMPLE001-0.15-iupac.cram.crai`
- `SAMPLE001-0.15-iupac.report.tsv`
- `SAMPLE001-0.15-iupac.fastanucfreq.tsv`
- `SAMPLE001-0.15-iupac.fasta.blast`
- `SAMPLE001-<subtype>.fasta`
- `SAMPLE001-<subtype>.cram` and `SAMPLE001-<subtype>.cram.crai`
- `SAMPLE001-<subtype>.report.tsv`
- `SAMPLE001-<subtype>.fastanucfreq.tsv`
- `SAMPLE001-<subtype>.vcf.gz`, `.csi`, `.stats`
- full `SAMPLE001-pilon-m*.vcf.gz`, `.csi`, `.stats` set
- `SAMPLE001-pilon-iupac.fasta.blast`
- `SAMPLE001-coverage.tsv`
- `SAMPLE001.vadr.pass_mod.gff`
- `SAMPLE001.vadr.bed`
- `SAMPLE001_resistance.tsv`
- `SAMPLE001_resistance.bed`
- `SAMPLE001_resistance.gff`
- `SAMPLE001_resistance_by_drug.tsv`
- `SAMPLE001_display_rug_kde_plot.png`
- `SAMPLE001_qc_summary.json`
- `hostile.json` when host filtering is enabled

The per-sample `*_qc_summary.json` is the machine-readable downstream contract for analysis tools.
It includes stable identifiers such as `sample_run_id`, pipeline metadata, extracted QC metrics, and
relative paths to key result files needed for tables, detail views, and IGV launchers.
All QC JSON `outputs` paths are relative to the directory containing the QC JSON. Standard sidecars
such as `<main_fasta>.fai`, `<main_cram>.crai`, and `<filtered_vcf>.csi` are published on disk but
not listed as separate JSON keys when their names are inferable. LID-specific duplicate FASTA
exports and `lid_2limsrs` are no longer published; Virtitta rewrites FASTA headers at export time.
When `sample_name` / `lid` is present, `display_rug_kde_plot` points to
`<sample_id>_display_rug_kde_plot.png`, titled `<LID> (<sample_id>)`.

The per-run QC summary is written to:

- `results/<run_name>/pipeline_info/qc_summary.json`
- `results/<run_name>/pipeline_info/qc_summary.jsonl`

If enabled with `--pipeline_info`, Nextflow run metadata is also written to:

- `results/<run_name>/pipeline_info/timeline.html`
- `results/<run_name>/pipeline_info/report.html`
- `results/<run_name>/pipeline_info/trace.txt`

If `--run_name` is not provided, the built-in metadata path falls back to `results/pipeline_info_<timestamp>/`. DAG rendering is disabled by default; add `--pipeline_info_dag` to request `dag.svg` when Graphviz is available.

If `--completion_log_dir <dir>` is provided, the pipeline also writes a legacy control-system
completion summary to `<dir>/<samplesheet-base>.complete`. This file is intended for external
launchers and is not written by default.
