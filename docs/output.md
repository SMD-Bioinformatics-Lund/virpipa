# Output

For each sample, the pipeline publishes a bash-compatible final results tree under:

- `results/<run_name>/<sample_id>/results/`

Example:

- `results/240101_A00000_0001_XXXXXX/SAMPLE001/results/`

Current final sample outputs include:

- `SAMPLE001.fasta` and `SAMPLE001.fasta.fai`
- `SAMPLE001.cram` and `SAMPLE001.cram.crai`
- `SAMPLE001.fasta.blast`
- `SAMPLE001-0.15-iupac.fasta`
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
- `SAMPLE001_qc_summary.json`
- `hostile.json` when host filtering is enabled

The per-sample `*_qc_summary.json` is the machine-readable downstream contract for analysis tools.
It includes stable identifiers such as `sample_run_id`, pipeline metadata, extracted QC metrics, and
relative paths to key result files needed for tables, detail views, and IGV launchers.
When LID-specific outputs exist under `results/lid/`, the QC JSON also exposes those paths together with
preferred export/display keys such as `export_fasta`, `export_iupac_fasta`, and `display_rug_kde_plot`.

Optional bash-style LID outputs are published when `sample_name` / `lid` is present:

- `results/<lid>.lid`
- `results/lid/<lid>.fasta`
- `results/lid/<lid>-0.15-iupac.fasta`
- `results/lid/<lid>_rug_kde_plot.png`
- `results/lid/<lid>-2limsrs.txt`

Pipeline run metadata is written to:

- `results/<run_name>/pipeline_info/timeline.html`
- `results/<run_name>/pipeline_info/report.html`
- `results/<run_name>/pipeline_info/trace.txt`
- `results/<run_name>/pipeline_info/dag.svg`
- `results/<run_name>/pipeline_info/qc_summary.json`
- `results/<run_name>/pipeline_info/qc_summary.jsonl`
