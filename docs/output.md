# Output

## Analytical IUPAC IGV review bundle

Finalization adds `<sample>-0.15-iupac.review.*` display copies of resistance,
VADR and coverage annotations and every filtered VCF, with fresh CSI indexes.
Original files remain unchanged. The genome/read pair is the existing
`<sample>-0.15-iupac.fasta` / `.cram`, with their indexes. Only contig identifiers
and applicable header declarations change; positions, alleles and evidence do not.

`<sample>-0.15-iupac.review.json` records schema version 1.0, source/display
reference identities and sequence checksums, shared coordinate provenance,
source/display track checksums and analytical roles. QC JSON embeds this as
`igv_review`, with additive flat `igv_iupac_*` output paths and
`igv_review_manifest`. Resistance remains calculated from IUPAC; VADR/coverage
and VCF tracks are scaffold-derived. VCF REF can differ from the displayed
IUPAC base: these are display copies, not newly called variants.

Missing or incompatible provenance, reference/read dictionaries, indexes or
source-track contigs fail preparation. Historical archives are not automatically
converted. Resistance re-calling regenerates existing review bundles within the
same backup/restore transaction.

For each sample, the pipeline publishes a flat Virtitta-compatible final archive under:

- `results/<run_name>/<sample_id>/`

Example:

- `results/240101_A00000_0001_XXXXXX/SAMPLE001/`

`--publish_mode routine` is the default. `--publish_mode debug` keeps the same final files at the
sample root and also publishes intermediate subdirectories such as `bam/`, `vcf/`, `fastq/`,
`fasta/`, `spades/`, `mummer/`, and `pilon/`.

Subsampled FASTQs are published under `fastq/subsample_<requested_pairs>/`, keeping
reference-selection and polishing subsets separate without changing their filenames.
The directory names record the requested limit, not the actual retained pair count
when fewer pairs are available. Hostile FASTQs remain directly under `fastq/`.
Routine mode publishes neither subsampling directory.

Current final sample outputs include:

- `SAMPLE001.fasta` and `SAMPLE001.fasta.fai`
- `SAMPLE001.fasta.consensus.json` (reference-filled mapping scaffold provenance)
- `SAMPLE001-majority.fasta`, `.fasta.fai` and `.fasta.consensus.json` (masked analytical majority)
- `SAMPLE001.cram` and `SAMPLE001.cram.crai`
- `SAMPLE001.fasta.blast`
- `SAMPLE001-0.15-iupac.fasta`
- `SAMPLE001-0.15-iupac.fasta.fai`
- `SAMPLE001-0.15-iupac.fasta.consensus.json`
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
- `SAMPLE001-coverage-1x.bed`
- `SAMPLE001.vadr.pass_mod.gff`
- `SAMPLE001.vadr.bed`
- `SAMPLE001_resistance.tsv`
- `SAMPLE001_resistance.bed`
- `SAMPLE001_resistance.gff`
- `SAMPLE001_resistance_by_drug.tsv`
- `SAMPLE001_resistance.json`
- `SAMPLE001_resistance_sites.gff3`
- `SAMPLE001_display_rug_kde_plot.png`
- `SAMPLE001_qc_summary.json`
- `hostile.json` when host filtering is enabled

Resistance JSON uses schema `3.1`; its outcome, low-depth qualification, site assessment, provenance, integrity and
gene-alignment data are carried through the per-sample QC JSON for Virtitta. The workflow's
resistance tracks target the main FASTA contig for a combined IGV view. Analytical IUPAC
FASTA headers contain identifiers only. Converter-4 provenance is in accompanying consensus JSON, bound to uppercase sequence checksums and shared coordinate origins. Analytical IUPAC and majority outputs mask zero coverage, retain observations at depth 1–6, and preserve coordinate length; `main_fasta` remains the reference-filled CRAM/IGV scaffold. The new `analytical_majority_fasta` output is preferred for exports and majority distance calculations. QC JSON contains compact `consensus` summaries and `main_consensus_metadata`, `iupac_consensus_metadata` and `analytical_majority_consensus_metadata` output paths. See `docs/resistance_pipeline.md` for assessment and recall safeguards.

The per-sample `*_qc_summary.json` is the machine-readable downstream contract for analysis tools.
It includes stable identifiers such as `sample_run_id`, pipeline metadata, extracted QC metrics, and
relative paths to key result files needed for tables, detail views, and IGV launchers.
`pipeline_version` records the producing Nextflow manifest version (`0.3.0` for the
three-FASTA analytical-consensus contract). It is additive to QC schema 4. Standalone
QC generation without `--pipeline-version` records null rather than inferring a version.
Resistance-only recall preserves the original producing pipeline version; caller,
consensus, rule and reference provenance describe the resistance evidence separately.
In `typing`, `main_blast_genotype` remains the historical scaffold top-hit field.
`best_mapping_subtype` is the subtype prefix of `report_reference`, the mapping-selected reference;
these can differ from the BLAST call.

New runs additionally publish `<sample>_typing.json` (`outputs.typing_json`) and embed
it as `typing.evidence` in QC schema 4. This is the authoritative analytical 15% IUPAC
typing result used for resistance. Consumers must prefer it when present, even when
its subtype is null; legacy fields must not fill in an uncertain new assignment.
Evidence schema `1.0` records `method`, `status` (`assigned`, `genotype_only`,
`ambiguous`, `unassigned`), nullable `subtype` and `selected_reference`, ranked `groups`,
`unrecognized_references`, descriptive `score_margin`, query lengths, sequence checksum
and role in `input`, and database filenames/checksums in `database`. Each group retains
its references and selected segments, with summed bit score, weighted identity, aligned
bases and whole-query/observed-query coverage. See the BLAST plan for metric definitions.
No-hit and equal-score cross-subtype ties still publish a complete sample; ambiguous or
unassigned typing produces 16 not-assessed drug outcomes with unknown (`null`) licensing.
Raw BLAST searches append exact `nident` as column 13; their first 12 columns retain
standard BLAST meanings. The target limit equals the database entry count, not a top-N
display limit. Mapping-reference report fields and filenames retain their existing role.
All QC JSON `outputs` paths are relative to the directory containing the QC JSON. Standard sidecars
such as `<main_fasta>.fai`, `<main_cram>.crai`, and `<filtered_vcf>.csi` are published on disk but
not listed as separate JSON keys when their names are inferable. LID-specific duplicate FASTA
exports and `lid_2limsrs` are no longer published; Virtitta rewrites FASTA headers at export time.
When `sample_name` / `lid` is present, `display_rug_kde_plot` points to
`<sample_id>_display_rug_kde_plot.png`, titled `<LID> (<sample_id>)`.

Consensus FASTA sequence lines are uniformly uppercase. This is an intentional case-only departure
from legacy output; hybrid-assembly and Pilon debug/intermediate FASTAs retain their provenance case.
The `<sample_id>-coverage-1x.bed` file contains merged, 0-based half-open intervals where the final
main CRAM has at least one aligned read base against the main FASTA. The existing coverage TSV remains
the aggregate percentage summary and is unchanged.

The per-run QC summary is written to:

- `results/<run_name>/pipeline_info/qc_summary.json`
- `results/<run_name>/pipeline_info/qc_summary.jsonl`

If enabled with `--pipeline_info`, Nextflow run metadata is also written to:

- `results/<run_name>/pipeline_info/timeline.html`
- `results/<run_name>/pipeline_info/report.html`
- `results/<run_name>/pipeline_info/trace.txt`

If `--run_name` is not provided, the built-in metadata path falls back to `results/pipeline_info_<timestamp>/`. DAG rendering is disabled by default; add `--pipeline_info_dag` to request `dag.svg` when Graphviz is available.

If `--completion_log_dir <dir>` is provided, the pipeline also writes a legacy control-system
completion summary to `<dir>/<run_name>-HCV.complete`. This file is intended for external
launchers and is not written by default.
