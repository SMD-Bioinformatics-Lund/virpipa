# Consensus and distance benchmark

This is an evaluation on branch `consensus-benchmark`, based on `78cd1e0`.
The production converter, resistance caller and published archives are unchanged.
Do not deploy a replacement based solely on tool agreement.

## Acceptance target

The primary resistance endpoint is observed amino-acid prevalence >=15%.
Sum synonymous codons before thresholding and retain their separate nucleotide
and codon counts. IUPAC translation can include unobserved amino acids; equality
to the current IUPAC output is therefore not an independent truth standard.

Compare the current converter with BCFtools 1.21 and iVar 1.4.4/SAMtools 1.21.
The native haploid BCFtools genotype and an ALT-application diagnostic are
separate arms, neither advertised as an allele-depth-based minority caller.
iVar thresholds 0, 0.85 and 0.90 are native cumulative consensus settings;
0.15 is a diagnostic, not a claim of equivalence to our per-allele 15% policy.

Sources: [HCV guidance](https://www.hcvguidelines.org/guidance/hcv-resistance-primer/),
[geno2pheno codon-frequency recommendation](https://academic.oup.com/nar/article/46/W1/W271/4990638),
[MiCall counts and consensus outputs](https://cfe-lab.github.io/MiCall/steps.html),
[BCFtools consensus](https://www.htslib.org/doc/bcftools.html#consensus),
[iVar manual](https://github.com/andersen-lab/ivar/blob/master/docs/MANUAL.md).

## Runnable probes

Activate `skrotis`. Tool commands are argument arrays, without shell evaluation.
For example, write a private JSON file containing:

```json
{
  "bcftools": ["apptainer", "exec", "assets/containers/bcftools_1.21.sif", "bcftools"],
  "samtools": ["apptainer", "exec", "assets/containers/samtools_1.21.sif", "samtools"],
  "ivar": ["apptainer", "exec", "/path/to/ivar-1.4.4.sif", "ivar"]
}
```

Use absolute image paths when running outside the repository. Container images
must never be committed. Omitting `ivar` permits a clearly labelled partial run.

```bash
python scripts/benchmark_consensus.py --tools /path/to/tools.json --output /tmp/consensus-probes-NEW --prepare-distances
python -m unittest discover -s tests -p test_consensus_benchmark.py -v
```

The output directory must not exist. The runner creates known molecules, indexed
BAMs, truth counts, native tool outputs, commands, versions, hashes and a JSON
scorecard. The 36 current controls cover 14/15/16% boundaries, synonymous
aggregation, an unobserved-codon mixture, multiallelic calls, missing/low depth,
quality exclusions, agreeing/conflicting templates, partial codons and 1/2/3-base
indels at 8/15/50/90% support. These BAMs have explicit CIGARs: they test calling
and conversion, not the accuracy of a read mapper. BAQ is explicitly disabled
for both tool paths in this controlled experiment; production retains its policy.
The raw-VCF baseline and normalized-VCF diagnostic are separate outputs: otherwise
normalization would silently change the baseline's masked reference spans.
`--prepare-distances` also writes a manifest and masks from the simulated BAMs.
Length-changing outputs are explicitly excluded from this diagnostic manifest,
not padded or assigned invented coverage.

## Archived cohort

Supply an explicit private text file of the 40 run directories, one per line.
Laptop paths use `/mnt/fs1`; Hopper uses `/fs1`. The inventory rejects missing
inputs, duplicate sample/run identities and a count other than 240 by default.

```bash
python scripts/benchmark_cohort.py --run-dirs /path/to/runs.txt --output /tmp/cohort-NEW --evidence
```

This hashes source FASTAs, alignments, indexes, coverage, resistance and QC files.
It compares complete-codon evidence from scaffold and IUPAC alignments using the
existing caller's extraction rules, and verifies that source hashes remain
unchanged. Different evidence does not by itself establish a changed drug call.
Run-level result directories and per-sample names remain in private output only.

## Distance comparison

Reuse Virtitta's actual distance functions, pinned by source hash and Git commit.
Supply a JSON manifest with one row per sample/method/mode:

```json
[
  {"sample_run_id": "ANON_RUN", "method": "baseline", "mode": "iupac", "subtype": "1a", "fasta": "/path/to/prepared.fa", "mask": "/path/to/prepared-mask.json"}
]
```

Include both `iupac` and `majority` rows for operational comparisons. Each FASTA
contains one unaligned sequence. A mask is a JSON array of booleans in that
sequence's coordinates. Apply the configured Virtitta trimming to both sequence
and evidence mask before constructing the manifest; record that preparation
command. Derive masks from read evidence, not FASTA letter case. For edited
indels, provide explicit projected coverage evidence; do not reuse scaffold
positions by index or fabricate coverage for inserted bases.

```bash
python scripts/benchmark_distances.py --manifest /path/to/prepared.json --output /tmp/distances-NEW --virtitta-source /path/to/virtitta --mafft /path/to/mafft
```

Operational majority comparisons reuse the IUPAC gap pattern, as Virtitta does.
Missing or length-incompatible pairs are reported in `incompatible.json`.
Pooled alignments include all versions within each subtype/mode. A third
comparison holds covered, valid nucleotide positions constant across both
methods for each pair; it deliberately excludes indel columns, which remain
represented in the native/pooled event counts. Uncertain subtypes form a separate
diagnostic group. Missing-overlap pairs never count as nearest identical samples.

Outputs include per-pair counts, bases compared, deltas, nearest-neighbour ties,
median/95th-percentile/maximum absolute changes and ranked TSV discrepancies.
Do not infer transmission or tree-topology equivalence from these counts.

## Remaining gates before a replacement recommendation

- Run all native arms against the same controls. iVar 1.4.4 is installed in the
  local `skrotis` environment; Hopper provisioning remains a separate step.
- Add equivalent shifted/padded and overlapping variant records, terminal
  absence, true stop/frameshift and mutation-linkage controls to the truth panel.
- Run raw synthetic reads through the real mapper; explicit-CIGAR probes do not
  establish mapper robustness. Check the two reported homopolymer loci directly.
- Freeze runnable minority-call preprocessing using documented tool operations;
  neither current native BCFtools arm supplies the full desired policy.
- Prepare positional coverage masks/trimming for actual candidate sequences;
  run operational and pooled distance comparisons on the cohort.
- Complete fresh Hopper runs for viable candidates and a baseline via site
  wrappers, Nextflow 26.04.3, isolated work/output paths and disabled import markers.
  Push committed code to the Hopper bare repository before using wrappers.
- Review every changed RAS/drug outcome, lost assessment, new integrity warning
  and nearest-neighbour change; explain larger matrix changes with callable bases.
- Produce full matrix TSV/JSON, difference visualizations and a recommendation.
  No IQ-TREE topology comparison or clinical-validation claim is included here.

Do not replace archive provenance with converter-4 metadata for another tool.
Any future production interface changes require explicit coordinate mappings and
new matching reference/alignment identities. Do not suppress true frameshifts to
make results look cleaner.

## Initial verification, 2026-11-13

- The pinned BCFtools 1.21/SAMtools 1.21 images and local iVar 1.4.4 completed
  all 36 explicit-CIGAR probes. The production codon extractor recovered all
  known molecule counts, including synonymous aggregation. This is not a claim
  that every native consensus met the desired resistance policy.
- At 85/15, iVar 0.85 omitted the minority; at 86/14, iVar 0.90 included it.
  Haploid BCFtools consensus omitted the minority, whereas unconditional ALT
  application selected it even below threshold. These are capability findings,
  not proposed production settings.
- The synonymous 84/8/8 control yielded 16% alanine by observed codon counts.
  Nucleotide-IUPAC translation also allowed unobserved serine and valine. The
  two-codon linkage control admitted unobserved leucine through IUPAC expansion.
- An 8% single-base indel caused the raw-VCF AWK path to mask eleven positions;
  normalization reduced the mask to one or two, but did not establish the
  desired frequency policy. Native majority outputs retained supported sequence.
- Five benchmark regression tests, ten existing consensus-policy tests and the
  fourteen Virtitta distance tests passed. The distance harness also completed
  real MAFFT operational, pooled and common-callable comparisons on the generated
  controls, with separate reports for unsupported length-changing coverage maps.
- All 240 archived analyses were inventoried and hashed. The read-evidence scan
  reached all 240: 207 had some codon-count difference between scaffold and
  IUPAC remapping (4,661 sites), but only two adequately covered sites crossed
  the >=15% amino-acid threshold. These counts are an independent diagnostic,
  not a fresh 240-sample pipeline rerun or proof of changed clinical predictions.
  Final source-integrity verification is represented by the COMPLETE file below.

Private output roots for this session: `/tmp/consensus-benchmark-all-02`,
`/tmp/consensus-pinned-probe-distances-01` and
`/tmp/consensus-benchmark-cohort-01`. The last directory is complete only when
its `COMPLETE` file exists. Earlier `*-01` probe directories include intermediate
harness versions; use `all-02` for the raw-versus-normalized baseline distinction.
The campaign CSV is `/fs1/jonas/hcv/bjorn-rerun/240samples.csv` on Hopper
(`/mnt/fs1` on the laptop). Its 240 identities exactly match the frozen archives
and its sample IDs are unique, which is required by the current workflow joins.
The original archive scan completed source-integrity verification successfully.
Matplotlib 3.11.2 is now available locally; set `MPLCONFIGDIR` to a writable
scratch directory and pass `--plots` for standalone difference heatmaps.

## Fresh Hopper baseline

`scripts/runme_consensus_benchmark.sh` is a dedicated site wrapper using
Java 23.0.2, Nextflow 26.04.3, strict syntax and the existing HPC/Slurm/Apptainer
profiles. Submit it from a pinned checkout after publishing the commit to the
Hopper bare repository:

```bash
sbatch scripts/runme_consensus_benchmark.sh /fs1/jonas/hcv/bjorn-rerun/240samples.csv /fs1/jonas/virpipa/consensus_benchmark_UNIQUE
```

It requires a new output directory, records the source commit, copies and hashes
the CSV, validates all 240 rows, hashes the unique FASTQs with four bounded
workers, and retains debug intermediates for upstream comparisons. Import
markers and completion-log publication are disabled explicitly. The source
FASTQs are read, never restored or modified by this wrapper. The wrapper does
not switch the production converter or implement a new resistance policy.

The verified BioContainers iVar image tag is
`quay.io/biocontainers/ivar:1.4.4--h077b44d_0`, registry manifest digest
`sha256:1de5c4e08c00c037797137b903f0849b3a4d55d750a6e3a4efeda0fcfe795cb7`.
Record the actual SIF checksum separately after transfer to Hopper. The user is
provisioning that image; local iVar testing already uses version 1.4.4.
