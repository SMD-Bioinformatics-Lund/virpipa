# Consensus benchmark results

Evaluation on `consensus-benchmark`. The pipeline's production converter,
resistance caller, output schemas and result archives are unchanged. Local
cohort distance analysis is still running; this report is not a release decision.

## Findings established so far

The reported opposing indels are reproducible alignment artifacts in the masked
IUPAC reference. Standard BCFtools and iVar reconstruction removes those patterns
in both reported windows when the same eligible read pool is remapped. This does
not establish that every apparent indel is artificial or may be repaired.

A FASTA consensus is insufficient for the requested >=15% amino-acid policy.
Read codons retain linkage and synonymous changes that IUPAC expansion loses.
BCFtools haploid genotypes and iVar cumulative thresholds also do not implement
an interchangeable per-amino-acid 15% minority policy.

The count-first design has an established reference: MiCall documents nucleotide
and amino-acid counts separately, with its HCV resistance step consuming
`amino.csv`. That is a useful independent target for validating a production
read-codon implementation; MiCall itself was not benchmarked here.
[MiCall pipeline and file documentation](https://cfe-lab.github.io/MiCall/steps.html).
BCFtools 1.21 also explicitly documents that consensus does not use AD fields to
choose alleles by depth.
[BCFtools 1.21 manual source](https://raw.githubusercontent.com/samtools/bcftools/1.21/doc/bcftools.txt).

## Design and provenance

- 240 analyses from 40 runs; isolated fresh baseline followed by nine native arms.
- Nextflow 26.04.3, BCFtools/SAMtools 1.21, iVar 1.4.4, Sentieon 202308.03.
  Local truth-panel and Hopper cohort container checksums match for all five images.
- Native arms: scaffold control, current IUPAC and analytical-majority controls,
  BCFtools haploid genotype consensus, iVar majority/0.85/0.90, and AWK
  IUPAC/majority controls on the same BAQ-disabled calls.
- One frozen pool of mapped primary, QC-passing, nonduplicate UMI-consensus reads
  supplies both callers and all native remappings. UMI consensus is not repeated.
- Complete read codons require base quality >=13. Synonymous codons are summed
  before the integer test `count * 100 >= depth * 15`. Counts remain available
  below the assessment gate. Seven-complete-codon and one-complete-codon gates
  are reported separately; choosing a production depth policy remains separate.
- Typing and resistance catalogue are held constant. Existing H77 numbering,
  codon extraction, integrity checks and rule evaluation are reused. Archived
  VADR boundaries are projected conservatively rather than freshly annotated.

The November 11 archive omitted Hostile `--reorder` in all 240 actual task
scripts. The fresh baseline includes the fix committed on November 6. In the
analysis with changed archived-versus-fresh drug outcomes, filtered paired records
are identical as a multiset but reordered, replacing 522,139 of one million
polishing pairs after seeded subsampling. Thus archive-versus-fresh differences
are confounded by this upstream change, not evidence of a new converter effect
or an identical-input reproducibility failure.

Both real workflows completed successfully. Baseline validation covered all 240
identities/layouts, consensus sidecars and FASTA/CRAM identities; all 9,400 traced
tasks completed. The native workflow completed 240 tasks without failures.
Container manifests, frozen source-input hashes and 2,160 native FASTA/CRAM
reference identities were checked locally.

## Known-molecule controls

All 41 fixtures completed calling, reconstruction and local remapping through
all nine arms: **369/369** recovered the expected complete-codon counts and
retained every template. Linked and unlinked mixtures retained their distinct
known haplotypes. Controls cover 14/15/16%, synonymous aggregation, unobserved
IUPAC combinations, multiallelic calls, missing/low depth, poor base quality,
mate agreement/conflict, partial codons, terminal absence, stop codons and
1/2/3-base insertions/deletions at 8/15/50/90% support.

This validates extraction/remapping for the defined molecules, not a claim that
all native FASTAs satisfy the desired minority policy. At 85/15, iVar 0.85 omits
the minority; iVar 0.90 includes a 14% minority. Haploid BCFtools omits minorities;
unconditional ALT application can select them below threshold. Synonymous
84/8/8 mixtures yield 16% amino-acid support despite each alternate codon being
below 15%. IUPAC combinations can include amino acids absent from the molecules.

BCFtools normalization/consensus produced the same known 239-base deletion
sequence from three equivalent shifted/padded records. A deliberately conflicting
SNP inside that deletion generated an overlap warning and was skipped. That
inconsistent callset requires review, not silent record merging.

In the reported windows, remapped current IUPAC references had 137 opposing
9-base-gap records and 47 opposing 10-base-gap records, respectively. Neither
pattern occurred with the BCFtools/iVar reconstructions. The corresponding
qualifying read amino acids and drug outcomes agree across the pilot arms.

## Resistance comparison

Against the remapped current-IUPAC control, using the same observed read-codon
policy and the seven-complete-codon gate:

| Native arm | Changed assessed amino acids | Lost/gained assessed sites | Changed drug outcomes |
| --- | ---: | ---: | ---: |
| BCFtools genotype | 0 | 23 / 1 | 1 |
| iVar majority | 0 | 1 / 0 | 0 |
| iVar 0.85 | 0 | 1 / 0 | 0 |
| iVar 0.90 | 1 | 1 / 0 | 0 |
| AWK, BAQ disabled, IUPAC | 1 | 0 / 1 | 0 |

The 23 lost BCFtools assessments all belong to one NS5B interval with an internal
stop; the scaffold control has the same exclusion. Sofosbuvir changes from
susceptible to not assessed under this conservative adapter. The iVar exclusion
is in a separate low-depth interval with a non-triplet change and downstream
integrity warnings. These require review, not removal of integrity safeguards.
The adequately covered changed amino-acid set is a near-threshold NS5A 28 mixture;
its drug outcomes remain unchanged. Deletion frequencies are not quantified by
this adapter, and distant compound mutations are not phased.

Comparing the production FASTA interpretation with read-codon interpretation also
changes read pool/remapping and depth policy. At seven complete codons, 17
analyses have 76 changed drug outcomes and 296 lost assessments. These must not
be attributed to replacing AWK or to amino-acid phasing alone. At one complete
codon, matching production's permissive depth floor more closely, the comparison
has 20 changed assessed amino-acid sets and seven changed drug outcomes in five
analyses, with no lost or gained assessments. Most of those analyses have
low-depth evidence; this sensitivity analysis is not a depth recommendation.

One well-covered example is decisive independently of the read-pool change:
IUPAC `WST` permits C/S/T, but the fresh production alignment records 220 `TCT`
and 108 `AGT` serine codons plus seven `ACT` threonine codons. Observed serine is
328/335 (97.9%), threonine 7/335 (2.1%), and cysteine zero. Only serine meets 15%.
The archived run has the same IUPAC codon and only 5/320 (1.6%) threonine. Thus
FASTA expansion includes a resistance-associated amino acid below the specified
threshold and another amino acid not observed in complete codons.

An early experimental qualifier overwrote some upstream exclusion reasons with
an insufficient-depth label. The qualifier now preserves integrity reasons.
Local requalification restored 943 reason labels while asserting unchanged
amino-acid calls and assessment states for every site. Original native scorecards
are retained; the qualified local reports are authoritative for reason labels.

## Distances

Pending completion of local alignments and matrix comparison. The primary
reference is `baseline_remap`, holding the eligible read pool constant. A separate
published-baseline matrix captures current operational coverage behavior.

Prepared sequences/masks use the selected Virtitta configuration: 50-base 5-prime
trim and its poly-T handling. Each candidate has coverage in its own reference
coordinates. The actual Virtitta distance engine is pinned by commit/source hash.
Operational majority matrices reuse the matching IUPAC alignment; pooled matrices
and common-callable nucleotide comparisons distinguish alignment/missingness
changes from substitutions. Raw indel-event counts and compared-base denominators
remain available. No tree-topology or transmission inference is made.

Sixteen native iVar IUPAC/majority pairs have different lengths. Affected
operational majority comparisons are unsupported by the present Virtitta
interface, and are flagged rather than forced into one coordinate system. Their
pooled comparisons use each sequence's actual reference/coverage coordinates.

## Artifacts and remaining decision

Private inputs, identities, matrices and figures are under
`test_output_consensus_benchmark/session_20261116/`; they are ignored by Git.
`CHECKPOINT.json` identifies current output roots. Key completed artifacts are
`cohort-scorecard-qualified/`, `depth-sensitivity.json`,
`mapper-truth-validation.json`, `record-probes-02/` and
`amino-acid-evidence.png`. Full raw results/work/cache remain on Hopper in the
isolated benchmark campaign directories.

The evidence supports separating genome consensus from observed amino-acid
prevalence. It does not justify selecting native iVar thresholds or a haploid
BCFtools consensus as a drop-in 15% resistance caller. A production implementation
still needs an explicit depth policy, deletion-frequency handling and validated
annotation/coordinate behavior for edited genomes. No production switch or
automatic re-import is part of this evaluation.
