# Resistance release checklist

## Agreed scope

Complete BLAST ranking and shared analytical-15%-IUPAC subtype selection for QC,
Virtitta and resistance, then validate and release the combined resistance redesign.
Preserve raw searches and mapping-reference metadata. Do not change consensus policy,
assembly, polishing, mapping, resistance rules, H77 AF011751.1 or affine alignment scores.

## Implementation checkpoints

- Ranking helper, synthetic tests and installed-database audit: implemented and committed.
- Implemented: retrieve all database subjects with explicit search limits and exact
  identical-base counts; publish versioned typing evidence and database/input checksums.
- Implemented: feed one authoritative assignment to resistance and QC; use the same
  ranking in re-call. Keep samples with ambiguous/unassigned typing, with explicit
  not-assessed drug results rather than arbitrary subtype selection. Genotype-only
  assignments may evaluate only explicitly genotype-applicable rules.
- Implemented QC evidence without changing historical field meanings. Virtitta prefers
  the new assignment, including explicit null assignments, over legacy fields;
  legacy import and app-owned review state are covered by regression tests.
- Implemented: compact subtype/support summary and expandable alternatives (up to five
  subtype groups, three references per group), with identity, aligned bases, coverage
  and descriptive bit scores. Secondary hits do not establish mixed infection.

## Release gates

- Reproducibility investigation found reordered Hostile output feeding position-based
  seeded subsampling, different retained read subsets and separate iterative Pilon
  cycling. Fixed-input Pilon was repeatable in the tested cases. See
  [the investigation](../scaffold_reproducibility.md). The tested native Hostile
  `--reorder` option is enabled in pipeline code; pinned aaa2633 passed two fresh
  full-sample analyses of case A and all 38 semantic comparisons of the complete
  chain. Accepted release decision: retain the existing Pilon stopping behaviour
  (equal consecutive change counts or the configured round limit, default ten).
  Cycle detection/reporting and early cycle termination are deferred; early breaking
  could select a different final sequence and requires a separate reviewed policy.
- Debug FASTQ publication collision corrected: subsampling limits now publish into
  separate `fastq/subsample_<requested_pairs>/` directories. The tiny two-alias fixture
  reproduced the old overwrite and passed on Nextflow 26.04.3 with both 25- and
  50-pair sets intact, matching mates and unchanged read payloads. Routine mode
  published neither set; strict Nextflow lint passed. This changes publication only;
  the completed full-sample analytical validation remains applicable.

- Local ranking, uncertainty, workflow, publication, import and display regressions;
  strict Nextflow validation and container-backed module checks.
- Fresh representative Hopper batch through the site launcher in isolated dev paths;
  explain subtype changes and every resulting resistance change.
- Validate desktop/webIGV, clustering, distances and repeat import. User has confirmed
  comments and QC status survive repeat import; additionally verify category, group and
  override preservation before collection-wide refresh.
- Review changes and suggest separate commits per checkpoint/repository. Human-controlled
  tags, publication and deployment follow validation; do not update a running checkout.
- Reprocess the collection to separate outputs and retain original archives until
  acceptance. An unexplained wrong resistance call, shifted coordinate or false clear
  assessment blocks release.

## Explicit deferrals

Accept documented susceptible versus susceptible* baseline/presentation differences;
do not claim complete online parity. Defer independent NS3/NS5A/NS5B typing, clades,
reference expansion and affine-score calibration. Interactive alignments, exploratory
cutoffs, curator overrides, richer gene-wide IGV annotation and coverage-BED retirement
are future work, not release prerequisites. Existing alignment artifacts can support
later UI work; do not promise all future features need no regenerated evidence.
