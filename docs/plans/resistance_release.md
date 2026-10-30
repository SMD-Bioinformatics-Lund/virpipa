# Resistance release checklist

## Agreed scope

Complete BLAST ranking and shared analytical-15%-IUPAC subtype selection for QC,
Virtitta and resistance, then validate and release the combined resistance redesign.
Preserve raw searches and mapping-reference metadata. Do not change consensus policy,
assembly, polishing, mapping, resistance rules, H77 AF011751.1 or affine alignment scores.

## Implementation checkpoints

- Ranking helper, synthetic tests and installed-database audit: implemented. Suggest
  committing this isolated checkpoint before wiring it into the workflow.
- Implemented: retrieve all database subjects with explicit search limits and exact
  identical-base counts; publish versioned typing evidence and database/input checksums.
- Implemented: feed one authoritative assignment to resistance and QC; use the same
  ranking in re-call. Keep samples with ambiguous/unassigned typing, with explicit
  not-assessed drug results rather than arbitrary subtype selection. Genotype-only
  assignments may evaluate only explicitly genotype-applicable rules.
- Implemented QC evidence without changing historical field meanings. Pending: Virtitta must
  prefer the new assignment, including explicit null assignments, over legacy fields;
  preserve old-archive import behaviour and app-owned review state.
- Pending: compact subtype/support summary and expandable alternatives (up to five
  subtype groups, three references per group), with identity, aligned bases, coverage
  and descriptive bit scores. Secondary hits do not establish mixed infection.

## Release gates

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
