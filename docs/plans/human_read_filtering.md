# Future human-read filtering replacement

The user intends to retire Hostile in favour of Deacon. Migration timing and the
implementation are not yet decided. Current work remains the bounded Hostile
ordering experiment; no Deacon installation or pipeline substitution is included.

Evaluate the replacement against the reproducibility findings in
[the scaffold investigation](../scaffold_reproducibility.md):

- Preserve paired-read identity, mate synchronization, quality scores and relevant
  header/UMI information. Explain intended changes in read retention rather than
  assuming identical filtering decisions between different tools.
- Verify deterministic read order, or explicitly address the downstream sampler's
  dependence on order. A fixed random seed alone is insufficient.
- Test the complete filter → seeded subsampling → assembly/polishing chain on the
  same frozen inputs, including mixtures, low coverage and below-cutoff controls.
- Establish offline Hopper containers/reference staging and record tool/reference
  provenance. Review performance with the site's actual CPU and memory allocations.
- Preserve the published QC/import interface through a deliberate adapter or document
  and validate an additive contract change before replacing the existing producer.

Keep diagnostic data and case-identification keys outside Git. Treat replacement of
the human-read filter separately from the resistance caller and Pilon cycle policy.
