# BLAST subtyping implementation and validation

Status: reference ranking, workflow/resistance/QC integration and Virtitta import/display
are implemented and locally tested. Fresh Hopper/import release validation remains
pending. See [the release checklist](resistance_release.md).

## Goal and boundaries

Replace opaque first-row reporting with reviewable subtype evidence, without a
large default table or an unvalidated confidence score. Keep analytical 15% IUPAC
as the intended primary input. BLAST treats ambiguous query characters as mismatches,
so preserving mixtures in FASTA does not make BLAST mixture-aware. Majority can be
a diagnostic comparison, not an assumed biologically superior replacement.

BLAST behaviour and presentation changes are now the final implementation work for
the resistance release. The current scaffold/report versus IUPAC/resistance typing
inconsistency is resolved for new archives through the shared analytical assignment;
the historical scaffold fields retain their original meaning. Do not silently change
subtype-specific resistance rules in Virtitta or infer strain linkage from IUPAC.

## Pipeline implementation sequence

### Ranking checkpoint

`scripts/rank_blast.py` reads standard 12-column BLAST output, optionally with
`nident` appended as column 13. It selects the maximum summed-bit-score chain per
reference using strictly non-overlapping query and subject intervals, consistent
orientation and subject order. This conservative policy discards overlapping segments
rather than attempting to prorate their scores. It never pools separate references.
Equal-score chains are selected deterministically; equal-score winning subtype
groups remain ambiguous. Unknown top-scoring labels prevent a confident assignment.
Genotype-only labels are retained explicitly. No usable hits yield unassigned typing;
malformed rows fail rather than becoming a no-hit result.

Identity uses total identical bases divided by total alignment columns, including
gaps. Legacy 12-column inputs approximate identical counts from rounded BLAST identity;
new searches will append exact `nident`. Whole-query coverage counts selected query
intervals once; observed coverage counts non-N characters in those intervals divided
by all non-N query characters. IUPAC ambiguity characters other than N remain observed.
All-N queries have null observed coverage. These metrics are not confidence values.

The existing subtype fixture retains 3a / D17763: three selected segments have summed
bit score 13888 and 99.7673% whole-query coverage. This is a fixture check, not a new
real-sample or biological validation result.

The workflow now calls the helper after analytical IUPAC BLAST and supplies the complete
database-file set. The same result feeds resistance and QC, and is published as typing
JSON. Sequence checksums bind the declared result to its analytical input. Re-call checks
archived evidence against the raw search; legacy searches without evidence retain empty
database provenance and are not assumed to contain every possible alternative.

### Implementation requirements (completed; fresh-run validation pending)

1. Audit reference labels, duplicates, incomplete/padded sequences and unequal subtype
   representation. Database absence of a runner-up is not proof of confidence.
2. Specify and test reference-level ranking of multiple local alignment segments:
   compatible ordering/orientation, overlap handling and non-duplicated query coverage.
   Do not assume the first segment is the complete reference match. Do not sum scores
   across different references belonging to a subtype.
3. Make coverage denominators explicit: distinguish whole-query coverage from the
   proportion of observed/non-N sequence represented. Keep standard BLAST identity
   distinct from any future IUPAC-compatible similarity metric.
4. Publish the evidence through additive QC JSON fields. Use a shared typing result
   for reporting and resistance input once the new ranking is validated. Preserve
   raw BLAST files and record input role, reference database and ranking method.
5. Document ties, absent hits, genotype-only labels and unsupported subtypes before
   replacing the existing first-row caller. Score-gap thresholds and no-call policy
   require a separate reviewed decision, not a hidden implementation default.

## Virtitta presentation

- Future extension: classify subtype 1a into clades I/II using a validated method
  and reference set, with method/reference provenance and an explicit unassigned
  result when unsupported. Present this below subtype, not as a resistance class;
  do not infer clades from a single BLAST accession or change resistance rules.
- Keep the main summary compact: closest subtype, identity and sequence support.
- In sample details, show up to five strongest subtype groups, represented by their
  best reference; expand each to up to three available reference hits.
- Show accession, identity, aligned bases, query coverage and bit score. Group by
  subtype, not only genotype, so several references of one subtype cannot obscure
  a competing subtype. Counts are display limits, not evidence requirements.
- Keep score margins in technical details, labelled descriptive rather than a
  confidence probability. Do not label a secondary hit as a mixed infection.
- Consider a later query-location strip to distinguish overlapping competitors from
  different subtype matches in different regions. Regional typing, mixture resolution
  and phylogenetic confirmation are separate extensions, not implied by this table.

## Validation and design gate

Before coding the classifier, settle segment-ranking, coverage and tie/no-call
semantics against fixtures. Include split alignments around Ns, overlapping segments,
duplicate references, tied scores, few/no hits, unknown subtype labels, partial
genomes and mixed samples. Compare IUPAC and majority on the same samples without
equating a higher score with greater biological accuracy.

Use anonymized comparisons with known subtype evidence, checking both the displayed
alternatives and downstream resistance-rule selection. Record reference-set limitations
and unresolved cases. Do not add arbitrary bit-score confidence thresholds.

Precedents for the design discussion:

- [NCBI ambiguity-code behaviour](https://blast.ncbi.nlm.nih.gov/blast/BLAST_guide.pdf)
- [NCBI BLAST report design](https://pmc.ncbi.nlm.nih.gov/articles/PMC3692093/)
- [Geno2pheno HCV methods and alternative-reference display](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155869)
