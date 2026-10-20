# Deferred BLAST subtyping implementation plan

Status: deferred until resistance calling has been checked on representative data.
This document records the agreed direction, not an implemented classifier.

## Goal and boundaries

Replace opaque first-row reporting with reviewable subtype evidence, without a
large default table or an unvalidated confidence score. Keep analytical 15% IUPAC
as the intended primary input. BLAST treats ambiguous query characters as mismatches,
so preserving mixtures in FASTA does not make BLAST mixture-aware. Majority can be
a diagnostic comparison, not an assumed biologically superior replacement.

Do not implement any BLAST behaviour or presentation changes as part of resistance
completion. The current scaffold/report versus IUPAC/resistance typing inconsistency
remains a known limitation until this work is implemented. Do not silently change
subtype-specific resistance rules in Virtitta or infer strain linkage from IUPAC.

## Pipeline implementation sequence

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
