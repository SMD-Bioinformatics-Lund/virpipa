# Resistance reference audit

Audit date: 2026-10-06. Database: `refgenomes/hcvglue/hcvgluerefs`, built 2024-10-23. The audit inspected the installed database rather than assuming it matched the current upstream HCV-GLUE project.

## Installed database

- 213 entries, 2,160,033 bases, with 93 distinct header labels. These include genotype-only labels and uncertain labels such as `4g?` and `6e?`; they are not 93 confidently assigned subtypes.
- Every entry is exactly 10,141 bases long. The export contains 153,580 `N` characters and other IUPAC ambiguity codes.
- Removing all `N` characters from `1a_AF009606` exactly recovers the repository's 9,646-base AF009606.1 genome. This supports the interpretation that at least some `N` characters represent alignment gaps/padding, rather than unknown biological bases.
- AF009606 is present; AF011751 is absent. Resistance numbering now uses AF011751.1 from a separate FASTA.

The database can support the existing BLAST subtype search, but its padded sequence export should not be used directly for read mapping or codon translation. Globally deleting `N` characters is also unsafe: the export alone does not distinguish alignment gaps from genuinely unknown or missing sequence. A resistance panel should use the original accession sequences and a verified alignment to H77.

## References for the current additional subtypes

All ten requested subtypes are represented. Gene intervals were located through the padded H77 entry: NS3 columns 3707–5599, NS5A 6545–8095, NS5B 8096–9868 (1-based inclusive). Counts below describe the aligned export, not validated GenBank CDS annotations.

| Subtype | Accessions in database | Observations in resistance gene intervals |
| --- | --- | --- |
| 2c | D50409, JX227949 | D50409 has 1773 non-N NS5B bases; JX227949 has 1681. No other ambiguity codes in these gene intervals. |
| 2i | DQ155561 | NS5B has 1734 non-N bases; its H77-aligned interval contains 39 Ns. No other ambiguity codes in these gene intervals. |
| 2k | AB031663, JX227953 | Both have 1773 non-N NS5B bases and no other ambiguity codes in these gene intervals. |
| 3h | JF735126, JF735121 | Both have 1773 non-N NS5B bases. JF735121 has ambiguous bases in NS3/NS5A; JF735126 in NS3/NS5B. |
| 4c | FJ462436 | NS3 has three ambiguous bases. NS5B has 1773 non-N bases. |
| 4d | DQ418786, FJ462437, EU392172 | DQ418786 has no non-N ambiguity codes in these gene intervals. FJ462437 has only 1503 non-N NS5B bases. EU392172 has ambiguities in NS5A/NS5B. |
| 4k | EU392173, FJ462438, EU392171 | All have 1773 non-N NS5B bases. All have some non-N ambiguity codes; EU392173/EU392171 fewer than FJ462438. |
| 4m | FJ462433, JX227972 | Both have 1773 non-N NS5B bases. JX227972 has no non-N ambiguity codes in these gene intervals; FJ462433 has ambiguities in all three. |
| 4r | FJ462439, JX227976 | Both have 1773 non-N NS5B bases. FJ462439 has ambiguities in NS3/NS5A; JX227976 in NS5B. |
| 6e | DQ314805, EU246932 | Both have 1773 non-N NS5B bases and ambiguities in NS3/NS5A; DQ314805 also in NS5B. |

All these entries have 1893 non-N bases in NS3. Their NS5A non-N lengths vary from 1332 to 1410 bases. Such differences may represent real subtype indels; they require alignment and original sequence verification rather than fixed H77 offsets.

Extraction command (uses the existing local container):

```bash
apptainer exec -B /home,/tmp assets/containers/blast_2.16.0.sif \
  blastdbcmd -db refgenomes/hcvglue/hcvgluerefs -entry all \
  -outfmt '%f' -out /tmp/hcvgluerefs-audit.fasta
```

Database sequence-volume SHA-256 (`hcvgluerefs.nsq`): `9838ab27f9226aa3e3db57e6f57dbfa74fdeb24a870a2d659a86778e4806babb`.

## Subtype handling to verify

[The geno2pheno methods](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155869), sections 2.1.1–2.1.3, document separate subtype assignment for NS3, NS5A, and NS5B, a default reporting reference for each region/subtype, H77 position numbering, and reporting matching rules from the closest common subtype for rare subtypes. An example uses related subtype 4d for a 4k query. That example is not a universal lookup table.

The committed 2026 rules snapshot includes genotype-wide selectors and explicit lists of rarer subtypes. Therefore rare-subtype handling must respect each rule's actual selector: an applicable genotype-wide rule is different from extrapolating a rule restricted to another subtype. The older paper alone does not establish how today's web service combines these cases.

Recommended next step: compare anonymized online reports for the ten subtypes against the local region sequences. Record the input checksum, web-service/rules version when available, subtype and selected reference accession separately for each region, reference used for substitution labels, scored substitutions, and drug outcomes. Extract only the matching NS3/NS5A/NS5B sequences from verified public accessions for a versioned resistance panel. Include validation cases with mixtures, missing sequence, indels, and discordant region subtypes. Keep any re-identification key outside Git.

The main workflow and historical re-call currently take one subtype from the first 15% whole-consensus BLAST data hit. The historical re-call bug that took the assembly-reference label from the report has been corrected. Neither entry point implements region-specific assignment or the related-subtype procedure yet.
