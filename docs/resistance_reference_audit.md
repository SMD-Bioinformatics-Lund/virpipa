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

## Observed online results

Fourteen supplied geno2pheno PDF reports dated 2026-10-06 cover all ten additional subtypes. Their title states version 1.0, but they do not identify a resistance-rule snapshot or reference-panel revision. Original identifiers and PDFs were not added to Git. Subsequently supplied sequences have been anonymized in `assets/test_data/resistance/rare_cases/`. The temporary re-identification key was lost during a WSL restart and recovered from restored original sequences on 2026-10-07; it is now stored persistently outside Git. The fixture README describes checksum-based recovery without guessing identifiers.

All 42 region assignments and all 224 drug-table rows were extracted. Within every report, NS3, NS5A, and NS5B agree on subtype and use the same accession for the subtype-relative substitution list. All observed accessions are present in the installed BLAST database.

| Subtype | Reports | Reference displayed for all three regions |
| --- | --- | --- |
| 2c | 1 | D50409 |
| 2i | 1 | DQ155561 |
| 2k | 1 | AB031663 |
| 3h | 1 | JF735121 |
| 4c | 1 | FJ462436 |
| 4d | 4 | DQ418786 |
| 4k | 1 | EU392173 |
| 4m | 1 | FJ462433 |
| 4r | 2 | FJ462439 |
| 6e | 1 | DQ314805 |

These are observed reporting references, not proof that every possible query of a subtype selects the same accession, nor a table translating rare subtypes into common subtypes for rule evaluation.

The subsequently supplied combined CSV contains 15 rows and 95 columns. All 14 PDF cases match the CSV on regional subtype, similarity, reporting reference, and all 224 drug predictions/scored-mutation lists. PDF `susceptible*` is explicitly exported as `substitution on scored position`. The CSV also supplies per-drug fold and evidence fields; for example, the 4m daclatasvir and ledipasvir rows contain fold `150` and evidence code `CR`. These fields should retain their exported meanings until the service's definitions are verified.

One additional CSV-only case is labelled `1*2` in all three regions, with reporting reference AJ851228. The user confirms that this sample was assigned genotype 1 without a letter subtype. Record genotype `1`, subtype unassigned, and the raw service label `1*2` separately. Its drug results are susceptible or substitution on scored position. The precise meaning of the service's `*2` suffix has not been established; do not convert it to 1a or 1b for rule selection.

The regional `all_resistance_muts_*` fields do not include every per-drug scored mutation. In the 3h case, the NS5A aggregate contains `93H`, while the pibrentasvir scored list contains `28I,93H`. Future comparison fixtures and IGV scored-position annotations should therefore use the per-drug lists, not only the regional aggregates. The identifier-bearing CSV remains outside Git.

The reports give direct resistant predictions for some rare subtypes:

- 3h: NS5A 93H is scored as resistant for daclatasvir, elbasvir, ledipasvir, and velpatasvir. Pibrentasvir is susceptible* with scored amino acids 28I and 93H; the subtype-relative list contains the 28IM mixture.
- 4m: daclatasvir and ledipasvir are resistant with scored amino acids 30S and 31L; elbasvir has reduced susceptibility with 30S; ombitasvir is susceptible* with 31L.
- 4r: one report is resistant to several NS5A inhibitors with 93H, while the other has susceptible* calls for daclatasvir (30R) and ombitasvir (28M).

The genotype-wide rules in the committed CSV support direct resistant predictions for these 93H and 30S observations. Automatically downgrading all rare-subtype hits to a related-subtype warning would contradict the supplied web results. The PDFs do not reveal the website's internal common-subtype selection procedure.

The susceptible* observations expose a limitation of comparing only against H77:

- 4m NS5A 31L matches H77 L31 but differs from the displayed subtype reference (M31L); the website scores 31L for ombitasvir as susceptible*.
- 4r NS5A 28M matches H77 M28 but differs from the displayed subtype reference (I28M); the website scores 28M for ombitasvir as susceptible*.
- Conversely, NS5A 30R appears in the H77-relative list but not the subtype-relative list for the supplied 4c, 4k, and 4r reports, yet is scored for daclatasvir as susceptible*. The four 4d reports also have H77-relative Q30R, but daclatasvir is susceptible with no scored mutations.

Thus neither a universal H77 difference check nor a universal difference check against the displayed subtype reference reproduces the website's scored-mutation presentation. The subtype/drug scoring baseline needs further verification; these observations should become regression cases once the corresponding input sequences are supplied.

The drug-table coverage issue has been corrected: the caller now constructs outcomes for all drugs in the catalogue and determines licensing from the full catalogue's genotype metadata, separately from subtype-specific mutation rules. For 2c, for example, the committed rules select no glecaprevir, pibrentasvir, sofosbuvir, or voxilaprevir rows, while the online report lists all four as susceptible. Without applicable rules, susceptibility requires assessed catalogue positions and does not extrapolate rules from another subtype; missing positions produce not assessed.

The original public accessions have now been retrieved and verified. `refgenomes/resistance/panel.json` stages one reporting genome per observed assigned subtype, with NS3/NS5A/NS5B gene coordinates and linking alignments to AF011751.1. The eleven supplied sample genomes each match one CSV row and cover each additional subtype plus the unassigned genotype-1 case once. Their anonymized expectations retain per-drug predictions and scored amino acids.

Next steps: validate subtype selection and scored-mutation presentation against these cases, resolve the susceptible* baseline, and expand the candidate panel's subtype coverage before activating it in the workflow. Include mixtures, missing sequence, indels, and discordant region subtypes in validation. Keep any re-identification key outside Git.

After local VADR annotation and boundary review on 2026-10-07, regression fixtures retain the original resistance-gene coordinates, including three non-triplet intervals with documented caveats. All 176 expected drug rows are now emitted. All seven resistant/reduced-susceptibility predictions match; the 23 remaining outcome differences are susceptible versus susceptible*. This does not establish per-drug scored-mutation-list parity. Five supplied cases are user-identified mixtures and are retained as IUPAC comparison cases, not phased strain-level truth.

Deferred pipeline work includes coding-integrity warnings (definite/possible stops, non-triplet indels and uncertain boundaries), site-level mapping reliability, and target-gene nucleotide/codon/amino-acid alignment artifacts. Mixed-population and compound-mutation linkage limitations must remain separate from integrity warnings. Virtitta presentation dependencies are tracked in its `docs/PLANS.md`.

The main workflow and historical re-call currently take one subtype from the first 15% whole-consensus BLAST data hit. The historical re-call bug that took the assembly-reference label from the report has been corrected. Neither entry point implements region-specific assignment or the related-subtype procedure yet.
