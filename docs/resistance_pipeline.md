# HCV resistance calling

VirPipa calls resistance from the 15% IUPAC consensus. The VCF is not used for resistance calling: it is relative to the sample's polished consensus and therefore cannot represent majority substitutions reliably.

`scripts/annotate_vcf_resistance.py` evaluates the committed geno2pheno rules in H77 amino-acid numbering (the pinned AF011751.1 sequence). It uses ambiguity-aware nucleotide alignment to an exact-subtype bridge from `refgenomes/resistance/panel.json` when available, then projects through its verified H77 linking alignment. Otherwise it records a direct H77 fallback. No related-subtype rule extrapolation is performed. Common-subtype bridges verified against the installed HCV-GLUE alignment include 1b, 2a, 2b, 3a, 5a and 6a; their equivalence to the web service's reporting-reference selection has not been established.

Alignment uses Biopython affine-gap alignment (compatible IUPAC bases +2, mismatch -3, unknown 0, gap opening -10, extension -1), partitioned by unique exact nine-nucleotide anchors agreeing with the whole-gene alignment. Each interval retains only mappings shared by all optimal alignments, up to 256. Exceeding that bound leaves the interval and downstream positions unassessed. Codons must project completely and contiguously. Terminal absence is not a deletion; an internal codon deletion requires two reliable flanks and is represented by a one-base genomic anchor. Unresolved non-triplet indels, uncertain starting phase and definite internal stops prevent downstream calls. Valid upstream sites remain assessable. This is conservative assessment, not automatic assembly repair.

The authoritative baseline remains the 15% IUPAC consensus. At every nonzero usable depth the converter selects every canonical allele whose AD frequency is **at least** 15%, including REF only when it qualifies. Usable depth is the sum of valid canonical-base AD counts; missing FORMAT/DP does not discard valid AD, while an explicitly insufficient DP makes the record uninterpretable. Zero usable coverage becomes `N`; depths 1–6 preserve observed alleles and qualify resistance assessments with a low-depth curator-review warning, including susceptible assessments. Indels are not applied and ambiguous indel records are masked separately from zero coverage. A four-allele IUPAC `N` is also unassessed, but is not reported as missing coverage.

The pipeline publishes three coordinate-compatible sequences: `<sample>.fasta` is the reference-filled mapping/CRAM/IGV scaffold, `<sample>-0.15-iupac.fasta` is analytical IUPAC, and `<sample>-majority.fasta` is analytical majority. Both analytical outputs come from the same unfiltered all-sites VCF with minimum usable depth one. Analytical majority chooses the most abundant nucleotide; ties prefer REF only when tied for maximum, otherwise A/C/G/T order. Scaffold generation retains its seven-observation threshold and reference fallback. The new analytical majority is for export and distance comparison, not an alternative CRAM reference.

All headers contain only the identifier. `<fasta>.consensus.json` (schema `1.0`, converter `4.0`) stores identity, uppercase sequence checksum, length, mode, cutoff/depth policy, source-reference checksum, coordinate origin, usable-depth array, masking reasons and diagnostic counts. The sidecar is validated against the FASTA. Compact provenance summaries and sidecar paths are included in QC JSON; SQLite is an imported copy, not the source of truth. Codon depth qualifications use this generating evidence rather than a subsequent remapping. Existing complete-read codon and indel-evidence minimum depths remain seven.

When the matching CRAM is supplied, the caller also counts complete codons from usable primary alignments. Bases require quality 13, sites require seven complete codons, overlapping mates count once, and conflicting mate codons are discarded. These frequencies are supporting evidence for exploratory 2–50% thresholds; they do not replace the authoritative 15% result.

Codon evidence cannot cross an insertion/deletion or skip. CRAM contig identity, length and available reference MD5 must match the analytical FASTA. Net non-triplet indels are also summarized in 30-nt windows using fragments spanning both flanks, quality-filtered indel flanks and discarded mate conflicts. At least seven fragments and 15% support produce a separate read/assembly integrity finding; nearby CIGAR placements can contribute to the same window. This signal does not prove a pipeline error or separate strains. Mapping-quality policy remains unchanged (no additional MAPQ cutoff).

Rules support `and`, same-position `or` alternatives, explicit deletions, and `168any` (a non-reference amino-acid substitution, not wild-type D168). Unsupported syntax fails closed. Each component is matched, disproved or unknown. A disproved component rules out a compound rule even if another component is unknown. A full match determines the drug result; otherwise any unresolved potentially defining rule prevents a susceptible result. The special partial-compound coverage warning applies **only** to the sole applicable rule for a licensed antiviral when some required mutations are observed, none are disproved, and missing coverage prevents evaluation of the remaining components. Coding-integrity and linkage limitations are separate, not partial-rule resistance calls.

## Alignment parameter status

The affine scores above are **candidate engineering parameters**, not geno2pheno
settings or biologically calibrated HCV values. They are unchanged in this update.
Keeping them fixed isolates this cleanup from changes in coordinate-alignment behaviour;
it does not validate the gap penalties, unknown-base handling or anchoring policy.
These scores affect coordinate projection, not antiviral resistance scores directly.
Regression tests establish implementation behaviour, not biological accuracy.
Calibration and parameter-sensitivity analysis remain future work after representative
data review. The NS5A linking alignment has had 114 universally gapped columns removed;
ungapped reference sequences and every H77-to-reference coordinate map are unchanged.

Pipeline `0.3.0` identifies the analytical-consensus contract change. Start a fresh
run after updating VirPipa; resuming a run across code versions is not a supported
validation strategy. There are no new external-tool version checks or cross-version
resume checks. See [the validation checklist](resistance_validation.md) for the next
local and Hopper checks, and [the deferred BLAST plan](plans/blast_subtyping.md) for
typing work deliberately excluded from this update.

## Outputs

- `*_resistance.tsv`: resistance-associated substitutions supporting fully matched rules.
- `*_resistance.bed` and `*_resistance.gff`: compatibility tracks containing called resistance substitutions.
- `*_resistance_by_drug.tsv`: every drug in the rules catalogue with an explicit outcome (`resistant`, `reduced susceptibility`, `susceptible`, `susceptible with scored substitutions`, `not licensed`, or `not assessed`). Licensing uses the catalogue's genotype metadata independently of subtype-specific mutation rules. A licensed drug with no applicable mutation rules is susceptible only when all its catalogue scoring positions are assessed; otherwise it is not assessed. This does not extrapolate another subtype's resistance rules.
- `*_resistance.json`: schema `3.1`, explicit drug outcomes, component states/reasons, reference/rule checksums, all resistance sites, nucleotide depths and low-depth qualifications, codon frequencies, counting parameters, coding-integrity findings and per-gene nucleotide alignments/translations. `consensus_provenance_verified` confirms the converter protocol, **not** clinical validation. Only explicit `--allow-legacy-consensus` comparison runs may omit sidecars; these receive an unverified-consensus warning and cannot rebase tracks.
- `*_resistance_sites.gff3`: all resistance-rule positions for IGV. Popup attributes include H77 and sample amino-acid positions, codons, observed amino acids, drugs, rules, depth, and frequencies. `observed_aa` includes all observed mixture alleles; `resistance_aa` includes only those supporting fully matched rules.

Site tracks use concise `Name` labels and distinguish scored substitutions from resistant calls. Sites without reliable genomic coordinates remain in JSON, not fabricated GFF intervals. Gene-level integrity findings are separate GFF features. Standalone tracks use the actual analytical FASTA contig. The workflow supplies `--track-reference-fasta` to target the main scaffold for a single combined IGV view; both FASTAs must have valid sidecars, a shared `coordinate_source_sha256` and matching lengths. Scaffold generation establishes a fresh coordinate origin from its input reference sequence checksum; analytical consensus inherits that origin through the scaffold sidecar. Optional `--consensus-metadata` and `--track-reference-metadata` arguments identify explicitly staged sidecars; otherwise `<fasta>.consensus.json` is used. Suffix stripping is not a coordinate conversion. The resistance module publishes all six required files under `<outdir>/<run>/<sample>/` in debug mode.

The per-sample QC JSON carries schema version, explicit `by_drug` records (including mutations and warnings), provenance, site assessment, integrity findings and gene alignments. Virtitta displays these rather than recalculating rules. JSON download provides alignment evidence; an interactive alignment viewer and cutoff controls remain deferred.

Track coordinate validation uses shared coordinate origin and equal lengths, not per-base agreement: observed analytical substitutions can legitimately differ from the reference-filled scaffold, particularly at depth 1–6. Sequence checksums validate each sidecar independently.

The rules snapshot is `assets/hcv_geno2pheno_rules.csv`; refresh it outside Hopper with `scripts/update_geno2pheno_rules.py` and review the diff before use.

Each JSON drug outcome includes `applicable_rule_count`, so an outcome with no applicable mutation rules remains distinguishable from evaluated rules with no matching mutations. The rules snapshot currently contains 16 drugs. Susceptible* presentation remains under investigation; see the reference audit.

## Historical re-call

Preview a lightweight re-call without changing archives:

```bash
python scripts/recall_resistance.py /path/to/sample-archives \
  --rules assets/hcv_geno2pheno_rules.csv \
  --h77-fasta refgenomes/resistance/AF011751.1.fasta
```

Preview reports eligibility and skip reasons. **Old consensus files are not eligible for resistance-only recall**: regenerate consensus and dependent outputs or rerun the pipeline first. Eligible samples require matching converter-4 JSON sidecars for analytical IUPAC and any supplied scaffold, its BLAST output, per-sample QC JSON and one VADR GFF (normalized `_mod.gff` is preferred over raw annotations). Do not synthesize sidecars to approve old consensus sequences. Flat and historical nested `results/` sample archives are supported. Without CRAM, baseline assessment and low-depth qualification are still possible from sidecars, but read evidence is unavailable.

Add `--apply` to stage and validate the entire replacement set before archive changes. All six resistance files, per-sample QC summaries and existing run aggregate JSON/JSONL are backed up first. Metadata unrelated to resistance is preserved. A non-overwriting backup, atomic journal and exclusive archive lock support rollback; caught failures restore old files and remove newly created outputs. File replacement is atomic individually, **not** an atomic filesystem-wide transaction. Sample QC completion markers are published after resistance artifacts; re-import into Virtitta only after completion. Do not run another archive writer or importer concurrently.

After a hard interruption, stop any surviving recall process and inspect `<archives>/.resistance-recall-lock/backup.txt`. Restore explicitly before retrying:

```bash
python scripts/recall_resistance.py /path/to/sample-archives \
  --restore-backup /path/to/sample-archives-resistance-backup-TIMESTAMP
```

Backups contain a `manifest.json` and originals under `files/<archive-relative-path>`. Existing output permissions are preserved. Symlinked replacement targets are rejected. A failure during backup preparation changes no archive files; inspect the incomplete backup before removing an abandoned lock. No real archives are changed by tests.

Both Python entry points default to AF011751.1 when `--h77-fasta` is omitted. The resistance JSON records the reference identifier and file checksum. The workflow and historical re-call select subtype from the first data hit in the 15% IUPAC BLAST output; the re-call preview displays that subtype. The assembly-reference label in the report is not used for rule selection. Region-specific subtype selection and related-subtype interpretation are not yet implemented; see [the reference audit](resistance_reference_audit.md).

## Validation

`python -m unittest discover -s tests -v` covers allele selection, coverage masking, H77 projection, deletions versus missing sequence, frameshift/stop safeguards, rule truth tables, archive rollback/QC refresh, track identifiers, VADR failure propagation, nine original anonymized cases and eleven additional subtype cases submitted to online geno2pheno. The additional cases check all 16 drug rows and resistant/reduced-susceptibility predictions, but do not establish susceptible* or scored-mutation-list parity. Five rare cases are user-identified mixtures, not phased strain-level truth. The private fixture re-identification keys must remain outside Git. Local module tests precede a representative Hopper run through the site wrappers; collection-wide reprocessing and deployment are separate operator actions.
