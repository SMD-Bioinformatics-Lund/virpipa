# Resistance validation handoff

## IUPAC IGV review gate

After publishing the reviewed commits to Hopper and updating the working checkout,
run the `finalize_results` fixture through the site test wrapper in routine and
debug modes, using separate output/work directories outside the source tree.
Validate each archive with `scripts/check_publish_layout.py` and re-import one
fresh real-sample run into the test Virtitta instance before processing more strains.

In both desktop IGV and webIGV verify:

- Genome and reads are the `-0.15-iupac` pair, with no missing-contig/index errors.
- All-site, resistance and gene annotations load; mutation links land on the right codon.
- Resistance popups retain amino acids, codons, depth and compound-rule details.
- All six displayed VCF thresholds remain available and labelled scaffold-derived.
  REF can differ from the IUPAC genome; these are not newly called variants.
- An older archive remains an explicitly labelled scaffold view.

Then run additional strains, including positive resistance, mixtures, low coverage
and coding-integrity findings. Compare drug outcomes to the online reports using
the exact submitted IUPAC FASTA. Keep any fixtures anonymous and keys outside git.

This is the next-stage checklist for the `0.3.0` analytical-consensus contract.
Passing regression tests is not biological or clinical validation. BLAST improvements
are deferred in [their own plan](plans/blast_subtyping.md).

## Local gates

### Fresh Hopper/import validation and reproducibility investigation, 2026-11-06

- The fresh 84e4283 Hopper run completed under Nextflow 26.04.3: 236 tasks completed,
  all exit 0. All six archives passed published-layout, ranked-typing, database/input
  provenance and CRAM/reference checks. Ranked typing is identical across typing JSON,
  QC and resistance evidence. All subtypes and drug outcomes match the prior batch.
- The user imported all six into Virtitta cbf7f0c and confirmed typing, expandable
  alternatives, preservation of review state/categories/comments and successful IGV
  opening. This closes the fresh-run/import check, not the reproducibility question.
- Retained scaffolds differ between runs. The subsequent investigation compared three
  cases and completed two fixed-input Pilon replay jobs and one repeated Hostile job.
  Pilon was repeatable for those fixed inputs. Hostile retained identical read pairs
  in different orders, and seeded subsampling selected different read subsets. The
  control below the sampling cutoff retained all pairs and its final sequence agreed.
- Read-support review also found low-depth terminal differences and internal high-depth
  differences outside the resistance targets. Case A changed VADR fail/pass state;
  case B shifted target-gene coordinates by two bases. Resistance-site codons, amino
  acids and assessment states agreed in all three investigated cases. Full private
  evidence remains outside Git. See [the investigation](scaffold_reproducibility.md).
- At the investigation stage no pipeline sampling/polishing behaviour had been changed;
  the subsequent ordering correction and acceptance decision are recorded below.
- Subsequent native Hostile `--reorder` experiments passed for both differing cases:
  identical ordered filtered pairs across repetitions, unchanged pair membership versus
  baseline, and identical seeded 250,000/1,000,000-pair subsets. The diagnostic wrapper
  supports the option, and the pipeline now enables it in both execution paths.
  The tiny Hostile module fixture completed with 100 input pairs and eight retained
  pairs; mate IDs, output filenames and the existing QC parser were verified. Strict
  Nextflow lint passed on all 25 files using locally installed 25.10.4. The subsequent
  full-sample gate passed: pinned aaa2633 was run twice on case A under 26.04.3 in
  independent fresh directories, with 41 completed tasks per run and no failed/cached
  analysis tasks. All 38 semantic comparisons agreed: ordered filtering, both sampled
  inputs, assemblies, every Pilon round, canonical BAM/CRAM records, final FASTAs and
  sidecars, variant payloads, VADR, typing/resistance JSON and semantic QC. Both followed
  21 → 3 → 2 → 0 → 0 changes. Layout and actual CRAM/reference MD5 checks passed;
  source and raw inputs remained unchanged. All 16 drug outcomes matched the earlier
  uncorrected case. Initial coordinator startup failures were recovered before analysis
  by explicit module initialization in private script copies. This is a one-case
  reproducibility result, not a universal claim. The user accepted retaining the current
  Pilon stopping behaviour for this release; cycle detection/reporting and early cycle
  termination are deferred because breaking earlier may select a different sequence.
  Deacon replacement is
  tracked separately in [the filtering plan](plans/human_read_filtering.md).

### BLAST integration checkpoint, 2026-10-30

- Python regression suite: 87 tests run, 6 skipped (Nextflow was not on the Python
  test runner's PATH), no failures. New contract tests cover unknown typing, QC
  evidence, stale re-call assignments and analytical sequence binding.
- Strict entry-point lint passed separately on locally installed Nextflow 25.10.4:
  25 files without errors. Production 26.04.3 verification remains a Hopper gate.
- Container-backed `subtype` and connected `resistance_typing` fixtures completed.
  The connected fixture selected 3a / D17763, emitted 16 drug outcomes and carried
  identical typing evidence into resistance, including all six database checksums.
  This fixture uses the existing explicitly unverified legacy comparison consensus;
  it does not establish biological concordance or production readiness.
- Virtitta integration subsequently passed all 187 tests, including seven new contract
  tests for assigned/uncertain/genotype-only typing, legacy imports, alternative display,
  invalid declarations and review/override preservation. Its evidence reader also
  accepted the actual connected-pipeline typing JSON (3a / D17763, six database hashes).
- Fresh Hopper/import validation remains pending. Local fixture results do not establish
  biological concordance. Keep old and new archives separate for release comparisons.

### Local verification recorded 2026-10-20

- VirPipa: 60 regression tests passed; Virtitta: 163 passed, including actual tree
  TSV/GrapeTree JSON outcome checks and pipeline-version import compatibility.
- Container-backed `consensus_contract`, `resistance`, `qc_summary` and `finalize_results`
  checks passed; the final published layout passed `check_publish_layout.py`.
- Both connected runs had three valid consensus sidecars with a shared coordinate
  origin and six resistance artifacts. The standalone resistance fixture had all
  16 drug outcomes and codon counting parameters. Sample and aggregate QC JSON
  recorded pipeline version `0.3.0`.
- NS5A retained all record identifiers, ungapped sequences and H77 coordinate maps;
  NS3 and NS5B were byte-identical. Python compilation and diff hygiene passed.
- These are fixture/regression checks, not a fresh full analysis from reads. The final
  publication fixture includes legacy artifacts; the connected consensus fixture
  separately checks the new contract. No Hopper run, new online comparison,
  production deployment or collection-wide reprocessing has been performed.

### Commands

Use the `skrotis` environment and repo-local containers for VirPipa. Run Virtitta's
suite in its own environment. Follow [the module reference](module_test_reference.md).

```bash
nextflow lint main.nf test_module.nf
python -m unittest discover -s tests -v
NXF_OFFLINE=true nextflow run test_module.nf -profile local_containers --module consensus_contract --publish_mode debug --outdir test_output_resistance_completion_contract
NXF_OFFLINE=true nextflow run test_module.nf -profile local_containers --module resistance --publish_mode debug --outdir test_output_resistance_completion_calls
NXF_OFFLINE=true nextflow run test_module.nf -profile local_containers --module qc_summary --virtitta_import_dir test_output_resistance_completion_qc/import_markers --outdir test_output_resistance_completion_qc
NXF_OFFLINE=true nextflow run test_module.nf -profile local_containers --module finalize_results --virtitta_import_dir test_output_resistance_completion_publish/import_markers --outdir test_output_resistance_completion_publish
python scripts/check_publish_layout.py --sample-dir test_output_resistance_completion_publish/fixture_run/SAMPLE001 --sample-id SAMPLE001 --mode routine
```

Run lint without `-format`, using the production Nextflow version before a Hopper
full-run submission. Module fixtures do not compile `main.nf`; they cannot replace
this gate. Completion-handler regression tests require Nextflow on PATH and run
only an isolated stub, including one intentional local task failure.

The QC/final-publication checks route import markers into their local test output
directories; these fixtures must not publish into `/fs1/results/cron/virtitta`.
Do not pass an empty CLI argument for this check: Nextflow 25.10.4 treats the bare
option as boolean true rather than disabling it. The marker's normal deployment
path is deliberately unchanged by this resistance update.

Check actual published outputs, not only work directories or successful task counts:

- Connected fixture: the same sample ID in two different runs retains separate
  sidecars and outputs; all three FASTA roles retain their intended semantics.
- Resistance: six artifacts, all 16 drug rows, coordinates, warnings and codon evidence.
  The sparse standard fixture is not a fully covered negative control.
- QC/final publication: `pipeline_version` is `0.3.0`, sidecars and relative output
  paths resolve, and run aggregates retain sample versions. Standalone QC generation
  without an explicit version records null.
- NS5A: all-gap columns are absent; reference sequences, every H77 coordinate map,
  calls and assessment states are unchanged. Alignment-file hashes necessarily differ.
- Virtitta: incomplete assessment is not displayed as absence of resistance in tree
  metadata or GrapeTree; scored substitutions alone are not positive resistance calls.

## Representative Hopper run: operator actions

1. Review and commit the tested changes, then publish the exact tested ref to the
   Hopper bare repository. Follow [deployment guidance](deployment.md); no deployment
   or remote writes are performed by local tests.
2. Use the existing site wrappers, not ad hoc Nextflow commands. `runme_test.sh`
   pulls updates, but the currently checked-in `runme.sh` does not: explicitly verify
   that its working checkout has the tested commit before launching it.
3. For the real-data run, prepare a small approved samplesheet and a new output
   location in an operator-controlled copy of the wrapper. The checked-in `runme.sh`
   hard-codes its samplesheet and uses the default output directory; do not assume it
   accepts input/output arguments. Retain its site module/environment setup.
4. Start fresh without `resume`. Regenerate consensus and dependent outputs from reads;
   old analytical FASTAs are not approved for resistance-only recall. Do not manufacture
   sidecars to make them eligible. Leave existing archives untouched.
5. Cover well-supported positive/negative examples, rare subtypes, low coverage and
   mixed populations. Use the nine supplied online cases and available rare cases
   for comparison. Original identifiers and re-identification keys stay outside Git.

## Review and acceptance

Compare each sample's 15% analytical sequence with the actual sequence submitted to
geno2pheno. If sequences differ, distinguish input differences from caller differences.
Record the online report date, available service version, per-region subtype/reference,
drug outcomes, supporting substitutions and positions. Do not imply full susceptible*
or scored-substitution parity from the existing resistant-call regression tests.

Inspect low-depth qualifications, missing sites, frameshifts/stops and indel evidence.
Mixed-population samples are challenge cases, not phased strain-level truth. Check
that uncertainty remains uncertainty rather than becoming a resistant or susceptible
call. Inspect the combined desktop/webIGV tracks against the scaffold and confirm
codon/amino-acid tooltips and genomic coordinates.

Import completed outputs into a disposable Virtitta database first. Verify sample
details, tree/GrapeTree metadata, FASTA export roles and review-preserving re-import.
Store an anonymized comparison table with discrepancies classified as input, typing,
coordinate mapping, rule interpretation or presentation; keep private mappings outside
Git. No clinical-accuracy claim follows from a small concordance check.

Proceed to collection-wide reprocessing only after all local gates pass and every
representative discrepancy is explained or explicitly accepted as a documented
limitation. An unexplained wrong call, shifted coordinate or false clear summary
blocks expansion; retain the original archives and investigate first.
