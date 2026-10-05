# HCV resistance calling

VirPipa calls resistance from the 15% IUPAC consensus. The VCF is not used for resistance calling: it is relative to the sample's polished consensus and therefore cannot represent majority substitutions reliably.

`scripts/annotate_vcf_resistance.py` aligns translated NS3, NS5A, and NS5B genes to H77 and evaluates the committed geno2pheno rules in H77 amino-acid numbering. IUPAC codons are expanded before translation, simple and compound rules are evaluated explicitly, and missing sequence is reported as not assessed. The default clinical result remains the 15% consensus result.

When the matching CRAM is supplied, the caller also counts complete codons from usable primary alignments. Bases require quality 13, sites require seven complete codons, overlapping mates count once, and conflicting mate codons are discarded. These frequencies are supporting evidence for exploratory 2–50% thresholds; they do not replace the authoritative 15% result.

## Outputs

- `*_resistance.tsv`: resistance-associated substitutions supporting fully matched rules.
- `*_resistance.bed` and `*_resistance.gff`: compatibility tracks containing called resistance substitutions.
- `*_resistance_by_drug.tsv`: every applicable drug with an explicit outcome (`resistant`, `reduced susceptibility`, `susceptible`, `susceptible with scored substitutions`, `not licensed`, or `not assessed`).
- `*_resistance.json`: versioned baseline calls, rule provenance/checksum, compound-rule state, all resistance sites, codon frequencies, counting parameters, and consistency warnings.
- `*_resistance_sites.gff3`: all resistance-rule positions for IGV. Popup attributes include H77 and sample amino-acid positions, codons, observed amino acids, drugs, rules, depth, and frequencies.

The rules snapshot is `assets/hcv_geno2pheno_rules.csv`; refresh it outside Hopper with `scripts/update_geno2pheno_rules.py` and review the diff before use.

## Historical re-call

Preview a lightweight re-call without changing archives:

```bash
python scripts/recall_resistance.py /path/to/sample-archives \
  --rules assets/hcv_geno2pheno_rules.csv \
  --h77-fasta refgenomes/1a-AF009606.fa
```

Add `--apply` to generate outputs in temporary directories and atomically replace resistance files. Existing resistance files are copied to a timestamped backup first. Samples lacking the 15% FASTA, its report, or a single VADR GFF are skipped. If a CRAM is absent, the baseline call is still produced but dynamic frequency evidence is unavailable.

## Validation

`python -m unittest discover -s tests -v` covers IUPAC translation, H77 indel mapping, compound rules, and nine anonymized consensus cases previously submitted to online geno2pheno. The private fixture re-identification key must remain outside Git.
