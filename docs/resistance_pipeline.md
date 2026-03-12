# HCV Drug Resistance Annotation Pipeline

## Overview

This pipeline annotates HCV VCF variants with drug resistance information from the geno2pheno database. It matches variants against resistance rules based on genotype subtype.

## Files

### 1. `scripts/download_geno2pheno_rules.py`
Downloads the resistance rules table from https://hcv.geno2pheno.org/index.php?page=Rules

Usage:
```bash
python scripts/download_geno2pheno_rules.py -o hbv_result_rules.csv
```

### 2. `scripts/parse_geno2pheno_rules.py`
Parses the rules CSV into a more usable format. Handles:
- Subtype matching: "1" matches all genotype 1 subtypes (1a, 1b, etc.)
- Compound rules: "445F and 451S" means ALL variants must be present

Usage:
```bash
python scripts/parse_geno2pheno_rules.py --subtype 3a
```

### 3. `scripts/annotate_vcf_resistance.py` (Main Script)
Annotates VCF variants with drug resistance information.

**Required arguments:**
- `--vcf` - VCF file (e.g., CMD1065A189-pilon.vcf.gz)
- `--gff` - VADR GFF file (e.g., CMD1065A189.vadr.pass_mod.gff)
- `--fasta` - IUPAC FASTA (e.g., CMD1065A189-0.15-iupac.fasta)
- `--subtype` - HCV subtype (e.g., 3a, 1b)

**Optional arguments:**
- `--sample-name` - Sample ID (default: derived from VCF filename)
- `--output-dir` - Output directory (default: results)
- `--rules` - Rules CSV (default: hbv_result_rules.csv)

**Example:**
```bash
python scripts/annotate_vcf_resistance.py \
    --vcf /path/to/CMD1065A189-pilon.vcf.gz \
    --gff /path/to/CMD1065A189.vadr.pass_mod.gff \
    --fasta /path/to/CMD1065A189-0.15-iupac.fasta \
    --subtype 3a \
    --sample-name CMD1065A189 \
    --output-dir results/
```

**Environment:**
Requires the `skrotis` mamba environment with:
- pysam
- pandas
- biopython

Run with: `mamba run -n skrotis python scripts/annotate_vcf_resistance.py ...`

## Output Files

Output goes to `results/` directory by default.

### 1. `sample_resistance.tsv`
Variant-focused results sorted by genomic position.

Columns:
- sample, gene, genomic_start, genomic_end, ref_nuc, alt_nuc, aa_pos, ref_aa, alt_aa, rule_definition, drugs, prediction, reference, strand

### 2. `sample_resistance.bed`
BED file for IGV visualization. One entry per unique amino acid change. Coordinates cover the full codon.

### 3. `sample_resistance_by_drug.tsv`
Drug-focused results with sections per drug.

## Pipeline Integration

To integrate into the Nextflow pipeline:

1. Add as a process in `modules/local/`
2. Input: tuple of (sample_id, vcf, gff, fasta, subtype)
3. Output: resistance TSV, BED, and by_drug TSV files

## Data Locations (example sample CMD1065A189)

- VCF: `results/vcf/CMD1065A189-pilon.vcf.gz`
- GFF: `results/CMD1065A189.vadr.pass_mod.gff`
- FASTA: `fasta/CMD1065A189-0.15-iupac.fasta`
- Subtype: extracted from blast file or lid file (stored in `$subtype` variable)

## Testing

Tested with sample CMD1065A189 (subtype 3a):
- Results match the PDF from geno2pheno
- BED coordinates correctly span the full codon (e.g., 3900-3903 for NS3:156)

## Notes

- The rules CSV (hbv_result_rules.csv) should be downloaded occasionally to get updated rules
- BED coordinates use codon boundaries (not just variant position)
- Results default to the `results/` folder

## Branch

`hcv-resistance-annotation` - Contains all changes
