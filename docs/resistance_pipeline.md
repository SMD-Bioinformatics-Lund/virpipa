# HCV Drug Resistance Annotation Pipeline

## Overview

This pipeline annotates HCV VCF variants with drug resistance information from the geno2pheno database. It matches variants against resistance rules based on genotype subtype.

## Files

### 1. `scripts/update_geno2pheno_rules.py`
Downloads the resistance rules table from https://hcv.geno2pheno.org/index.php?page=Rules and can also rebuild normalized data from an existing CSV snapshot.

Usage:
```bash
python scripts/update_geno2pheno_rules.py --output-csv assets/hcv_geno2pheno_rules.csv
```

### 2. `scripts/annotate_vcf_resistance.py` (Main Script)
Annotates VCF variants with drug resistance information.

**Required arguments:**
- `--vcf` - VCF file (e.g., SAMPLE001-pilon.vcf.gz)
- `--gff` - VADR GFF file (e.g., SAMPLE001.vadr.pass_mod.gff)
- `--fasta` - IUPAC FASTA (e.g., SAMPLE001-0.15-iupac.fasta)
- `--subtype` - HCV subtype (e.g., 3a, 1b)

**Optional arguments:**
- `--sample-name` - Sample ID (default: derived from VCF filename)
- `--output-dir` - Output directory (default: results subfolder of sample)
- `--rules` - Rules JSON or CSV (default: assets/hcv_geno2pheno_rules.csv)
- `--assets-dir` - Directory for reference files (default: assets)
- `--ref-bed` - Also generate reference BED with all resistance positions

**Example:**
```bash
python scripts/annotate_vcf_resistance.py \
    --vcf /path/to/vcf/SAMPLE001-pilon.vcf.gz \
    --gff /path/to/results/SAMPLE001.vadr.pass_mod.gff \
    --fasta /path/to/fasta/SAMPLE001-0.15-iupac.fasta \
    --subtype 3a \
    --sample-name SAMPLE001 \
    --rules assets/hcv_geno2pheno_rules.csv
```

**Output:**
- Results go to sample's `results/` folder (e.g., `results/SAMPLE001_resistance.tsv`)
- Reference files go to `assets/` folder (e.g., `assets/hcv_geno2pheno_rules.csv`, `assets/resistance_reference.bed`)

**Environment:**
Use the `skrotis` environment or the pipeline python container.

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

The resistance module is now wired in the Nextflow pipeline. It consumes:
- `SAMPLE001-pilon-m0.15.vcf.gz` for the resistance-call set
- `SAMPLE001.vadr.pass_mod.gff` for gene coordinates
- `SAMPLE001-0.15-iupac.fasta` for codon translation
- subtype parsed from `SAMPLE001-0.15-iupac.fasta.blast`

## Data Locations (example sample SAMPLE001)

- VCF: `results/SAMPLE001-pilon-m0.15.vcf.gz`
- GFF: `results/SAMPLE001.vadr.pass_mod.gff`
- FASTA: `results/SAMPLE001-0.15-iupac.fasta`
- Subtype: extracted from the first hit in `results/SAMPLE001-0.15-iupac.fasta.blast`

## Testing

Fixture-backed positive test data lives in `assets/test_data/resistance/` and uses a synthetic `NS5A 93H` variant for subtype `3a`.
Run with:

```bash
nextflow run test_module.nf -profile local_containers --module resistance
```

## Notes

- The geno2pheno rules table should be refreshed manually outside the HPC when needed
- BED coordinates use codon boundaries (not just variant position)
- Results default to the `results/` folder

## Notes

- The pipeline currently consumes the committed rules CSV by default through `params.resistance_rules`
- A normalized JSON artifact is still supported if you want a stricter machine-facing rules contract later
