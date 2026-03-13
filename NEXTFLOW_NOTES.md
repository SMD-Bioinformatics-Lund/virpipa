# virpipa Nextflow Pipeline Notes

## Current Status (2026-03-13)

Pipeline runs end-to-end for single sample with ref_dir.

### Completed Steps
- SUBSAMPLE_READS ✅
- REMOVE_HOSTILE ✅  
- MAP_READS (sentieon UMI workflow) ✅
- SELECT_BEST_REFERENCE ✅
- ASSEMBLE_SPADES ✅
- ASSEMBLE_HYBRID (mummer) ✅
- VARIANT_CALLING (bcftools) ✅
- CREATE_CONSENSUS (bcftools + awk) ✅
- SUBTYPE_BLAST ✅
- ANNOTATE_VADR ✅

### Missing vs Original Bash Pipeline
1. **Pilon polishing** - 10 iteration loop (critical for quality!)
2. **Mapping to hybrid/pilon** - need to map reads to hybrid assembly
3. **Multiple VCF filtering** - different min fractions (0.01, 0.05, 0.1, etc.)
4. **CRAM file creation**
5. **Coverage logging**
6. **ANNOTATE_RESISTANCE** module exists but not connected

### Key Technical Notes

**Channel patterns that work:**
- Use `cross()` + `filter()` instead of `join()` for matching samples
- Use `tuple()` not `[ ]` for Nextflow tuples
- Use `def` for variables in closures, not `params`
- Use `file().toAbsolutePath()` for absolute paths
- Pass ref_dir as part of tuple, not separate val channel

**Container issues:**
- VADR uses `scripts/vadr_annotate.sh` not direct vadr command
- BLAST db at `${ref_dir}/hcvglue/` (directory, not file)
- samtools faidx, not bcftools faidx

### Test Data Locations
- New results: `/mnt/fs1/jonas/src/virpipa/results/test_run/SAMPLE001/`
- Old bash: `/mnt/fs1/jonas/hcv/results/test_run_bash_original/`

### To Run
```bash
cd /fs1/jonas/src/virpipa
nextflow run main.nf -profile hpc --input assets/test_samplesheet.csv --ref_dir refgenomes
```

### For Pilon Polishing
Look at hcvpipe.sh lines 481-510 for the polishing function:
- Map reads to hybrid assembly
- Run pilon 10 times
- Stop when changes plateau
- Create final consensus with iupac codes
