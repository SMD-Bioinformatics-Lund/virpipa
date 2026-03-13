# virpipa Nextflow Pipeline Notes

## Current Status (2026-03-13)

Pipeline runs end-to-end with verified identical output to bash pipeline.

### Completed Steps
- SUBSAMPLE_READS ✅
- REMOVE_HOSTILE ✅  
- MAP_READS (sentieon UMI workflow) ✅
- SELECT_BEST_REFERENCE ✅
- ASSEMBLE_SPADES ✅
- ASSEMBLE_HYBRID (mummer) ✅
- POLISH_PILON_LOOP (10 iterations, convergence check) ✅
- VARIANT_CALLING (bcftools) ✅
- CREATE_CONSENSUS (bcftools + awk) ✅
- SUBTYPE_BLAST ✅
- ANNOTATE_VADR ✅

### Verified Output (vs bash pipeline)
- Polished assembly: 9456 bp (identical)
- Pilon rounds: converges correctly
- BAM, VCF, VADR, BLAST outputs: all present ✅

### Missing vs Original Bash Pipeline
1. ~~**Pilon polishing**~~ - Now implemented!
2. **Multiple VCF filtering** - different min fractions (0.01, 0.05, 0.1, etc.)
3. **CRAM file creation**
4. **Coverage logging**
5. **ANNOTATE_RESISTANCE** module exists but not connected

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
