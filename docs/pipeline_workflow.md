# HCV Pipeline Workflow Documentation

This document describes the complete bash pipeline workflow for reference when implementing the Nextflow version.

## Test Data
- Sample: SAMPLE001
- Input: `/fs1/jonas/hcv/test_data/SAMPLE001_122-634521_S26_R1_001.fastq.gz`
- Reference genomes: 12 HCV subtypes in `/fs1/jonas/hcv/refgenomes/`
- Container dir: `/fs1/resources/containers/`

## Bash Pipeline Complete Workflow (hcvpipe.sh)

### Step 0: Setup & Initialization
- Parse command line arguments
- Set up container tools (sentieon, samtools, bcftools, hostile, spades, mummer, blast, python, pilon, etc.)
- Create output directories (fastq, bam, bcf, vcf, fasta, spades, mummer, tmp, results)

### Step 1: Remove Human Reads (if -H not set)
- Command: `hostile clean --offline --fastq1 {r1} --fastq2 {r2} --out-dir {outdir}/fastq`
- Output: `${outdir}/fastq/*_clean_*.fastq.gz`
- Also creates: `${outdir}/results/hostile.json`

### Step 1b: Optional Subsampling (if -s specified, default 500000)
- Command: `seqtk sample -s100 {read} {nofreads} | pigz > {output}`
- Output: `${outdir}/fastq/*_sub.fastq.gz`

### Step 2: Subsample for Reference Selection (250k reads, always)
- Command: `seqtk sample -s100 {read} 250000 | pigz > {output}`
- Output: `${outdir}/tmp/*_sub.fastq.gz`

### Step 3: Map to All Reference Genomes & Create FASTAs
For each of 12 reference genomes (1a-AF009606, 1a-M62321, 1b-D90208, 2a-D00944, 2b-D10988, 3a-D17763, 3k-HPCJK049E1, 4a-GU814265, 5a-Y13184, 6a-Y12083, 6g-HPCJK046E2, 7a-EF108306):

#### 3a: umimap() - Map with UMI consensus
```
sentieon umi extract -d 3M2S+T,3M2S+T {r1} {r2} | \
sentieon bwa mem \
    -R "@RG\\tID:{id}\\tSM:{id}\\tLB:{id}\\tPL:illumina" \
    -t {cpus} -k 11 -B 2 -L 25 \
    -p -C {ref} - | \
sentieon umi consensus --copy_tags XR,RX,MI,XZ -o {consensus_fastq}

sentieon bwa mem \
    -R "@RG\\tID:{id}\\tSM:{id}\\tLB:{id}\\tPL:illumina" \
    -t {cpus} -k 11 -B 2 -L 25 \
    -p -C {ref} {consensus_fastq} | \
sentieon util sort -i - --sam2bam --umi_post_process -o {bamfile}

samtools view -@ 30 {bamfile} -e "sclen < 30" --with-header -b -o {filterbam}
samtools index {filterbam}
samtools stats {filterbam} > {filterbam}.stats
```
- Output: `${outdir}/bam/${id}-${genome}.r11b2L25.bwa.umi.filter.sort.bam`

#### 3b: bam2fasta() - Create consensus FASTA
```
bcftools mpileup -Ob -f {ref} -d 1000000 -a AD,DP -o {outdir}/bcf/{id}-{genome}.bcf {bam}
bcftools call -Oz -m -A --ploidy 1 -o {outdir}/vcf/{id}-{genome}.vcf.gz {outdir}/bcf/{id}-{genome}.bcf
bcftools index {outdir}/vcf/{id}-{genome}.vcf.gz
zcat {outdir}/vcf/{id}-{genome}.vcf.gz > {outdir}/vcf/{id}-{genome}.vcf
awk -v MIN_AF={ambiguity} -v MIN_DP=7 -f {scripts}/vcf_to_iupac.awk {vcf} {ref} > {outdir}/fasta/{id}-{genome}.fasta
sed -i "s/>.*/>{id}-{genome}/" {outdir}/fasta/{id}-{genome}.fasta
bcftools stats {outdir}/vcf/{id}-{genome}.vcf.gz > {outdir}/vcf/{id}-{genome}.vcf.gz.stats
```
- Output: `${outdir}/fasta/${id}-${genome}.fasta`

### Step 4: Select Best Reference Subtype
- Command: Parse BAM stats files, find lowest error rate
- Output: `subtype` variable (e.g., "3a-D17763")

### Step 5: Map to Best Reference with All Reads
- Copy best reference FASTA to results/
- Create CRAM from best reference BAM

### Step 6: De Novo Assembly (SPAdes)
- Command: `spades.py -1 {r1} -2 {r2} --rnaviral -o {outdir}/spades --threads {cpus}`
- Output: `${outdir}/spades/contigs.fasta` -> renamed to `${outdir}/spades/${id}.spades`

### Step 7: Build Hybrid Reference (MUMmer + custom script)
```
mkdir -p {outdir}/mummer
nucmer --maxmatch -p {id} {ref}/{subtype}.fa {outdir}/spades/{id}.spades
delta-filter -q {outdir}/mummer/{id}.delta > {id}.delta-filter
show-tiling {outdir}/mummer/{id}.delta-filter > {id}.tiling
python {build_hybrid_reference.py} {ref}/{subtype}.fa {outdir}/spades/{id}.spades {outdir}/mummer/{id}.tiling {id}
samtools faidx {outdir}/mummer/{id}.hybrid.fasta
sentieon bwa index {outdir}/mummer/{id}.hybrid.fasta
```
- Output: `${outdir}/mummer/${id}.hybrid.fasta`

### Step 8: Pilon Polishing (Loop, max 10 iterations)
For each iteration:
```
# Map to current pilon FASTA
umimap() or umimapnoopt() - similar to Step 3a but with pilon FASTA

# Run pilon
pilon --fix all --mindepth 5 --changes --genome {ref} --frags {paired.bam} --unpaired {unpaired.bam} --outdir {outdir}/pilon --output {id}-pilon-{round}
sed -i -e 's/_pilon//' {outdir}/pilon/{id}-pilon-{round}.fasta
samtools faidx {outdir}/pilon/{id}-pilon-{round}.fasta
sentieon bwa index {outdir}/pilon/{id}-pilon-{round}.fasta

# Check for convergence - stop if no changes
```
- Early stop if `${id}-pilon-{i}.changes` == `${id}-pilon-{i-1}.changes`
- Final outputs:
  - `${outdir}/pilon/${id}.fasta` (final polished FASTA)
  - `${outdir}/pilon/${id}-iupac.fasta`

### Step 9: Map to Final Pilon FASTA (umimapnoopt)
```
sentieon bwa index {pilon_fasta}  # if not exists
sentieon umi extract -d 3M2S+T,3M2S+T {r1} {r2} | \
sentieon bwa mem -R "@RG\\tID:{id}\\tSM:{id}\\tLB:{id}\\tPL:illumina" -t {cpus} -p -C {ref} - | \
sentieon umi consensus --copy_tags XR,RX,MI,XZ -o {consensus_fastq}
sentieon bwa mem -R "@RG\\tID:{id}\\tSM:{id}\\tLB:{id}\\tPL:illumina" -t {cpus} -p -C {ref} {consensus_fastq} | \
sentieon util sort -i - --sam2bam --umi_post_process -o {bamfile}
samtools view -@ 30 {bamfile} -e "sclen < 30" --with-header -b -o {filterbam}
samtools index {filterbam}
samtools stats {filterbam} > {filterbam}.stats
```
- Output: `${outdir}/bam/${id}-pilon.r11b2L25.bwa.umi.filter.sort.bam`

### Step 10: Create bam2fasta Consensus (Majority Call Correction)
- Same as Step 3b but with pilon BAM and pilon FASTA
- Parameters: MIN_AF=1.0
- Then: mv original pilon FASTA to .org, cp bam2fasta output to pilon FASTA, fix header
- Output: `${outdir}/pilon/${id}.fasta` (replaced with 1.0-iupac)

### Step 11: Variant Calling (on polished genome)
```
bcftools mpileup -Ob -f {pilon_fasta} -d 1000000 -a AD,DP -o {outdir}/bcf/{id}-pilon.bcf {bam}
bcftools call -m -A --ploidy 1 -Oz -o {outdir}/vcf/{id}-pilon.vcf.gz {outdir}/bcf/{id}-pilon.bcf
bcftools index {outdir}/vcf/{id}-pilon.vcf.gz
```
- Output: `${outdir}/vcf/${id}-pilon.vcf.gz`

### Step 12: Filter VCF at Multiple Thresholds
For minfrac in 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4:
```
bcftools filter -i "FMT/AD[0:1] / FMT/DP[0] >= {minfrac} | FMT/AD[0:2] / FMT/DP[0] >= {minfrac} | FMT/AD[0:3] / FMT/DP[0] >= {minfrac}" \
    {vcf} -Oz -o {outdir}/vcf/{id}-pilon-m{minfrac}.vcf.gz
bcftools index {outdir}/vcf/{id}-pilon-m{minfrac}.vcf.gz
bcftools stats {outdir}/vcf/{id}-pilon-m{minfrac}.vcf.gz > {outdir}/vcf/{id}-pilon-m{minfrac}.vcf.gz.stats
```
- Output: `${outdir}/vcf/${id}-pilon-m{minfrac}.vcf.gz`

### Step 13: Create 15% IUPAC Consensus
- Same as Step 3b but with MIN_AF=0.15
- Output: `${outdir}/fasta/${id}-0.15-iupac.fasta`

### Step 14: Create CRAM Files
```
samtools view -O cram,embed_ref -T {ref_fasta} {bam} -o {cram}
samtools index {cram}
```
CRAMs created:
- `${outdir}/results/${id}.cram` (pilon FASTA)
- `${outdir}/results/${id}-0.15-iupac.cram` (0.15-iupac FASTA)
- `${outdir}/results/${id}-${subtype}.cram` (best reference)

### Step 15: Copy Results to results/ Directory
Files copied to `${outdir}/results/`:
- `${id}.fasta` (pilon FASTA)
- `${id}.fasta.fai`
- `${id}-0.15-iupac.fasta`
- `${id}-pilon-iupac.fasta`
- `${id}-pilon-m*.vcf.gz*` (all filtered VCFs)
- `${id}.cram`, `${id}-0.15-iupac.cram`, `${id}-${subtype}.cram`

### Step 16: Create Reports (createreport)
- Extract VCF stats: SNP count, multiallelic, AF0, AF99
- samtools coverage -> datamash transpose
- Nucleotide frequency from FASTA
- Output: `${outdir}/results/${id}.report.tsv`

### Step 17: BLAST Subtyping
```
blastn -query {fasta} -db {refdir}/hcvglue/hcvglue -outfmt 6 > {fasta}.blast
```
- Output: `${outdir}/results/*.blast`

### Step 18: VADR Annotation
- Command: `vadr_annotate.sh {fasta} {outdir} {id}`
- Output: `${outdir}/results/${id}*_mod.gff`, `${outdir}/results/${id}*bed`

### Step 19: Coverage Statistics
- Command: `samtools coverage` with KDE/RUG plots

## Implementation Notes for Nextflow

1. `container` directive: Use Nextflow's `container` directive for processes when practical.
2. Bash scripting: Use bash in process scripts instead of Groovy-heavy command assembly.
3. File naming: Keep the same naming conventions as bash for parity checks.
4. File comparison: Do not rely solely on checksums because some outputs embed working paths.
5. Output directories: Mimic the bash directory structure where it matters for comparison.
