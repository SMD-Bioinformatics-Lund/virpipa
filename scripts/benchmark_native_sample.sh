#!/bin/bash
# Experimental native-tool reconstruction and fixed-consensus-read remapping.
# Run only in a new task directory through benchmark_native.nf.
set -euo pipefail
[[ $# -eq 8 || $# -eq 9 ]] || { echo 'Expected REF CRAM IUPAC MAJORITY SAMPLE CONTAINERS SCRIPTS CPUS [RESISTANCE_JSON]' >&2; exit 2; }
source_ref=$(readlink -f "$1")
source_cram=$(readlink -f "$2")
source_iupac=$(readlink -f "$3")
source_majority=$(readlink -f "$4")
sample=$5
tools_dir=$6
scripts_dir=$7
threads=$8
runtime=${VIRPIPA_BENCHMARK_RUNTIME:-apptainer}
binds=${VIRPIPA_BENCHMARK_BINDS:-/fs1,/fs2,/local,/mnt/beegfs}
samtools() { "$runtime" exec -B "$binds" "$tools_dir/samtools_1.21.sif" samtools "$@"; }
bcftools() { "$runtime" exec -B "$binds" "$tools_dir/bcftools_1.21.sif" bcftools "$@"; }
ivar() { "$runtime" exec -B "$binds" "$tools_dir/ivar_1.4.4.sif" ivar "$@"; }
bwa() { "$runtime" exec -B "$binds" "$tools_dir/sentieon_202308.03.sif" sentieon bwa "$@"; }
python_tool() { "$runtime" exec -B "$binds" "$tools_dir/python_hcvpipe.sif" python "$@"; }

mkdir native
cp -L "$source_ref" native/scaffold_control.fasta
cp -L "$source_iupac" native/baseline_iupac.fasta
cp -L "$source_majority" native/baseline_majority.fasta
sha256sum "$source_ref" "$source_cram" "$source_iupac" "$source_majority" > native/inputs.sha256
sha256sum "$tools_dir"/{bcftools_1.21,samtools_1.21,ivar_1.4.4,sentieon_202308.03,python_hcvpipe}.sif > native/containers.sha256
{ bcftools --version; samtools --version; ivar version; } > native/versions.txt
samtools faidx native/scaffold_control.fasta
# Freeze exactly the eligibility used by the read-codon counter: mapped primary,
# nonduplicate, QC-passing reads. Both native callers and all remappings use it.
samtools view -b -F 0xF04 -T native/scaffold_control.fasta -o native/source.primary.bam "$source_cram"
samtools index native/source.primary.bam
bcftools mpileup -A -B -Q 13 -d 1000000 -a FORMAT/AD,FORMAT/DP \
    -f native/scaffold_control.fasta -Ob -o native/calls.bcf native/source.primary.bam
bcftools call -m -A --ploidy 1 -Oz -o native/calls.vcf.gz native/calls.bcf
bcftools view native/calls.vcf.gz > native/calls.vcf
bcftools norm -f native/scaffold_control.fasta -Oz -o native/normalized.vcf.gz native/calls.vcf.gz
bcftools index native/normalized.vcf.gz
called_sample=$(bcftools query -l native/normalized.vcf.gz)
bcftools consensus -f native/scaffold_control.fasta -s "$called_sample" -H I \
    --chain native/bcftools_genotype.chain native/normalized.vcf.gz > native/bcftools_genotype.fasta

samtools mpileup -aa -A -B -Q 0 -d 1000000 -f native/scaffold_control.fasta native/source.primary.bam > native/reads.pileup
for entry in '0:ivar_majority' '0.85:ivar_085' '0.90:ivar_090'; do
    cutoff=${entry%%:*}
    name=${entry#*:}
    ivar consensus -q 13 -m 1 -n N -c 0.8 -t "$cutoff" -i "${sample}_${name}" \
        -p "native/$name" < native/reads.pileup
    cp "native/$name.fa" "native/$name.fasta"
done
python_tool "$scripts_dir/consensus_metadata.py" --vcf native/calls.vcf \
    --reference native/scaffold_control.fasta --output native/awk_baqoff_iupac.fasta \
    --identifier "${sample}_awk_baqoff_iupac" --mode iupac --min-af 0.15
python_tool "$scripts_dir/consensus_metadata.py" --vcf native/calls.vcf \
    --reference native/scaffold_control.fasta --output native/awk_baqoff_majority.fasta \
    --identifier "${sample}_awk_baqoff_majority" --mode analytical_majority

# These are already UMI consensus reads. Do not run UMI consensus a second time.
samtools collate -u -O native/source.primary.bam |
    samtools fastq -n -1 native/reads.R1.fastq -2 native/reads.R2.fastq \
        -0 native/reads.other.fastq -s native/reads.single.fastq -
cat native/reads.other.fastq native/reads.single.fastq > native/reads.unpaired.fastq
sha256sum native/reads.*fastq > native/fixed_reads.sha256
rg="@RG\\tID:$sample\\tSM:$sample\\tPL:illumina"
for name in scaffold_control baseline_iupac baseline_majority bcftools_genotype \
            ivar_majority ivar_085 ivar_090 awk_baqoff_iupac awk_baqoff_majority; do
    ref="native/$name.fasta"
    samtools faidx "$ref"
    bwa index "$ref"
    chunks=()
    if [[ -s native/reads.R1.fastq && -s native/reads.R2.fastq ]]; then
        bwa mem -t "$threads" -R "$rg" "$ref" native/reads.R1.fastq native/reads.R2.fastq |
            samtools sort -o "native/$name.paired.bam"
        chunks+=("native/$name.paired.bam")
    fi
    if [[ -s native/reads.unpaired.fastq ]]; then
        bwa mem -t "$threads" -R "$rg" "$ref" native/reads.unpaired.fastq |
            samtools sort -o "native/$name.unpaired.bam"
        chunks+=("native/$name.unpaired.bam")
    fi
    if [[ ${#chunks[@]} -eq 0 ]]; then
        samtools dict "$ref" | samtools view -b -o "native/$name.merged.bam" -
    else
        samtools merge -f "native/$name.merged.bam" "${chunks[@]}"
    fi
    samtools view -e 'sclen < 30' -b -o "native/$name.bam" "native/$name.merged.bam"
    samtools index "native/$name.bam"
    samtools view -O cram,embed_ref -T "$ref" -o "native/$name.cram" "native/$name.bam"
    samtools index "native/$name.cram"
    samtools depth -aa "native/$name.bam" > "native/$name.depth.tsv"
done
if [[ $# -eq 9 ]]; then
    python_tool "$scripts_dir/benchmark_native_score.py" --native-dir native \
        --archived-resistance "$9" --output native/scores.json
fi
sha256sum -c native/inputs.sha256 > native/source_verification.log
touch native/COMPLETE
