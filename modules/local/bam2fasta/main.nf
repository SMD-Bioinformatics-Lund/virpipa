process BAM2FASTA {
    tag { "${sample_id}:${ref_name}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/pilon", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(bam), path(bai), path(ref_fasta), val(ref_name)
        val ambiguity_threshold
    
    output:
        tuple val(run_name), val(sample_id), path("${sample_id}-${ref_name}.fasta"), path("${sample_id}-${ref_name}.fasta.fai"), emit: fasta
        path "${sample_id}-${ref_name}.vcf.gz", emit: vcf
        path "${sample_id}-${ref_name}.vcf.gz.csi", emit: vcf_index
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    def bcftools = container_dir ? 
        "apptainer exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools" :
        "bcftools"
    
    def samtools = container_dir ? 
        "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"
    
    """
    # Generate pileup and call variants
    ${bcftools} mpileup -Ob -f ${ref_fasta} -d 1000000 -a AD,DP ${bam} > ${sample_id}.bcf
    ${bcftools} call -Oz -m -A --ploidy 1 -o ${sample_id}-${ref_name}.vcf.gz ${sample_id}.bcf
    ${bcftools} index ${sample_id}-${ref_name}.vcf.gz
    
    # Decompress VCF for awk
    ${bcftools} view -O v ${sample_id}-${ref_name}.vcf.gz > input.vcf
    
    # Create consensus with IUPAC codes
    awk -v MIN_AF=${ambiguity_threshold} -v MIN_DP=7 -f ${projectDir}/scripts/vcf_to_iupac.awk input.vcf ${ref_fasta} > ${sample_id}-${ref_name}.fasta
    
    # Fix header
    sed -i 's/>.*/>${sample_id}/' ${sample_id}-${ref_name}.fasta
    
    # Index
    ${samtools} faidx ${sample_id}-${ref_name}.fasta
    """
}
