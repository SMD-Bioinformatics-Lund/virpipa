process CREATE_CONSENSUS {
    tag { "${sample_id}:${vcf}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf), path(fasta), path(fai)
        val min_freq
    
    output:
        path "*iupac.fasta", emit: consensus
        path "*.fasta", emit: fasta
        path "*.fai", emit: fai
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def scripts_dir = params.scripts_dir ?: '${projectDir}/scripts'
    
    def bcftools = container_dir ? 
        "apptainer exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools" :
        "bcftools"
    
    def samtools = container_dir ? 
        "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"
    
    """
    # Decompress VCF to uncompressed format
    ${bcftools} view -O v ${vcf} > ${vcf.baseName}
    
    # Use AWK script to create IUPAC consensus
    awk -v MIN_AF=${min_freq} -v MIN_DP=7 -f ${scripts_dir}/vcf_to_iupac.awk ${vcf.baseName} ${fasta} > ${sample_id}_consensus_iupac.fasta
    
    # Index the fasta using samtools
    ${samtools} faidx ${sample_id}_consensus_iupac.fasta
    """
}
