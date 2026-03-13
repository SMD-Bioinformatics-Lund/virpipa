process VARIANT_CALLING {
    tag { "${sample_id}:${ref_name}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/vcf", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(bam), path(bai), path(ref_fasta), val(ref_name)
    
    output:
        tuple val(run_name), val(sample_id), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
        path "*.vcf.gz.stats", emit: stats
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def freebayes = "apptainer exec -B ${bind_paths} ${container_dir}/freebayes_1.3.8.sif freebayes"
        def bcftools = "apptainer exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools"
        
        """
        # Variant calling with freebayes
        ${freebayes} -f ${ref_fasta} --ploidy 1 ${bam} > ${sample_id}-${ref_name}.vcf
        
        # Compress and index
        ${bcftools} view -Oz ${sample_id}-${ref_name}.vcf > ${sample_id}-${ref_name}.vcf.gz
        ${bcftools} index -t ${sample_id}-${ref_name}.vcf.gz
        
        # Stats
        ${bcftools} stats ${sample_id}-${ref_name}.vcf.gz > ${sample_id}-${ref_name}.vcf.gz.stats
        """
    } else {
        """
        echo "VARIANT_CALLING requires container with freebayes and bcftools"
        exit 1
        """
    }
}
