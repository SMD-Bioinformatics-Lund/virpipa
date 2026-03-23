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
        tuple val(run_name), val(sample_id), path("*.vcf.gz"), path("*.vcf.gz.csi"), emit: vcf
        path "*.vcf.gz.stats", emit: stats
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def bcftools = "apptainer exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools"
        
        """
        # Generate pileup (matching bash pipeline: -Ob with separate call step)
        ${bcftools} mpileup -Ob -f ${ref_fasta} -d 1000000 -a AD,DP -o ${sample_id}-${ref_name}.bcf ${bam}
        ${bcftools} call -Oz -m -A --ploidy 1 -o ${sample_id}-${ref_name}.vcf.gz ${sample_id}-${ref_name}.bcf
        
        # Index
        ${bcftools} index ${sample_id}-${ref_name}.vcf.gz
        
        # Stats
        ${bcftools} stats ${sample_id}-${ref_name}.vcf.gz > ${sample_id}-${ref_name}.vcf.gz.stats
        """
    } else {
        """
        echo "VARIANT_CALLING requires container with bcftools"
        exit 1
        """
    }
}
