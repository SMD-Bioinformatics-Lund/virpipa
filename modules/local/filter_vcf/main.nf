process FILTER_VCF {
    tag { "${sample_id}" }
    label 'process_medium'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/vcf", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf), path(tbi), path(ref_fasta), val(ref_name)
    
    output:
        path "*.vcf.gz", emit: filtered_vcfs
        path "*.vcf.gz.csi", emit: indices
        path "*.vcf.gz.stats", emit: stats
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def bcftools = "apptainer exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools"
        
        """
        for minfrac in 0.01 0.05 0.1 0.15 0.2 0.3 0.4; do
            ${bcftools} filter -i "FMT/AD[0:1] / FMT/DP[0] >= \${minfrac} | FMT/AD[0:2] / FMT/DP[0] >= \${minfrac} | FMT/AD[0:3] / FMT/DP[0] >= \${minfrac}" \\
                ${vcf} -Oz -o ${sample_id}-${ref_name}-m\${minfrac}.vcf.gz
            
            ${bcftools} index ${sample_id}-${ref_name}-m\${minfrac}.vcf.gz
            ${bcftools} stats ${sample_id}-${ref_name}-m\${minfrac}.vcf.gz > ${sample_id}-${ref_name}-m\${minfrac}.vcf.gz.stats
        done
        """
    } else {
        """
        echo "FILTER_VCF requires container with bcftools"
        exit 1
        """
    }
}
