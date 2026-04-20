process FILTER_VCF {
    tag { "${sample_id}" }
    label 'process_medium'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/vcf", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf), path(vcf_index), val(vcf_prefix)
    
    output:
        tuple val(run_name), val(sample_id), path("${vcf_prefix}-m*.vcf.gz"), emit: filtered_vcfs
        tuple val(run_name), val(sample_id), path("${vcf_prefix}-m*.vcf.gz.csi"), emit: indices
        tuple val(run_name), val(sample_id), path("${vcf_prefix}-m*.vcf.gz.stats"), emit: stats
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'

    def bcftools = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools" :
        "bcftools"

    """
    set -euo pipefail

    for minfrac in 0.01 0.05 0.1 0.15 0.2 0.3 0.4; do
        ${bcftools} filter -i "FMT/AD[0:1] / FMT/DP[0] >= \${minfrac} | FMT/AD[0:2] / FMT/DP[0] >= \${minfrac} | FMT/AD[0:3] / FMT/DP[0] >= \${minfrac}" \\
            ${vcf} -Oz -o ${vcf_prefix}-m\${minfrac}.vcf.gz

        ${bcftools} index ${vcf_prefix}-m\${minfrac}.vcf.gz
        ${bcftools} stats ${vcf_prefix}-m\${minfrac}.vcf.gz > ${vcf_prefix}-m\${minfrac}.vcf.gz.stats
    done
    """
}
