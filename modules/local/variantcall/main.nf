process VARIANT_CALLING {
    tag { "${sample_id}:${vcf_prefix}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/vcf", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(bam), path(bai), path(ref_fasta), val(vcf_prefix)
    
    output:
        tuple val(run_name), val(sample_id), path("${sample_id}-${vcf_prefix}.vcf.gz"), path("${sample_id}-${vcf_prefix}.vcf.gz.csi"), emit: vcf
        tuple val(run_name), val(sample_id), path("${sample_id}-${vcf_prefix}.vcf.gz.stats"), emit: stats
        tuple val(run_name), val(sample_id), path("${sample_id}-${vcf_prefix}.bcf"), emit: bcf
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    def bcftools = container_dir ?
        "apptainer exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools" :
        "bcftools"

    def ref_copy = ref_fasta.getName()
    def bcf_name = "${sample_id}-${vcf_prefix}.bcf"
    def vcf_name = "${sample_id}-${vcf_prefix}.vcf.gz"
    def stats_name = "${sample_id}-${vcf_prefix}.vcf.gz.stats"

    """
    set -euo pipefail

    if [[ "${ref_fasta}" != "${ref_copy}" ]]; then
        cp -L ${ref_fasta} ${ref_copy}
    fi

    ${bcftools} mpileup -Ob -f ${ref_copy} -d 1000000 -a AD,DP -o ${bcf_name} ${bam}
    ${bcftools} call -Oz -m -A --ploidy 1 -o ${vcf_name} ${bcf_name}
    ${bcftools} index ${vcf_name}
    ${bcftools} stats ${vcf_name} > ${stats_name}
    """
}
