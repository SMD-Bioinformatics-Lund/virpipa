process BAM2FASTA {
    tag { "${sample_id}:${ref_name}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/fasta", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(bam), path(bai), path(ref_fasta), val(ref_name)
        val ambiguity_threshold
    
    output:
        tuple val(run_name), val(sample_id), path("${sample_id}-${ref_name}.fasta"), path("${sample_id}-${ref_name}.fasta.fai"), emit: fasta
        tuple val(run_name), val(sample_id), path("${sample_id}.fasta"), path("${sample_id}.fasta.fai"), emit: replacement_fasta
        path "${sample_id}-${ref_name}.vcf.gz", emit: vcf
        path "${sample_id}-${ref_name}.vcf.gz.csi", emit: vcf_index
        path "${sample_id}-${ref_name}.vcf.gz.stats", emit: stats
        path "${sample_id}-${ref_name}.bcf", emit: bcf
        tuple val(run_name), val(sample_id), path("${sample_id}-${ref_name}.vcf.gz"), emit: vcf_with_meta
        tuple val(run_name), val(sample_id), path("${sample_id}-${ref_name}.vcf.gz.csi"), emit: vcf_index_with_meta
        tuple val(run_name), val(sample_id), path("${sample_id}-${ref_name}.vcf.gz.stats"), emit: stats_with_meta
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    
    def bcftools = container_dir ? 
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools" :
        "bcftools"
    
    def samtools = container_dir ? 
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"
    
    def ref_copy = ref_fasta.getName()
    def bcf_name = "${sample_id}-${ref_name}.bcf"
    def vcf_name = "${sample_id}-${ref_name}.vcf.gz"
    def fasta_name = "${sample_id}-${ref_name}.fasta"
    def stats_name = "${sample_id}-${ref_name}.vcf.gz.stats"
    def replacement_fasta = "${sample_id}.fasta"

    """
    set -euo pipefail

    if [[ "${ref_fasta}" != "${ref_copy}" ]]; then
        cp -L ${ref_fasta} ${ref_copy}
    fi

    ${bcftools} mpileup -Ob -f ${ref_copy} -d 1000000 -a AD,DP -o ${bcf_name} ${bam}
    ${bcftools} call -Oz -m -A --ploidy 1 -o ${vcf_name} ${bcf_name}
    ${bcftools} index ${vcf_name}
    ${bcftools} stats ${vcf_name} > ${stats_name}

    zcat ${vcf_name} > input.vcf

    awk -v MIN_AF=${ambiguity_threshold} -v MIN_DP=7 -f ${projectDir}/scripts/vcf_to_iupac.awk input.vcf ${ref_copy} > ${fasta_name}

    sed -i "s/>.*/>${sample_id}-${ref_name}/" ${fasta_name}
    ${samtools} faidx ${fasta_name}

    cp ${fasta_name} ${replacement_fasta}
    sed -i "s/>.*/>${sample_id}/" ${replacement_fasta}
    ${samtools} faidx ${replacement_fasta}
    """
}
