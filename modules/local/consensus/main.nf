process CREATE_CONSENSUS {
    tag { "${sample_id}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy', pattern: "*.fasta"
    
    input:
        tuple val(run_name), val(sample_id), path(vcf), path(reference_fasta)
        val min_freq
    
    output:
        tuple val(run_name), val(sample_id), path("${sample_id}-${min_freq}-iupac.fasta"), path("${sample_id}-${min_freq}-iupac.fasta.fai"), emit: consensus
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    
    def samtools = container_dir ? 
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"

    def ref_copy = reference_fasta.getName()
    def output_fasta = "${sample_id}-${min_freq}-iupac.fasta"
    
    """
    set -euo pipefail

    if [[ "${reference_fasta}" != "${ref_copy}" ]]; then
        cp -L ${reference_fasta} ${ref_copy}
    fi

    zcat ${vcf} > input.vcf

    awk -v MIN_AF=${min_freq} -v MIN_DP=7 -f ${projectDir}/scripts/vcf_to_iupac.awk input.vcf ${ref_copy} > ${output_fasta}

    sed -i "s/>.*/>${sample_id}-${min_freq}-iupac/" ${output_fasta}

    ${samtools} faidx ${output_fasta}
    """
}
