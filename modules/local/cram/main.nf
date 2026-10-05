process CREATE_CRAM {
    tag { "${sample_id}:${output_name}" }
    label 'process_medium'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', enabled: params.publish_mode == 'debug'
    
    input:
        tuple val(run_name), val(sample_id), path(bam), path(bai), path(ref_fasta), val(output_name)
    
    output:
        tuple val(run_name), val(sample_id), path("${output_name}.cram"), path("${output_name}.cram.crai"), emit: cram_with_index
    
    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'

    def samtools = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"

    def ref_copy = ref_fasta.getName()

    """
    set -euo pipefail

    if [[ "${ref_fasta}" != "${ref_copy}" ]]; then
        cp -L ${ref_fasta} ${ref_copy}
    fi

    ${samtools} view -O cram,embed_ref -T ${ref_copy} ${bam} -o ${output_name}.cram
    ${samtools} index ${output_name}.cram
    """
}
