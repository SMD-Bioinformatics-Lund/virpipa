process CREATE_CONSENSUS {
    tag { "${sample_id}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', pattern: "*.fasta*", enabled: params.publish_mode == 'debug'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf), path(reference_fasta)
        val min_freq
        path reference_metadata
    
    output:
        tuple val(run_name), val(sample_id), path("${sample_id}-${min_freq}-iupac.fasta"), path("${sample_id}-${min_freq}-iupac.fasta.fai"), emit: consensus
        tuple val(run_name), val(sample_id), path("${sample_id}-${min_freq}-iupac.fasta.consensus.json"), emit: metadata
        tuple val(run_name), val(sample_id), path("${sample_id}-majority.fasta"), path("${sample_id}-majority.fasta.fai"), path("${sample_id}-majority.fasta.consensus.json"), emit: analytical_majority
    
    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    
    def samtools = container_dir ? 
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"
    def python = container_dir ? "${container_runtime} exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" : 'python3'

    def ref_copy = reference_fasta.getName()
    def output_fasta = "${sample_id}-${min_freq}-iupac.fasta"
    def metadata_arg = reference_metadata ? "--reference-metadata ${reference_metadata}" : ''
    
    """
    set -euo pipefail

    if [[ "${reference_fasta}" != "${ref_copy}" ]]; then
        cp -L ${reference_fasta} ${ref_copy}
    fi

    zcat ${vcf} > input.vcf

    ${python} ${projectDir}/scripts/consensus_metadata.py --vcf input.vcf --reference ${ref_copy} ${metadata_arg} --output ${output_fasta} --identifier ${sample_id}-${min_freq}-iupac --mode iupac --min-af ${min_freq}
    ${python} ${projectDir}/scripts/consensus_metadata.py --vcf input.vcf --reference ${ref_copy} ${metadata_arg} --output ${sample_id}-majority.fasta --identifier ${sample_id}-majority --mode analytical_majority

    ${samtools} faidx ${output_fasta}
    ${samtools} faidx ${sample_id}-majority.fasta
    """
}
