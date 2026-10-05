process BUILD_QC_SUMMARY {
    tag { sample_id }
    label 'process_low'

    cpus 2
    memory '4 GB'
    time '30m'

    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', overwrite: true

    input:
        tuple val(run_name), val(sample_id), val(lid), path(results_dir)
        val sample_info_json_path

    output:
        path "*_qc_summary.json", emit: json
        tuple val(run_name), val(sample_id), path("*_qc_summary.json"), emit: json_with_meta

    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    def scripts_dir = params.scripts_dir != 'scripts' ? params.scripts_dir : (active_profiles.contains('hpc') ? '/fs1/jonas/src/virpipa/scripts' : "${projectDir}/scripts")
    def python = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" :
        "python3"

    def sampleInfoArg = sample_info_json_path ? "--sample-info '${sample_info_json_path}'" : ''

    """
    set -euo pipefail

    ${python} ${scripts_dir}/build_qc_summary.py \\
        --results-dir ${results_dir} \\
        --run-name '${run_name}' \\
        --sample-id '${sample_id}' \\
        --lid '${lid}' \\
        ${sampleInfoArg} \\
        --output ${sample_id}_qc_summary.json
    """
}
