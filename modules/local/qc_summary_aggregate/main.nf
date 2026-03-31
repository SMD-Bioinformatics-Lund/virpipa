process AGGREGATE_QC_SUMMARY {
    tag { run_name }
    label 'process_low'

    cpus 2
    memory '4 GB'
    time '30m'

    publishDir "${params.outdir}/${run_name}/pipeline_info", mode: 'copy', overwrite: true

    input:
        tuple val(run_name), path(qc_jsons)

    output:
        path "qc_summary.json", emit: json
        path "qc_summary.jsonl", emit: jsonl

    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def scripts_dir = params.scripts_dir ?: "${projectDir}/scripts"
    def python = container_dir ?
        "apptainer exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" :
        "python3"

    """
    set -euo pipefail

    ${python} ${scripts_dir}/aggregate_qc_summaries.py \\
        --output-json qc_summary.json \\
        --output-jsonl qc_summary.jsonl \\
        ${qc_jsons}
    """
}
