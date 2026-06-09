process WRITE_VIRTITTA_IMPORT_MARKER {
    tag { run_name }
    label 'process_low'

    publishDir { params.virtitta_import_dir?.toString()?.trim() ?: '.' }, mode: 'copy', overwrite: true, enabled: params.virtitta_import_dir != null && params.virtitta_import_dir.toString().trim() != ''

    input:
        tuple val(run_name), val(run_dir), path(qc_json), path(qc_jsonl)

    output:
        path "${run_name}.sqlimport", emit: marker

    when:
        params.virtitta_import_dir != null && params.virtitta_import_dir.toString().trim() != ''

    script:
    def normalized_run_dir = run_dir.endsWith('/') ? run_dir : "${run_dir}/"
    """
    set -euo pipefail

    printf '%s\n' '--run-dir ${normalized_run_dir}' > "${run_name}.sqlimport"
    """
}
