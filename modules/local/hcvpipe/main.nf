process RUN_HCVPIPE {
    tag { "${run_name}:${sample_id}" }
    label 'process_high'

    cpus params.cpus
    memory params.memory
    time params.time

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(run_name), val(sample_id), val(read1), val(read2), val(lid)

    output:
    path("${run_name}/${sample_id}"), emit: sample_dirs

    script:
    def lidArg = lid ? "-l '${lid}'" : ''
    def remove_human = params.remove_human instanceof Boolean ? params.remove_human : params.remove_human?.toString()?.toBoolean()
    def subsample_reads = params.subsample_reads?.toString()?.trim() ? params.subsample_reads.toString().toInteger() : 0
    def hostileArg = remove_human ? '' : '-H'
    def containerArg = params.container_dir ? "--container-dir '${params.container_dir}'" : ''
    def bindArg = params.bind_paths ? "--bind-paths '${params.bind_paths}'" : ''
    def hostileCacheArg = params.hostile_cache_dir ? "--hostile-cache-dir '${params.hostile_cache_dir}'" : ''
    def readArgs = read2 ? "'${read1}' '${read2}'" : "'${read1}'"

    """
    work_root=\$(pwd)
    mkdir -p "\${work_root}/${run_name}"
    bash '${params.scripts_dir}/hcvpipe.sh' \\
      --scripts-dir '${params.scripts_dir}' \\
      --ref-dir '${params.ref_dir}' \\
      ${containerArg} \\
      ${bindArg} \\
      ${hostileCacheArg} \\
      -o "\${work_root}/${run_name}" \\
      --outname '${sample_id}' \\
      -s ${subsample_reads} \\
      -c ${task.cpus} \\
      ${lidArg} \\
      ${hostileArg} \\
      ${readArgs}
    """
}
