process RUN_HCVPIPE {
    tag { sample_id }
    label 'process_high'

    cpus params.cpus
    memory params.memory
    time params.time

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), val(read1), val(read2), val(lid)

    output:
    path("${sample_id}"), emit: sample_dirs

    script:
    def lidArg = lid ? "-l '${lid}'" : ''
    def hostileArg = params.remove_human ? '' : '-H'
    def containerArg = params.container_dir ? "--container-dir '${params.container_dir}'" : ''
    def bindArg = params.bind_paths ? "--bind-paths '${params.bind_paths}'" : ''
    def hostileCacheArg = params.hostile_cache_dir ? "--hostile-cache-dir '${params.hostile_cache_dir}'" : ''
    def readArgs = read2 ? "'${read1}' '${read2}'" : "'${read1}'"

    """
    bash '${params.scripts_dir}/hcvpipe.sh' \\
      --scripts-dir '${params.scripts_dir}' \\
      --ref-dir '${params.ref_dir}' \\
      ${containerArg} \\
      ${bindArg} \\
      ${hostileCacheArg} \\
      -o . \\
      --outname '${sample_id}' \\
      -s ${params.subsample_reads} \\
      -c ${task.cpus} \\
      ${lidArg} \\
      ${hostileArg} \\
      ${readArgs}
    """
}
