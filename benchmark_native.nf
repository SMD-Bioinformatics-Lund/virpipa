#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process NATIVE_CONSENSUS_EXPERIMENT {
    tag { "${run_name}:${sample_id}" }
    cpus 8
    memory '8 GB'
    time '2h'
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy'
    input:
        tuple val(run_name), val(sample_id), path(scaffold), path(cram), path(crai), path(iupac), path(majority), path(resistance)
    output:
        tuple val(run_name), val(sample_id), path('native'), emit: experiments
    script:
    def containers = params.container_dir ?: "${projectDir}/assets/containers"
    def binds = params.bind_paths
    def runtime = params.container_runtime ?: 'apptainer'
    """
    export VIRPIPA_BENCHMARK_RUNTIME='${runtime}'
    export VIRPIPA_BENCHMARK_BINDS='${binds}'
    bash ${projectDir}/scripts/benchmark_native_sample.sh \\
        '${scaffold}' '${cram}' '${iupac}' '${majority}' '${sample_id}' \\
        '${containers}' '${projectDir}/scripts' '${task.cpus}' '${resistance}'
    """
}

workflow {
    if (!params.input || !params.outdir) {
        error 'Supply --input <native-inputs.tsv> and an isolated --outdir'
    }
    inputs = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.run_name, row.sample_id, file(row.scaffold, checkIfExists: true),
                file(row.scaffold_cram, checkIfExists: true), file(row.scaffold_index, checkIfExists: true),
                file(row.iupac, checkIfExists: true), file(row.majority, checkIfExists: true), file(row.resistance, checkIfExists: true))
        }
    NATIVE_CONSENSUS_EXPERIMENT(inputs)
}
