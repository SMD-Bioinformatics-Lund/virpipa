#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HCVPIPE } from './workflows/hcvpipe'

workflow {
    main:
    HCVPIPE()

    onComplete:
    if (params.completion_log_dir?.toString()?.trim()) {
        def source = (params.csv ?: params.input ?: workflow.runName).toString()
        def base = source ? file(source).baseName : workflow.runName
        def logFile = file("${params.completion_log_dir}/${base}.complete")
        logFile.parent.toFile().mkdirs()

        def msg = """\
            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            scriptFile  : ${workflow.scriptFile}
            workDir     : ${workflow.workDir}
            csv         : ${params.csv ?: params.input}
            exit status : ${workflow.exitStatus}
            errorMessage: ${workflow.errorMessage}
            errorReport :
            ${workflow.errorReport ?: ''}
            """
            .stripIndent()

        logFile.text = msg
    }
}
