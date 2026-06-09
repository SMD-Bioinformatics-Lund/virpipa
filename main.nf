#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HCVPIPE } from './workflows/hcvpipe'

workflow {
    main:
    HCVPIPE()

    onComplete:
    if (params.completion_log_dir?.toString()?.trim()) {
        def runName = params.run_name?.toString()?.trim()
        def source = (params.csv ?: params.input)?.toString()

        if (!runName && source) {
            def sourceFile = file(source)
            if (sourceFile.exists()) {
                def rows = sourceFile.readLines().findAll { it.trim() && !it.trim().startsWith('#') }
                if (rows.size() > 1) {
                    def headers = rows[0].split(',', -1).collect { it.trim() }
                    def values = rows[1].split(',', -1).collect { it.trim() }
                    def runNameIndex = headers.findIndexOf { it == 'run_name' || it == 'sequencing_run' }
                    if (runNameIndex >= 0 && runNameIndex < values.size()) {
                        runName = values[runNameIndex]
                    }
                }
            }
        }

        if (!runName) {
            runName = workflow.runName
        }

        def logFile = file("${params.completion_log_dir}/${runName}-HCV.complete")
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
