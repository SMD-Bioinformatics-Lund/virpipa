#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def coerceIntegerParam(def raw) {
    def text = raw?.toString()?.trim()
    if (text && !text.equalsIgnoreCase('false')) {
        return text.toInteger()
    }
    return 0
}

include { SUBSAMPLE_READS } from '../modules/local/subsample/main'

workflow {
    if (params.partition && params.queue && params.partition.toString() != params.queue.toString()) {
        error "Provide either --partition or --queue for the SLURM partition, not both with different values"
    }

    if (params.input && params.csv && params.input.toString() != params.csv.toString()) {
        error "Provide either --input or --csv for the samplesheet, not both with different values"
    }

    def samplesheet_param = params.input ?: params.csv
    if (!samplesheet_param) {
        error "Missing required parameter: --input <samplesheet.csv> (alias: --csv)"
    }

    channel
        .fromPath(samplesheet_param, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample = (row.clarity_sample_id ?: row.sample ?: row.id ?: '').toString().trim()
            if (!sample) {
                error "Samplesheet row is missing sample id"
            }

            def read1 = (row.read1 ?: row.fastq_1 ?: '').toString().trim()
            if (!read1) {
                error "Samplesheet row for sample '${sample}' is missing read1/fastq_1"
            }

            def read2 = (row.read2 ?: row.fastq_2 ?: '').toString().trim()
            def run_name = (row.run_name ?: params.run_name ?: 'test').toString().trim()
            
            // Return tuple with sample_id and reads
            return [run_name, sample, file(read1), read2 ? file(read2) : []]
        }
        .set { ch_samples }

    // Pass nreads as separate parameter
    SUBSAMPLE_READS(ch_samples, coerceIntegerParam(params.subsample_reads))
}
