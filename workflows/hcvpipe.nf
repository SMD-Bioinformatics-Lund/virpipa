include { RUN_HCVPIPE } from '../modules/local/hcvpipe/main'

workflow HCVPIPE {
    if (!params.input) {
        error "Missing required parameter: --input <samplesheet.csv>"
    }

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample = (row.clarity_sample_id ?: row.sample ?: row.id ?: '').toString().trim()
            if (!sample) {
                error "Samplesheet row is missing sample id (expected clarity_sample_id, sample, or id)"
            }

            def read1 = (row.read1 ?: row.fastq_1 ?: '').toString().trim()
            if (!read1) {
                error "Samplesheet row for sample '${sample}' is missing read1/fastq_1"
            }

            def read2 = (row.read2 ?: row.fastq_2 ?: '').toString().trim()
            def lid = (row.sample_name ?: row.lid ?: '').toString().trim()
            tuple(sample, read1, read2, lid)
        }
        .set { ch_samples }

    RUN_HCVPIPE(ch_samples)

    emit:
    RUN_HCVPIPE.out.sample_dirs
}
