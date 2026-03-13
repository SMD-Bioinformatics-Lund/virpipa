#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBSAMPLE_READS } from '../modules/local/subsample/main'
include { REMOVE_HOSTILE } from '../modules/local/hostile/main'
include { MAP_READS } from '../modules/local/mapping/main'

workflow HCVPIPE {
    if (!params.input) {
        error "Missing required parameter: --input <samplesheet.csv>"
    }

    // Channel: read samplesheet
    Channel
        .fromPath(params.input, checkIfExists: true)
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
            
            return [run_name, sample, file(read1), read2 ? file(read2) : []]
        }
        .set { ch_samples }

    // Step 1: Subsample reads
    SUBSAMPLE_READS(ch_samples, params.subsample_reads)
    ch_subsampled = SUBSAMPLE_READS.out.reads

    // Step 2: Remove human reads (optional)
    if (params.remove_human) {
        def cache_dir = params.hostile_cache_dir ?: ''
        REMOVE_HOSTILE(ch_subsampled, cache_dir)
        ch_prepped = REMOVE_HOSTILE.out.reads
    } else {
        ch_prepped = ch_subsampled
    }

    // Step 3: Map reads
    if (params.genome) {
        def genome_file = file(params.genome)
        def genome_name = genome_file.simpleName
        // Add genome to each sample tuple
        ch_for_mapping = ch_prepped.map { run_name, sample_id, read1, read2 ->
            [run_name, sample_id, read1, read2, genome_file, genome_name]
        }
        MAP_READS(ch_for_mapping)
        ch_mapped = MAP_READS.out.bams
    } else {
        ch_mapped = ch_prepped
    }

    // Output - just publish the preprocessed reads
    ch_mapped.view { "Final bams: $it" }
}
