#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBSAMPLE_READS } from '../modules/local/subsample/main'
include { REMOVE_HOSTILE } from '../modules/local/hostile/main'
include { MAP_READS } from '../modules/local/mapping/main'
include { ASSEMBLE_SPADES } from '../modules/local/assembly_spades/main'
include { ASSEMBLE_HYBRID } from '../modules/local/assembly_hybrid/main'
include { CREATE_CONSENSUS } from '../modules/local/consensus/main'
include { ANNOTATE_VADR } from '../modules/local/annotate_vadr/main'
include { SUBTYPE_BLAST } from '../modules/local/subtype/main'
include { ANNOTATE_RESISTANCE } from '../modules/local/resistance/main'
include { VARIANT_CALLING } from '../modules/local/variantcall/main'

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

    // Determine references: single genome or all in ref_dir
    def ch_references
    if (params.genome) {
        // Single genome specified
        def genome_file = file(params.genome)
        ch_references = Channel.value([tuple(genome_file.simpleName, genome_file)])
    } else if (params.ref_dir) {
        // Use all genomes in ref_dir
        def ref_dir = file(params.ref_dir)
        ch_references = Channel.fromPath("${ref_dir}/*.fa")
            .map { it -> [it.simpleName, it] }
    } else {
        error "Must specify either --genome or --ref_dir"
    }

    // Step 3: Map to each reference
    ch_mapped_all = ch_references.combine(ch_prepped).map { ref_name, ref_file, run_name, sample_id, read1, read2 ->
        [run_name, sample_id, read1, read2, ref_file, ref_name]
    }
    MAP_READS(ch_mapped_all)
    ch_mapped = MAP_READS.out.bams

    // Step 4: Spades assembly
    ASSEMBLE_SPADES(ch_prepped)
    ch_assembly = ASSEMBLE_SPADES.out.contigs

    // Step 5: Hybrid assembly with mummer
    // Map each assembly to each reference
    ch_hybrid_all = ch_assembly.combine(ch_references).map { run_name, sample_id, contigs, ref_name, ref_file ->
        [run_name, sample_id, contigs, ref_file, ref_name]
    }.combine(ch_references).map { run_name, sample_id, contigs, ref_file, ref_name, ref_file2, ref_name2 ->
        [run_name, sample_id, contigs, ref_file, ref_name]
    }
    
    ASSEMBLE_HYBRID(ch_hybrid_all)
    ch_hybrid = ASSEMBLE_HYBRID.out.hybrid_assembly

    // (Polishing - to be implemented)
    
    // Output
    ch_hybrid.view { "Final hybrid assemblies: $it" }
}
