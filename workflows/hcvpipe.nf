#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBSAMPLE_READS } from './modules/local/subsample/main'
include { REMOVE_HOSTILE } from './modules/local/hostile/main'
include { INDEX_REFERENCE } from './modules/local/index/main'
include { MAP_READS } from './modules/local/mapping/main'
include { ASSEMBLE_SPADES } from './modules/local/assembly_spades/main'
include { ASSEMBLE_HYBRID } from './modules/local/assembly_hybrid/main'
include { POLISH_PILON } from './modules/local/polish/main'
include { CREATE_CONSENSUS } from './modules/local/consensus/main'
include { ANNOTATE_VADR } from './modules/local/annotate_vadr/main'
include { SUBTYPE_BLAST } from './modules/local/subtype/main'
include { ANNOTATE_RESISTANCE } from './modules/local/resistance/main'

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
        REMOVE_HOSTILE(ch_subsampled, params.hostile_cache_dir)
        ch_prepped = REMOVE_HOSTILE.out.reads
    } else {
        ch_prepped = ch_subsampled
    }

    // Step 3: Get reference genomes
    Channel
        .fromPath("${params.ref_dir}/*.fa")
        .map { [it.baseName, it] }
        .set { ch_references }

    // Step 4: Map to each reference
    MAP_READS(ch_prepped, ch_references)
    ch_bams = MAP_READS.out.bams

    // Step 5: Assemble with SPAdes
    ASSEMBLE_SPADES(ch_prepped)
    ch_contigs = ASSEMBLE_SPADES.out.contigs

    // Step 6: Hybrid assembly
    ASSEMBLE_HYBRID(ch_contigs, ch_references)
    ch_hybrid = ASSEMBLE_HYBRID.out.hybrid_assembly

    // Step 7: Polish
    POLISH_PILON(ch_bams, ch_hybrid)
    ch_polished = POLISH_PILON.out.polished

    // Step 8: Create consensus (requires VCF - simplified for now)
    // This would need bcftools VCF generation first

    // Step 9: Annotate with VADR
    ANNOTATE_VADR(ch_polished)
    ch_vadr_gff = ANNOTATE_VADR.out.gff

    // Step 10: Subtype with BLAST
    SUBTYPE_BLAST(ch_polished)

    // Step 11: Annotate resistance (needs VCF, GFF, FASTA, subtype)
    // This would need VCF generation first
}
