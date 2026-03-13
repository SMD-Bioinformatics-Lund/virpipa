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
    if (params.genome) {
        def genome_file = file(params.genome)
        ch_references = Channel.value(tuple(genome_file.simpleName, genome_file))
    } else if (params.ref_dir) {
        def ref_dir = file(params.ref_dir)
        ch_references = Channel.fromPath("${ref_dir}/*.fa")
            .map { tuple(it.simpleName, it) }
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
    // Map each assembly to each reference - build hybrid for EACH ref
    // Use combine to create cartesian product of assembly x references
    ch_hybrid_tuple = ch_assembly.combine(ch_references).map { run_name, sample_id, contigs, ref_name, ref_file ->
        [ [run_name, sample_id, contigs], ref_file, ref_name ]
    }
    
    // Split the tuple into 3 separate channels for the process
    ch_hybrid_contigs = ch_hybrid_tuple.map { it[0] }
    ch_hybrid_ref = ch_hybrid_tuple.map { it[1] }
    ch_hybrid_refname = ch_hybrid_tuple.map { it[2] }
    
    ASSEMBLE_HYBRID(ch_hybrid_contigs, ch_hybrid_ref, ch_hybrid_refname)
    ch_hybrid = ASSEMBLE_HYBRID.out.hybrid_assembly

    // Step 6: Select best reference based on mapping stats (use first for now)
    // In production, would parse stats to find lowest error rate
    // ch_references is a value channel, so we can just use it directly
    
    // (Polishing loop - 10 iterations - to be implemented)
    
    // Step 7: Variant calling on mapped reads (using original mapping to refs)
    // Combine bams with references
    ch_variant_input = ch_mapped.combine(ch_references).map { run_name, sample_id, bam, bai, ref_name, ref_file ->
        [run_name, sample_id, bam, bai, ref_file, ref_name]
    }
    
    VARIANT_CALLING(ch_variant_input)
    ch_vcf = VARIANT_CALLING.out.vcf

    // Step 8: Create consensus from VCF
    // ch_vcf has: (run_name, sample_id, vcf, vcf_idx)
    // ch_references has: (ref_name, ref_file)
    // combine creates: (run_name, sample_id, vcf, vcf_idx, ref_name, ref_file)
    ch_consensus_input = ch_vcf.combine(ch_references).map { run_name, sample_id, vcf, vcf_idx, ref_name, ref_file ->
        [run_name, sample_id, vcf, ref_file, file(ref_file.toString() + '.fai')]
    }
    
    CREATE_CONSENSUS(ch_consensus_input, "0.15")
    ch_consensus = CREATE_CONSENSUS.out.fasta

    // Step 9: Subtype with BLAST
    // Get blast db from params or use default
    def blast_db = params.blast_db ? file(params.blast_db) : file("${params.ref_dir}/hcvglue/hcvgluerefs")
    ch_subtype_tuple = ch_consensus.map { run_name, sample_id, fasta ->
        [ [run_name, sample_id, fasta], blast_db ]
    }
    ch_subtype_fasta = ch_subtype_tuple.map { it[0] }
    ch_subtype_db = ch_subtype_tuple.map { it[1] }
    
    SUBTYPE_BLAST(ch_subtype_fasta, ch_subtype_db)

    // Step 10: Annotate with VADR
    def vadr_model = params.vadr_model ?: 'vadr-models-flavi'
    ch_vadr_tuple = ch_consensus.map { run_name, sample_id, fasta ->
        [ [run_name, sample_id, fasta], vadr_model ]
    }
    ch_vadr_fasta = ch_vadr_tuple.map { it[0] }
    ch_vadr_model = ch_vadr_tuple.map { it[1] }
    
    ANNOTATE_VADR(ch_vadr_fasta, ch_vadr_model)

    // Step 11: Annotate resistance (needs VCF + GFF + fasta + subtype)
    // For now, skip as it requires more complex input handling
    
    // Output final results
    ch_consensus.view { "Final consensus: $it" }
}
