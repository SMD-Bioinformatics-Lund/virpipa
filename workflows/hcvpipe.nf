#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBSAMPLE_READS } from '../modules/local/subsample/main'
include { REMOVE_HOSTILE } from '../modules/local/hostile/main'
include { MAP_READS } from '../modules/local/mapping/main'
include { SELECT_BEST_REFERENCE } from '../modules/local/bestref/main'
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
    ch_stats = MAP_READS.out.stats

    // Step 4: Spades assembly
    ASSEMBLE_SPADES(ch_prepped)
    ch_assembly = ASSEMBLE_SPADES.out.contigs

    // Step 5: Select best reference based on mapping stats
    // Group stats by sample
    ch_stats_per_sample = ch_stats.groupTuple(by: [0, 1])
    
    // Get ref_dir - convert to absolute path to handle both relative and absolute
    def ref_dir = file(params.ref_dir).absolutePath
    
    // Add ref_dir to each entry as part of the tuple  
    ch_best_ref_input = ch_stats_per_sample.map { run_name, sample_id, stats_list ->
        tuple(run_name, sample_id, stats_list, ref_dir)
    }
    
    SELECT_BEST_REFERENCE(ch_best_ref_input)
    ch_best_ref = SELECT_BEST_REFERENCE.out.best_ref

    // Step 5b: Hybrid assembly with best reference
    // Combine assembly with best ref
    ch_hybrid_input = ch_assembly.combine(ch_best_ref).map { run_name, sample_id, contigs, best_ref_name, best_ref_file ->
        [ [run_name, sample_id, contigs], best_ref_file, best_ref_name ]
    }
    
    // Split the tuple into 3 separate channels for the process
    ch_hybrid_contigs = ch_hybrid_input.map { it[0] }
    ch_hybrid_ref = ch_hybrid_input.map { it[1] }
    ch_hybrid_refname = ch_hybrid_input.map { it[2] }
    
    ASSEMBLE_HYBRID(ch_hybrid_contigs, ch_hybrid_ref, ch_hybrid_refname)
    ch_hybrid = ASSEMBLE_HYBRID.out.hybrid_assembly

    // (Polishing loop - 10 iterations - to be implemented)
    
    // Step 7: Variant calling on mapped reads - use ONLY best reference
    // Filter mapped bams to only keep those matching best reference
    ch_mapped_best = ch_mapped.join(ch_best_ref).map { run_name, sample_id, bam, bai, best_ref_name, best_ref_file ->
        [run_name, sample_id, bam, bai, best_ref_file, best_ref_name]
    }
    
    VARIANT_CALLING(ch_mapped_best)
    ch_vcf = VARIANT_CALLING.out.vcf

    // Step 8: Create consensus from VCF - only for best reference
    ch_consensus_input = ch_vcf.map { run_name, sample_id, vcf, vcf_idx ->
        // We need the best ref file - join with best_ref
        [run_name, sample_id, vcf]
    }.join(ch_best_ref).map { run_name, sample_id, vcf, best_ref_name, best_ref_file ->
        [run_name, sample_id, vcf, best_ref_file, file(best_ref_file.toString() + '.fai')]
    }
    
    CREATE_CONSENSUS(ch_consensus_input, "0.15")
    // Get consensus with sample metadata
    ch_consensus_with_meta = ch_consensus_input.map { run_name, sample_id, vcf, ref_file, fai ->
        [run_name, sample_id]
    }.combine(CREATE_CONSENSUS.out.fasta).map { run_name, sample_id, fasta ->
        [run_name, sample_id, fasta]
    }

    // Step 9: Subtype with BLAST
    // Get blast db from params or use default
    def blast_db = params.blast_db ? file(params.blast_db) : file("${params.ref_dir}/hcvglue/hcvgluerefs")
    ch_subtype_tuple = ch_consensus_with_meta.map { run_name, sample_id, fasta ->
        [ [run_name, sample_id, fasta], blast_db ]
    }
    ch_subtype_fasta = ch_subtype_tuple.map { it[0] }
    ch_subtype_db = ch_subtype_tuple.map { it[1] }
    
    SUBTYPE_BLAST(ch_subtype_fasta, ch_subtype_db)

    // Step 10: Annotate with VADR
    def vadr_model = params.vadr_model ?: 'vadr-models-flavi'
    ch_vadr_tuple = ch_consensus_with_meta.map { run_name, sample_id, fasta ->
        [ [run_name, sample_id, fasta], vadr_model ]
    }
    ch_vadr_fasta = ch_vadr_tuple.map { it[0] }
    ch_vadr_model = ch_vadr_tuple.map { it[1] }
    
    ANNOTATE_VADR(ch_vadr_fasta, ch_vadr_model)

    // Step 11: Annotate resistance (needs VCF + GFF + fasta + subtype)
    // For now, skip as it requires more complex input handling
    
    // Output final results
    // ch_consensus_with_meta.view { "Final consensus: $it" }
}
