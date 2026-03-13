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

include { POLISH_PILON_LOOP } from '../modules/local/polish/main'
include { FILTER_VCF } from '../modules/local/filter_vcf/main'
include { CREATE_CRAM } from '../modules/local/cram/main'
include { LOG_COVERAGE } from '../modules/local/coverage/main'

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
    def ref_dir = file(params.ref_dir).toAbsolutePath()
    
    // Add ref_dir to each entry as part of the tuple  
    ch_best_ref_input = ch_stats_per_sample.map { run_name, sample_id, stats_list ->
        tuple(run_name, sample_id, stats_list, ref_dir)
    }
    
    SELECT_BEST_REFERENCE(ch_best_ref_input)
    ch_best_ref = SELECT_BEST_REFERENCE.out.best_ref

    // Extract ref name from the fasta filename (e.g., 3a-D17763.fa -> 3a-D17763)
    ch_best_ref_with_name = ch_best_ref.map { run_name, sample_id, fasta_file ->
        def ref_name = fasta_file.baseName
        tuple(run_name, sample_id, ref_name, fasta_file)
    }

    // Step 5b: Hybrid assembly with best reference
    // Use cross to combine and then filter by matching sample
    ch_hybrid_input = ch_assembly.cross(ch_best_ref_with_name)
        .filter { assembly, best_ref -> assembly[1] == best_ref[1] }
        .map { assembly, best_ref -> 
            tuple(tuple(assembly[0], assembly[1], assembly[2]), best_ref[3], best_ref[2])
        }
    
    // Split the tuple into 3 separate channels for the process
    ch_hybrid_contigs = ch_hybrid_input.map { it[0] }
    ch_hybrid_ref = ch_hybrid_input.map { it[1] }
    ch_hybrid_refname = ch_hybrid_input.map { it[2] }
    
    ASSEMBLE_HYBRID(ch_hybrid_contigs, ch_hybrid_ref, ch_hybrid_refname)
    ch_hybrid = ASSEMBLE_HYBRID.out.hybrid_assembly

    // Step 6: Polishing loop (10 iterations with convergence check)
    // Prepare input: combine reads with hybrid assembly
    ch_polish_input = ch_hybrid.cross(ch_prepped)
        .filter { hybrid, reads -> hybrid[1] == reads[1] }
        .map { hybrid, reads ->
            tuple(hybrid[0], hybrid[1], reads[2], reads[3], hybrid[2])
        }
    
    POLISH_PILON_LOOP(ch_polish_input)
    ch_polished = POLISH_PILON_LOOP.out.polished
    ch_pilon_bam_with_index = POLISH_PILON_LOOP.out.final_bam_with_index
    
    // Step 6b: Create CRAM from polished BAM
    // Use cross and filter since join is not working well
    ch_polished_simple = ch_polished.map { run_name, sample_id, fasta, fai -> [sample_id, run_name, fasta] }
    ch_pilon_simple = ch_pilon_bam_with_index.map { run_name, sample_id, bam, bai -> [sample_id, run_name, bam, bai] }
    
    ch_cram_input = ch_polished_simple.cross(ch_pilon_simple)
        .filter { polish, bam -> polish[0] == bam[0] }
        .map { polish, bam ->
            def fasta_abs = polish[2].toAbsolutePath()
            [polish[1], polish[0], bam[2], bam[3], fasta_abs]
        }
    
    CREATE_CRAM(ch_cram_input, "pilon")
    ch_cram_output = CREATE_CRAM.out.cram_with_index
    
    // Step 6c: Log coverage from CRAM - use polished fasta as reference
    ch_cram_out_simple = ch_cram_output
        .map { run_name, sample_id, cram, crai -> [sample_id, run_name, cram, crai] }
    
    ch_coverage_input = ch_cram_out_simple.cross(ch_polished_simple)
        .filter { cram, polish -> cram[0] == polish[0] }
        .map { cram, polish ->
            def fasta_abs = polish[2].toAbsolutePath()
            [cram[1], cram[0], cram[2], cram[3], fasta_abs]
        }
    
    LOG_COVERAGE(ch_coverage_input)
    
    // Step 7: Variant calling on mapped reads - use ONLY best reference
    // Cross and filter by sample_id
    ch_mapped_best = ch_mapped.cross(ch_best_ref_with_name)
        .filter { bam, best_ref -> bam[1] == best_ref[1] }
        .map { bam, best_ref ->
            tuple(bam[0], bam[1], bam[2], bam[3], best_ref[3], best_ref[2])
        }
    
    VARIANT_CALLING(ch_mapped_best)
    ch_vcf = VARIANT_CALLING.out.vcf

    // Step 7b: Filter VCF at multiple min fractions
    ch_filter_input = ch_vcf.cross(ch_best_ref_with_name)
        .filter { vcf, best_ref -> vcf[1] == best_ref[1] }
        .map { vcf, best_ref ->
            tuple(vcf[0], vcf[1], vcf[2], vcf[3], best_ref[3], best_ref[2])
        }
    
    FILTER_VCF(ch_filter_input)

    // Step 8: Create consensus from VCF - only for best reference
    ch_consensus_input = ch_vcf.cross(ch_best_ref_with_name)
        .filter { vcf, best_ref -> vcf[1] == best_ref[1] }
        .map { vcf, best_ref ->
            tuple(vcf[0], vcf[1], vcf[2], best_ref[3], file(best_ref[3].toString() + '.fai'))
        }
    
    CREATE_CONSENSUS(ch_consensus_input, "0.15")
    // Get consensus with sample metadata
    ch_consensus_with_meta = ch_consensus_input.map { run_name, sample_id, vcf, ref_file, fai ->
        tuple(run_name, sample_id)
    }.combine(CREATE_CONSENSUS.out.fasta).map { run_name, sample_id, fasta ->
        tuple(run_name, sample_id, fasta)
    }

    // Step 9: Subtype with BLAST
    // Get blast db from params or use default (directory containing hcvgluerefs)
    def blast_db_path = params.blast_db ? file(params.blast_db) : file("${params.ref_dir}/hcvglue")
    
    if (blast_db_path.exists()) {
        ch_subtype_tuple = ch_consensus_with_meta.map { run_name, sample_id, fasta ->
            [ tuple(run_name, sample_id, fasta), blast_db_path ]
        }
        ch_subtype_fasta = ch_subtype_tuple.map { it[0] }
        ch_subtype_db = ch_subtype_tuple.map { it[1] }
        
        SUBTYPE_BLAST(ch_subtype_fasta, ch_subtype_db)
    } else {
        println "WARNING: BLAST database not found at ${blast_db_path}, skipping SUBTYPE_BLAST"
    }

    // Step 10: Annotate with VADR
    def vadr_model = params.vadr_model ?: 'vadr-models-flavi'
    ch_vadr_tuple = ch_consensus_with_meta.map { run_name, sample_id, fasta ->
        [ [run_name, sample_id, fasta], vadr_model ]
    }
    ch_vadr_fasta = ch_vadr_tuple.map { it[0] }
    ch_vadr_model = ch_vadr_tuple.map { it[1] }
    
    ANNOTATE_VADR(ch_vadr_fasta, ch_vadr_model)
    ch_vadr_gff = ANNOTATE_VADR.out.gff
    
    // Step 11: Annotate resistance (needs VCF + GFF + fasta + subtype + rules)
    // Use hbv_result_rules.csv from assets as default
    def rules_csv = params.resistance_rules ? file(params.resistance_rules) : file("${projectDir}/assets/hbv_result_rules.csv")
    
    if (rules_csv.exists()) {
        // Get subtype from BLAST results - for now use a placeholder or extract from consensus header
        // TODO: Connect BLAST subtype to resistance annotation
        ch_resistance_input = ch_vcf.cross(ch_vadr_gff)
            .filter { vcf, gff -> vcf[1] == gff[1] }
            .map { vcf, gff ->
                tuple(vcf[0], vcf[1], vcf[2], gff[2])
            }
        
        // Add consensus fasta and placeholder subtype
        ch_resistance_with_fasta = ch_resistance_input.combine(ch_consensus_with_meta)
            .map { run_name, sample_id, vcf, gff, consensus ->
                tuple(run_name, sample_id, vcf, gff, consensus)
            }
        
        // For now, use a default subtype - in production this would come from BLAST
        ch_resistance_full = ch_resistance_with_fasta.map { run_name, sample_id, vcf, gff, fasta ->
            tuple(tuple(run_name, sample_id, vcf, gff, fasta), '1a')
        }
        
        ch_resistance_vcf = ch_resistance_full.map { it[0] }
        ch_resistance_subtype = ch_resistance_full.map { it[1] }
        
        ANNOTATE_RESISTANCE(ch_resistance_vcf, ch_resistance_subtype, rules_csv)
    } else {
        println "WARNING: Resistance rules not found at ${rules_csv}, skipping ANNOTATE_RESISTANCE"
    }
    
    // Output final results
    // ch_consensus_with_meta.view { "Final consensus: $it" }
}
