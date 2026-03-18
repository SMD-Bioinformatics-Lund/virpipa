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
include { CREATE_REPORT } from '../modules/local/report/main'

include { BAM2FASTA } from '../modules/local/bam2fasta/main'

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

    // Step 1: Remove human reads first (to match bash pipeline)
    if (params.remove_human) {
        def cache_dir = params.hostile_cache_dir ?: ''
        REMOVE_HOSTILE(ch_samples, cache_dir)
        ch_prepped = REMOVE_HOSTILE.out.reads
    } else {
        ch_prepped = ch_samples
    }

    // Step 2: Subsample reads for pilon polishing and reference selection
    // If subsample_reads is set, use that for pilon (matching bash pipeline)
    // The 250k subsample for reference selection is taken from the same source
    if (params.subsample_reads) {
        SUBSAMPLE_READS(ch_prepped, params.subsample_reads)
        ch_pilon_reads = SUBSAMPLE_READS.out.reads
        // Also use the same subsampled reads for reference selection (close enough to 250k)
        ch_subsampled = ch_pilon_reads
    } else {
        ch_pilon_reads = ch_prepped
        // For reference selection, subsample 250k
        SUBSAMPLE_READS(ch_prepped, 250000)
        ch_subsampled = SUBSAMPLE_READS.out.reads
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
    ch_mapped_all = ch_references.combine(ch_subsampled).map { ref_name, ref_file, run_name, sample_id, read1, read2 ->
        [run_name, sample_id, read1, read2, ref_file, ref_name]
    }
    MAP_READS(ch_mapped_all)
    ch_mapped = MAP_READS.out.bams
    ch_stats = MAP_READS.out.stats

    // Step 4: Spades assembly
    ASSEMBLE_SPADES(ch_subsampled)
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
    // Use subsampled reads for pilon (params.subsample_reads) - matching bash pipeline
    ch_polish_input = ch_hybrid.cross(ch_pilon_reads)
        .filter { hybrid, reads -> hybrid[1] == reads[1] }
        .map { hybrid, reads ->
            tuple(hybrid[0], hybrid[1], reads[2], reads[3], hybrid[2])
        }
    
    POLISH_PILON_LOOP(ch_polish_input)
    ch_polished = POLISH_PILON_LOOP.out.polished
    ch_pilon_bam_with_index = POLISH_PILON_LOOP.out.final_bam_with_index
    
    // Step 6b: Regenerate consensus from pilon BAM (matching bash pipeline)
    // This uses bam2fasta to create a consensus at 100% IUPAC from the pilon BAM
    // Need to join with best_ref to get the actual ref_name (not just sample_id)
    // Use sample_id as join key, keep run_name separate
    ch_pilon_for_bam2fasta = ch_pilon_bam_with_index.map { run_name, sample_id, bam, bai ->
        [sample_id, run_name, bam, bai]
    }
    ch_polished_for_bam2fasta = ch_polished.map { run_name, sample_id, fasta, fai ->
        [sample_id, fasta, fai]
    }
    ch_best_ref_for_bam2fasta = ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta ->
        [sample_id, ref_name]
    }
    
    ch_bam2fasta_input = ch_pilon_for_bam2fasta
        .join(ch_polished_for_bam2fasta)
        .join(ch_best_ref_for_bam2fasta)
        .map { sample_id, run_name, bam, bai, fasta, fai, ref_name ->
            tuple(run_name, sample_id, bam, bai, fasta, ref_name)
        }
    
    BAM2FASTA(ch_bam2fasta_input, "1.0")
    ch_pilon_regenerated = BAM2FASTA.out.fasta
    
    // Step 6c: Create CRAM from polished BAM
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
    
    // Step 7: Variant calling on pilon-polished BAM (matching bash pipeline)
    // Use pilon BAM and bam2fasta-regenerated FASTA as reference
    ch_pilon_for_varcall = ch_pilon_bam_with_index.map { run_name, sample_id, bam, bai ->
        [sample_id, run_name, bam, bai]
    }
    
    // Get bam2fasta regenerated pilon fasta for variant calling
    ch_polished_for_varcall = ch_pilon_regenerated.map { run_name, sample_id, fasta, fai ->
        [sample_id, run_name, fasta]
    }
    
    // Join pilon BAM with regenerated pilon fasta
    // Also join with best reference to get the ref name
    ch_variant_input = ch_pilon_for_varcall
        .map { it -> [it[0], it] }  // key by sample_id
        .join(ch_polished_for_varcall.map { it -> [it[0], it] })
        .join(ch_best_ref_with_name.map { it -> [it[1], it] })  // key by sample_id
        .map { sample_id, pilon, polished, best_ref ->
            // best_ref = [run_name, sample_id, ref_name, fasta_file]
            tuple(pilon[1], pilon[0], pilon[2], pilon[3], polished[2], best_ref[2])
        }
    
    VARIANT_CALLING(ch_variant_input)
    ch_vcf = VARIANT_CALLING.out.vcf
    
    // Save copy for resistance before using
    ch_vcf_for_resistance = ch_vcf

    // Step 7b: Filter VCF at multiple min fractions
    ch_filter_input = ch_vcf.cross(ch_best_ref_with_name)
        .filter { vcf, best_ref -> vcf[1] == best_ref[1] }
        .map { vcf, best_ref ->
            tuple(vcf[0], vcf[1], vcf[2], vcf[3], best_ref[3], best_ref[2])
        }
    
    FILTER_VCF(ch_filter_input)

    // Step 8: Create consensus from VCF - use bam2fasta output (100% IUPAC) as reference
    // This matches bash pipeline which replaces pilon FASTA with bam2fasta output
    ch_vcf_for_consensus = ch_vcf.map { run_name, sample_id, vcf, vcf_idx -> [sample_id, run_name, vcf, vcf_idx] }
    ch_polished_for_consensus = ch_pilon_regenerated.map { run_name, sample_id, fasta, fai -> [sample_id, run_name, fasta, fai] }
    
    // ch_vcf_for_consensus: [sample_id, run_name, vcf, vcf_idx]
    // ch_polished_for_consensus: [sample_id, run_name, fasta, fai]
    ch_consensus_input = ch_vcf_for_consensus.cross(ch_polished_for_consensus)
        .filter { vcf_data, polished_data -> vcf_data[0] == polished_data[0] }
        .map { vcf_data, polished_data ->
            tuple(vcf_data[1], vcf_data[0], vcf_data[2], polished_data[2], polished_data[3])
        }
    
    CREATE_CONSENSUS(ch_consensus_input, "0.15")
    // Get consensus with sample metadata - use consensus output which emits *iupac.fasta files
    // Need to match back the consensus files to sample metadata
    ch_consensus_keyed = ch_consensus_input.map { run_name, sample_id, vcf, fasta, fai ->
        [sample_id, run_name, sample_id]
    }
    ch_consensus_files = CREATE_CONSENSUS.out.consensus.map { fasta ->
        def fname = fasta.baseName.replaceAll('-0.15-iupac$', '')
        [fname, fasta]
    }
    ch_consensus_with_meta = ch_consensus_keyed.cross(ch_consensus_files)
        .filter { key, fasta -> key[2] == fasta[0] }
        .map { key, fasta ->
            tuple(key[0], key[1], fasta[1])
        }
    
    // Save copy for resistance annotation (channels can only be used once)
    ch_consensus_for_resistance = ch_consensus_with_meta

    // Step 9: Subtype with BLAST
    // Get blast db from params or use default (directory containing hcvgluerefs)
    def blast_db_path = params.blast_db ? file(params.blast_db) : file("${params.ref_dir}/hcvglue")
    
    // Get subtype from BLAST (for report generation)
    // Note: SUBTYPE_BLAST outputs blast files, subtype extraction needs to be done separately
    // For now, use placeholder - can be enhanced later to parse subtype from BLAST results
    ch_subtype_for_report = Channel.of(['unknown'])
    
    // Step 9b: Create report files (matching bash pipeline)
    // Get the VCF stats for the unfiltered pilon VCF (from VARIANT_CALLING)
    ch_vcf_stats = ch_vcf.map { run_name, sample_id, vcf, tbi -> 
        [sample_id, run_name, vcf, tbi] 
    }
    
    // Get best ref name
    ch_best_ref_info = ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta ->
        [sample_id, ref_name]
    }
    
    // Run stats on the main VCF (not filtered)
    def container_dir = params.container_dir ?: ''
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        process STATS_VCF {
            tag { "${sample_id}" }
            label 'process_low'
            
            input:
                tuple val(run_name), val(sample_id), path(vcf), path(tbi), val(ref_name)
            
            output:
                tuple val(sample_id), path("*.vcf.gz.stats")
            
            script:
            def bcftools = "apptainer exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools"
            """
            ${bcftools} stats \${vcf} > \${sample_id}.vcf.gz.stats
            """
        }
        
        ch_vcf_stats = ch_vcf
            .combine(ch_best_ref_info, by: 1)
            .map { run_name, sample_id, vcf, tbi, ref_name ->
                tuple(run_name, sample_id, vcf, tbi, ref_name)
            }
            | STATS_VCF
    }
    
    // Build report input: combine cram, ref_fasta, vcf_stats, subtype
    // Format: (run_name, sample_id, vcf_stats, cram, crai, fasta, subtype)
    ch_report_input = ch_cram_output
        .map { run_name, sample_id, cram, crai -> [sample_id, run_name, cram, crai] }
        .combine(ch_pilon_regenerated.map { run_name, sample_id, fasta, fai -> [sample_id, fasta] }, by: 0)
        .combine(ch_vcf_stats, by: 0)
        .map { sample_id, run_name, cram, crai, fasta, vcf_stats ->
            tuple(run_name, sample_id, vcf_stats, cram, crai, fasta)
        }
        .combine(ch_subtype_for_report.ifEmpty(['unknown']), by: 0)
        .map { run_name, sample_id, vcf_stats, cram, crai, fasta, subtype ->
            tuple(run_name, sample_id, vcf_stats, cram, crai, fasta, subtype)
        }
    
    CREATE_REPORT(ch_report_input)
    
    // BLAST for subtype (for downstream analysis)
    ch_consensus_for_blast = ch_consensus_with_meta
    
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
    
    // Save GFF for resistance before using
    ch_vadr_gff_for_resistance = ch_vadr_gff
    
    // Step 11: Annotate resistance (needs VCF + GFF + fasta + subtype + rules)
    rules_path = params.resistance_rules ?: "${projectDir}/assets/hbv_result_rules.csv"
    rules_csv = file(rules_path)
    
    // Combine VCF with GFF by sample_id
    ch_vcf_gff = ch_vcf_for_resistance.cross(ch_vadr_gff_for_resistance)
        .filter { vcf, gff -> vcf[1] == gff[1] }
        .map { vcf, gff -> [vcf[0], vcf[1], vcf[2], gff[2]] }
    
    // Add consensus fasta 
    ch_resistance_input = ch_vcf_gff.combine(ch_consensus_for_resistance)
        .map { run_name, sample_id, vcf, gff, cons_run, cons_sample, cons_fasta ->
            [run_name, sample_id, vcf, gff, cons_fasta]
        }
    
    // Add placeholder subtype
    ch_resistance_full = ch_resistance_input.map { run_name, sample_id, vcf, gff, fasta ->
        [tuple(run_name, sample_id, vcf, gff, fasta), '1a']
    }
    
    ANNOTATE_RESISTANCE(ch_resistance_full.map { it[0] }, ch_resistance_full.map { it[1] }, rules_csv)
    
    // Output final results
    // ch_consensus_with_meta.view { "Final consensus: $it" }
}
