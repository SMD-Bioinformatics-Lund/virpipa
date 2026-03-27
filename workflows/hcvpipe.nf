#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBSAMPLE_READS } from '../modules/local/subsample/main'
include { REMOVE_HOSTILE } from '../modules/local/hostile/main'
include { MAP_READS } from '../modules/local/mapping/main'
include { MAP_READS_NOOPT } from '../modules/local/mapping_noopt/main'
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
include { CREATE_CRAM as CREATE_CRAM_BESTREF } from '../modules/local/cram/main'
include { CREATE_CRAM as CREATE_CRAM_PILON } from '../modules/local/cram/main'
include { CREATE_CRAM as CREATE_CRAM_IUPAC } from '../modules/local/cram/main'
include { LOG_COVERAGE } from '../modules/local/coverage/main'
include { CREATE_REPORT as CREATE_REPORT_BESTREF } from '../modules/local/report/main'
include { CREATE_REPORT as CREATE_REPORT_IUPAC } from '../modules/local/report/main'

include { BAM2FASTA as BAM2FASTA_BESTREF } from '../modules/local/bam2fasta/main'
include { BAM2FASTA as BAM2FASTA_PILON } from '../modules/local/bam2fasta/main'

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
    ch_mapped_with_ref = MAP_READS.out.bams_with_ref

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

    // Step 5c: Build the chosen best-reference consensus/VCF track from the
    // initial reference-mapping BAM. This feeds the `${sample}-${subtype}`
    // CRAM/report outputs that the bash pipeline keeps in `results/`.
    ch_best_ref_mapped = ch_mapped_with_ref
        .map { run_name, sample_id, ref_name, bam, bai ->
            [[sample_id, ref_name], [run_name, sample_id, ref_name, bam, bai]]
        }
        .join(ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta ->
            [[sample_id, ref_name], [run_name, sample_id, ref_name, fasta]]
        })
        .map { key, mapped, best_ref ->
            tuple(mapped[0], mapped[1], mapped[3], mapped[4], best_ref[3], mapped[2])
        }

    BAM2FASTA_BESTREF(ch_best_ref_mapped, "0.01")
    ch_best_ref_fasta = BAM2FASTA_BESTREF.out.fasta

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
    
    ch_bam2fasta_input = ch_pilon_for_bam2fasta
        .join(ch_polished_for_bam2fasta)
        .map { sample_id, run_name, bam, bai, fasta, fai ->
            tuple(run_name, sample_id, bam, bai, fasta, '1.0-iupac')
        }
    
    BAM2FASTA_PILON(ch_bam2fasta_input, "1.0")
    ch_pilon_regenerated = BAM2FASTA_PILON.out.replacement_fasta
    
    // Step 6c: Create CRAM from the regenerated pilon replacement FASTA.
    // This matches the bash pipeline, which replaces pilon/${sample}.fasta
    // with the 1.0-iupac majority-call FASTA before creating CRAM outputs.
    ch_regenerated_simple = ch_pilon_regenerated.map { run_name, sample_id, fasta, fai -> [sample_id, run_name, fasta] }
    ch_pilon_simple = ch_pilon_bam_with_index.map { run_name, sample_id, bam, bai -> [sample_id, run_name, bam, bai] }
    
    ch_cram_input = ch_regenerated_simple.cross(ch_pilon_simple)
        .filter { regenerated, bam -> regenerated[0] == bam[0] }
        .map { regenerated, bam ->
            def fasta_abs = regenerated[2].toAbsolutePath()
            [regenerated[1], regenerated[0], bam[2], bam[3], fasta_abs, regenerated[0]]
        }
    
    CREATE_CRAM_PILON(ch_cram_input)
    ch_cram_output = CREATE_CRAM_PILON.out.cram_with_index

    // Step 6c-b: Create the best-reference CRAM from the original reference
    // mapping outputs and the original chosen reference FASTA.
    ch_best_ref_cram_input = ch_mapped_with_ref
        .map { run_name, sample_id, ref_name, bam, bai ->
            [[sample_id, ref_name], [run_name, sample_id, ref_name, bam, bai]]
        }
        .join(ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta ->
            [[sample_id, ref_name], [run_name, sample_id, ref_name, fasta]]
        })
        .map { key, mapped, best_ref ->
            def ref_abs = file(best_ref[3]).toAbsolutePath()
            tuple(mapped[0], mapped[1], mapped[3], mapped[4], ref_abs, "${mapped[1]}-${mapped[2]}")
        }

    CREATE_CRAM_BESTREF(ch_best_ref_cram_input)
    ch_best_ref_cram_output = CREATE_CRAM_BESTREF.out.cram_with_index
    
    // Step 6d: Log coverage from CRAM using the regenerated replacement FASTA
    ch_coverage_input = ch_cram_output
        .map { run_name, sample_id, cram, crai -> 
            [sample_id, run_name, cram, crai]
        }
        .join(ch_pilon_regenerated.map { run_name, sample_id, fasta, fai -> [sample_id, fasta] })
        .map { sample_id, run_name, cram, crai, fasta ->
            def fasta_abs = file(fasta).toAbsolutePath()
            tuple(run_name, sample_id, cram, crai, fasta_abs)
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
    
    // Join pilon BAM with regenerated pilon replacement fasta
    ch_variant_input = ch_pilon_for_varcall
        .map { it -> [it[0], it] }  // key by sample_id
        .join(ch_polished_for_varcall.map { it -> [it[0], it] })
        .map { sample_id, pilon, polished ->
            tuple(pilon[1], pilon[0], pilon[2], pilon[3], polished[2], 'pilon')
        }
    
    VARIANT_CALLING(ch_variant_input)
    ch_vcf = VARIANT_CALLING.out.vcf
    
    // Save copy for resistance before using
    ch_vcf_for_resistance = ch_vcf

    // Step 7b: Filter VCF at multiple min fractions
    ch_filter_input = ch_vcf.map { run_name, sample_id, vcf, vcf_idx ->
            tuple(run_name, sample_id, vcf, vcf_idx, "${sample_id}-pilon")
        }
    
    FILTER_VCF(ch_filter_input)

    // Step 8: Create consensus from VCF - use bam2fasta output (100% IUPAC) as reference
    // This matches bash pipeline which replaces pilon FASTA with bam2fasta output
    ch_vcf_for_consensus = ch_vcf.map { run_name, sample_id, vcf, vcf_idx -> [sample_id, run_name, vcf, vcf_idx] }
    ch_polished_for_consensus = ch_pilon_regenerated.map { run_name, sample_id, fasta, fai -> [sample_id, run_name, fasta] }
    
    // ch_vcf_for_consensus: [sample_id, run_name, vcf, vcf_idx]
    // ch_polished_for_consensus: [sample_id, run_name, fasta]
    ch_consensus_input = ch_vcf_for_consensus.cross(ch_polished_for_consensus)
        .filter { vcf_data, polished_data -> vcf_data[0] == polished_data[0] }
        .map { vcf_data, polished_data ->
            tuple(vcf_data[1], vcf_data[0], vcf_data[2], polished_data[2])
        }
    
    CREATE_CONSENSUS(ch_consensus_input, "0.15")
    ch_consensus_with_meta = CREATE_CONSENSUS.out.consensus.map { run_name, sample_id, fasta, fai ->
        tuple(run_name, sample_id, fasta)
    }
    
    // Save copy for resistance annotation (channels can only be used once)
    ch_consensus_for_resistance = ch_consensus_with_meta

    // Step 8b: Map the 0.15-iupac consensus with the no-opt sentieon path
    ch_iupac_mapping_input = ch_consensus_with_meta.cross(ch_pilon_reads)
        .filter { consensus, reads -> consensus[1] == reads[1] }
        .map { consensus, reads ->
            tuple(consensus[0], consensus[1], reads[2], reads[3], consensus[2], '0.15-iupac')
        }

    MAP_READS_NOOPT(ch_iupac_mapping_input)
    ch_iupac_bams = MAP_READS_NOOPT.out.bams

    // Step 8c: Create the 0.15-iupac CRAM from the no-opt mapped BAM
    ch_iupac_cram_input = ch_iupac_bams
        .map { run_name, sample_id, bam, bai -> [sample_id, run_name, bam, bai] }
        .join(ch_consensus_with_meta.map { run_name, sample_id, fasta -> [sample_id, fasta] })
        .map { sample_id, run_name, bam, bai, fasta ->
            def fasta_abs = file(fasta).toAbsolutePath()
            tuple(run_name, sample_id, bam, bai, fasta_abs, "${sample_id}-0.15-iupac")
        }

    CREATE_CRAM_IUPAC(ch_iupac_cram_input)
    ch_iupac_cram_output = CREATE_CRAM_IUPAC.out.cram_with_index

    // Step 8d: Create the 0.15-iupac report from the filtered 0.15 VCF stats
    ch_iupac_report_input = FILTER_VCF.out.stats
        .filter { run_name, sample_id, stats -> stats.getName() == "${sample_id}-pilon-m0.15.vcf.gz.stats" }
        .map { run_name, sample_id, stats -> [sample_id, run_name, stats] }
        .join(ch_iupac_cram_output.map { run_name, sample_id, cram, crai -> [sample_id, run_name, cram, crai] })
        .join(ch_consensus_with_meta.map { run_name, sample_id, fasta -> [sample_id, fasta] })
        .join(ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta -> [sample_id, ref_name] })
        .map { sample_id, run_name, stats, cram, crai, fasta, ref_name ->
            tuple(run_name, sample_id, stats, cram, crai, fasta, ref_name, "${sample_id}-0.15-iupac")
        }

    CREATE_REPORT_IUPAC(ch_iupac_report_input)

    // Step 8e: Create the `${sample}-${subtype}` report from the chosen
    // best-reference BAM2FASTA outputs and CRAM.
    ch_best_ref_report_input = BAM2FASTA_BESTREF.out.stats
        .filter { stats -> !stats.getName().contains('1.0-iupac') }
        .map { stats ->
            def base = stats.getName().replaceFirst(/\.vcf\.gz\.stats$/, '')
            def parts = base.split('-', 2)
            tuple(parts[0], base, stats)
        }
        .join(ch_best_ref_cram_output.map { run_name, sample_id, cram, crai -> [sample_id, [run_name, cram, crai]] })
        .join(ch_best_ref_fasta.map { run_name, sample_id, fasta, fai ->
            def suffix = fasta.getBaseName().replaceFirst("^${sample_id}-", '')
            [sample_id, [suffix, fasta]]
        })
        .join(ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta -> [sample_id, ref_name] })
        .filter { sample_id, report_base, stats, cram_data, fasta_data, ref_name ->
            report_base == "${sample_id}-${ref_name}" && fasta_data[0] == ref_name
        }
        .map { sample_id, report_base, stats, cram_data, fasta_data, ref_name ->
            tuple(cram_data[0], sample_id, stats, cram_data[1], cram_data[2], fasta_data[1], ref_name, report_base)
        }

    CREATE_REPORT_BESTREF(ch_best_ref_report_input)

    // Step 9: Subtype with BLAST
    // Get blast db from params or use default (directory containing hcvgluerefs)
    def blast_db_path = params.blast_db ? file(params.blast_db) : file("${params.ref_dir}/hcvglue")
    
    // Get subtype from BLAST (for report generation)
    // Note: SUBTYPE_BLAST outputs blast files, subtype extraction needs to be done separately
    // For now, use placeholder - can be enhanced later to parse subtype from BLAST results
    ch_subtype_for_report = Channel.of(['unknown'])
    
    // Step 9b: Create report files (matching bash pipeline) - SKIPPING FOR NOW
    // Get the VCF stats for the unfiltered pilon VCF (from VARIANT_CALLING)
    // This step can be added back after debugging
    
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
    ch_vadr_input = BAM2FASTA_PILON.out.replacement_fasta.map { run_name, sample_id, fasta, fai ->
        tuple(run_name, sample_id, fasta)
    }
    ch_vadr_tuple = ch_vadr_input.map { run_name, sample_id, fasta ->
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
