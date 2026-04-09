#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def resolvePathFromBase(String rawPath, def baseDir) {
    def pathText = rawPath?.toString()?.trim()
    if (!pathText) {
        return null
    }

    def candidate = new File(pathText)
    if (candidate.isAbsolute()) {
        return file(candidate)
    }

    return file(new File(baseDir.toString(), pathText))
}

include { SUBSAMPLE_READS } from '../modules/local/subsample/main'
include { SUBSAMPLE_READS as SUBSAMPLE_READS_REFSEL } from '../modules/local/subsample/main'
include { REMOVE_HOSTILE } from '../modules/local/hostile/main'
include { MAP_READS } from '../modules/local/mapping/main'
include { MAP_READS_NOOPT } from '../modules/local/mapping_noopt/main'
include { SELECT_BEST_REFERENCE } from '../modules/local/bestref/main'
include { ASSEMBLE_SPADES } from '../modules/local/assembly_spades/main'
include { ASSEMBLE_HYBRID } from '../modules/local/assembly_hybrid/main'
include { CREATE_CONSENSUS } from '../modules/local/consensus/main'
include { ANNOTATE_VADR } from '../modules/local/annotate_vadr/main'
include { ANNOTATE_RESISTANCE } from '../modules/local/resistance/main'
include { VARIANT_CALLING } from '../modules/local/variantcall/main'
include { BUILD_QC_SUMMARY } from '../modules/local/qc_summary/main'
include { AGGREGATE_QC_SUMMARY } from '../modules/local/qc_summary_aggregate/main'

include { POLISH_PILON_LOOP } from '../modules/local/polish/main'
include { FILTER_VCF } from '../modules/local/filter_vcf/main'
include { CREATE_CRAM as CREATE_CRAM_BESTREF } from '../modules/local/cram/main'
include { CREATE_CRAM as CREATE_CRAM_PILON } from '../modules/local/cram/main'
include { CREATE_CRAM as CREATE_CRAM_IUPAC } from '../modules/local/cram/main'
include { LOG_COVERAGE } from '../modules/local/coverage/main'
include { CREATE_REPORT as CREATE_REPORT_BESTREF } from '../modules/local/report/main'
include { CREATE_REPORT as CREATE_REPORT_IUPAC } from '../modules/local/report/main'
include { FINALIZE_RESULTS } from '../modules/local/finalize_results/main'

include { BAM2FASTA as BAM2FASTA_BESTREF } from '../modules/local/bam2fasta/main'
include { BAM2FASTA as BAM2FASTA_PILON } from '../modules/local/bam2fasta/main'
include { SUBTYPE_BLAST as SUBTYPE_BLAST_MAIN } from '../modules/local/subtype/main'
include { SUBTYPE_BLAST as SUBTYPE_BLAST_IUPAC } from '../modules/local/subtype/main'
include { SUBTYPE_BLAST as SUBTYPE_BLAST_PILON_IUPAC } from '../modules/local/subtype/main'

workflow HCVPIPE {
    if (params.input && params.csv && params.input.toString() != params.csv.toString()) {
        error "Provide either --input or --csv for the samplesheet, not both with different values"
    }

    def samplesheet_param = params.input ?: params.csv
    if (!samplesheet_param) {
        error "Missing required parameter: --input <samplesheet.csv> (alias: --csv)"
    }

    def samplesheet_path = file(samplesheet_param).toAbsolutePath()
    def samplesheet_dir = samplesheet_path.parent
    def samplesheet_parent = samplesheet_dir?.parent
    def samplesheet_grandparent = samplesheet_parent?.parent
    def project_root = projectDir.toString()
    def resolved_ref_dir = params.ref_dir ? resolvePathFromBase(params.ref_dir, project_root) : null
    def resolved_genome = params.genome ? resolvePathFromBase(params.genome, project_root) : null
    def resolved_blast_db = params.blast_db ? resolvePathFromBase(params.blast_db, project_root) : null
    def resolved_rules_path = params.resistance_rules
        ? resolvePathFromBase(params.resistance_rules, project_root)
        : resolvePathFromBase('assets/hcv_geno2pheno_rules.csv', project_root)
    def expected_reference_count = resolved_genome
        ? 1
        : new File(resolved_ref_dir.toString()).listFiles()?.count { it.isFile() && it.name.endsWith('.fa') } ?: 0
    def resolved_sample_info_json = [
        file("${samplesheet_dir}/clarity_sample_info.json"),
        samplesheet_parent ? file("${samplesheet_parent}/clarity_sample_info.json") : null,
        samplesheet_grandparent ? file("${samplesheet_grandparent}/clarity_sample_info.json") : null
    ].find { it && it.exists() }

    // Channel: read samplesheet
    Channel
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
            def run_name = (row.run_name ?: row.sequencing_run ?: params.run_name ?: 'test').toString().trim()
            def lid = (row.sample_name ?: row.lid ?: '').toString().trim()

            def resolved_read1 = resolvePathFromBase(read1, samplesheet_dir)
            def resolved_read2 = read2 ? resolvePathFromBase(read2, samplesheet_dir) : []

            return [run_name, sample, lid, resolved_read1, resolved_read2]
        }
        .set { ch_sample_rows }

    ch_samples = ch_sample_rows.map { run_name, sample_id, lid, read1, read2 ->
        [run_name, sample_id, read1, read2]
    }

    ch_sample_lids = ch_sample_rows.map { run_name, sample_id, lid, read1, read2 ->
        [sample_id, [run_name, lid]]
    }

    // Step 1: Remove human reads first (to match bash pipeline)
    if (params.remove_human) {
        def cache_dir = params.hostile_cache_dir ?: ''
        REMOVE_HOSTILE(ch_samples, cache_dir)
        ch_prepped = REMOVE_HOSTILE.out.reads
        ch_hostile_json = REMOVE_HOSTILE.out.hostile_json_with_meta.map { run_name, sample_id, hostile_json ->
            [sample_id, hostile_json.toString()]
        }
    } else {
        ch_prepped = ch_samples
        ch_hostile_json = ch_sample_rows.map { run_name, sample_id, lid, read1, read2 ->
            [sample_id, '']
        }
    }

    // Step 2: Subsample reads for pilon polishing and reference selection
    // If subsample_reads is set, use that for pilon (matching bash pipeline)
    // Reference selection should never use more than 250k reads
    if (params.subsample_reads) {
        SUBSAMPLE_READS(ch_prepped, params.subsample_reads)
        ch_pilon_reads = SUBSAMPLE_READS.out.reads
        if (params.subsample_reads <= 250000) {
            ch_subsampled = ch_pilon_reads
        } else {
            SUBSAMPLE_READS_REFSEL(ch_prepped, 250000)
            ch_subsampled = SUBSAMPLE_READS_REFSEL.out.reads
        }
    } else {
        ch_pilon_reads = ch_prepped
        // For reference selection, subsample 250k
        SUBSAMPLE_READS_REFSEL(ch_prepped, 250000)
        ch_subsampled = SUBSAMPLE_READS_REFSEL.out.reads
    }

    // Determine references: single genome or all in ref_dir
    if (resolved_genome) {
        def genome_file = resolved_genome
        ch_references = Channel.value(tuple(genome_file.simpleName, genome_file))
    } else if (resolved_ref_dir) {
        def ref_dir = resolved_ref_dir
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
    ch_stats_per_sample = ch_stats.groupTuple(by: [0, 1], size: expected_reference_count)
    
    // Get ref_dir - convert to absolute path to handle both relative and absolute
    def ref_dir = resolved_ref_dir.toAbsolutePath()
    
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
    // Key explicitly by sample_id so batch runs keep all samples instead of
    // relying on cross/filter fan-in behavior.
    ch_hybrid_input = ch_assembly
        .map { run_name, sample_id, contigs ->
            [sample_id, [run_name, sample_id, contigs]]
        }
        .join(ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta ->
            [sample_id, [ref_name, fasta]]
        })
        .map { sample_id, assembly, best_ref ->
            tuple(tuple(assembly[0], assembly[1], assembly[2]), best_ref[1], best_ref[0])
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

    // Match the bash pipeline's initial best-reference bam2fasta call, which
    // uses the default majority-call threshold rather than the later quasispecies
    // thresholds used for the pilon-derived tracks.
    BAM2FASTA_BESTREF(ch_best_ref_mapped, "1.0")
    ch_best_ref_fasta = BAM2FASTA_BESTREF.out.fasta

    // Step 6: Polishing loop (10 iterations with convergence check)
    // Prepare input: combine reads with hybrid assembly
    // Use subsampled reads for pilon (params.subsample_reads) - matching bash pipeline
    ch_polish_input = ch_hybrid
        .map { run_name, sample_id, hybrid_fasta ->
            [sample_id, [run_name, sample_id, hybrid_fasta]]
        }
        .join(ch_pilon_reads.map { run_name, sample_id, read1, read2 ->
            [sample_id, [read1, read2]]
        })
        .map { sample_id, hybrid, reads ->
            tuple(hybrid[0], hybrid[1], reads[0], reads[1], hybrid[2])
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
    
    ch_cram_input = ch_regenerated_simple
        .join(ch_pilon_simple)
        .map { sample_id, regenerated_run_name, fasta, pilon_run_name, bam, bai ->
            def fasta_abs = fasta.toAbsolutePath()
            [regenerated_run_name, sample_id, bam, bai, fasta_abs, sample_id]
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
    
    def extractSubtypeFromBlast = { blast_file ->
        def blast_line = blast_file.readLines().find { line ->
            line && !line.startsWith('query acc.ver')
        }
        if (!blast_line) {
            error "Could not derive subtype from BLAST output: ${blast_file}"
        }

        def fields = blast_line.split('\t')
        if (fields.size() < 2 || !fields[1]) {
            error "Malformed BLAST output while deriving subtype: ${blast_file}"
        }

        def subtype = fields[1].split('_')[0]
        if (!subtype) {
            error "Could not parse subtype from BLAST subject '${fields[1]}' in ${blast_file}"
        }

        return subtype
    }

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
    ch_consensus_input = ch_vcf_for_consensus
        .join(ch_polished_for_consensus)
        .map { sample_id, run_name, vcf, vcf_idx, polished_run_name, fasta ->
            tuple(run_name, sample_id, vcf, fasta)
        }
    
    CREATE_CONSENSUS(ch_consensus_input, "0.15")
    ch_consensus_with_meta = CREATE_CONSENSUS.out.consensus.map { run_name, sample_id, fasta, fai ->
        tuple(run_name, sample_id, fasta)
    }
    
    // Step 8b: Map the 0.15-iupac consensus with the no-opt sentieon path
    ch_iupac_mapping_input = ch_consensus_with_meta
        .map { run_name, sample_id, fasta ->
            [sample_id, [run_name, sample_id, fasta]]
        }
        .join(ch_pilon_reads.map { run_name, sample_id, read1, read2 ->
            [sample_id, [read1, read2]]
        })
        .map { sample_id, consensus, reads ->
            tuple(consensus[0], consensus[1], reads[0], reads[1], consensus[2], '0.15-iupac')
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
        .map { run_name, sample_id, stats_files ->
            def stats = stats_files.find { it.getName() == "${sample_id}-pilon-m0.15.vcf.gz.stats" }
            stats ? tuple(sample_id, stats) : null
        }
        .filter { sample_id, stats -> stats != null }
        .join(ch_iupac_cram_output.map { run_name, sample_id, cram, crai -> [sample_id, [run_name, cram, crai]] })
        .join(ch_consensus_with_meta.map { run_name, sample_id, fasta -> [sample_id, fasta] })
        .join(ch_best_ref_with_name.map { run_name, sample_id, ref_name, fasta -> [sample_id, ref_name] })
        .map { sample_id, stats, cram_data, fasta, ref_name ->
            tuple(cram_data[0], sample_id, stats, cram_data[1], cram_data[2], fasta, ref_name, "${sample_id}-0.15-iupac")
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
    def blast_db_path = resolved_blast_db ?: file("${resolved_ref_dir}/hcvglue")
    
    // Run BLAST on the main polished FASTA, the 0.15-iupac FASTA, and the pilon-iupac FASTA.
    ch_main_blast_with_meta = Channel.empty()
    ch_iupac_blast_with_meta = Channel.empty()
    ch_pilon_iupac_blast_with_meta = Channel.empty()

    if (blast_db_path.exists()) {
        ch_main_blast_tuple = ch_pilon_regenerated.map { run_name, sample_id, fasta, fai ->
            [ tuple(run_name, sample_id, fasta), blast_db_path ]
        }
        SUBTYPE_BLAST_MAIN(ch_main_blast_tuple.map { it[0] }, ch_main_blast_tuple.map { it[1] })
        ch_main_blast_with_meta = SUBTYPE_BLAST_MAIN.out.blast_with_meta

        ch_iupac_blast_tuple = ch_consensus_with_meta.map { run_name, sample_id, fasta ->
            [ tuple(run_name, sample_id, fasta), blast_db_path ]
        }
        SUBTYPE_BLAST_IUPAC(ch_iupac_blast_tuple.map { it[0] }, ch_iupac_blast_tuple.map { it[1] })
        ch_iupac_blast_with_meta = SUBTYPE_BLAST_IUPAC.out.blast_with_meta

        ch_pilon_iupac_blast_tuple = POLISH_PILON_LOOP.out.polished_pilon_iupac.map { run_name, sample_id, fasta ->
            [ tuple(run_name, sample_id, fasta), blast_db_path ]
        }
        SUBTYPE_BLAST_PILON_IUPAC(ch_pilon_iupac_blast_tuple.map { it[0] }, ch_pilon_iupac_blast_tuple.map { it[1] })
        ch_pilon_iupac_blast_with_meta = SUBTYPE_BLAST_PILON_IUPAC.out.blast_with_meta
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
    
    // Step 11: Annotate resistance from the filtered m0.15 VCF and the 0.15 typing result
    def rules_json = resolved_rules_path

    ch_vcf_for_resistance = FILTER_VCF.out.filtered_vcfs
        .map { run_name, sample_id, vcfs ->
            def filtered_vcf = vcfs.find { it.getName() == "${sample_id}-pilon-m0.15.vcf.gz" }
            if (!filtered_vcf) {
                error "Could not find ${sample_id}-pilon-m0.15.vcf.gz in FILTER_VCF output"
            }
            [sample_id, [run_name, filtered_vcf]]
        }

    ch_gff_for_resistance = ANNOTATE_VADR.out.gff.map { run_name, sample_id, gff ->
        [sample_id, [run_name, gff]]
    }

    ch_iupac_fasta_for_resistance = ch_consensus_with_meta.map { run_name, sample_id, fasta ->
        [sample_id, [run_name, fasta]]
    }

    ch_subtype_for_resistance = ch_iupac_blast_with_meta.map { run_name, sample_id, blast ->
        [sample_id, [run_name, extractSubtypeFromBlast(blast)]]
    }

    ch_resistance_full = ch_vcf_for_resistance
        .join(ch_gff_for_resistance)
        .join(ch_iupac_fasta_for_resistance)
        .join(ch_subtype_for_resistance)
        .map { sample_id, vcf_meta, gff_meta, fasta_meta, subtype_meta ->
            [
                tuple(vcf_meta[0], sample_id, vcf_meta[1], gff_meta[1], fasta_meta[1]),
                subtype_meta[1]
            ]
        }

    ANNOTATE_RESISTANCE(ch_resistance_full.map { it[0] }, ch_resistance_full.map { it[1] }, rules_json)

    // Step 12: Assemble bash-style results contract in one place.
    ch_final_results_input = ch_sample_lids
        .join(ch_hostile_json)
        .join(BAM2FASTA_PILON.out.replacement_fasta.map { run_name, sample_id, fasta, fai -> [sample_id, [fasta, fai]] })
        .join(ch_cram_output.map { run_name, sample_id, cram, crai -> [sample_id, [cram, crai]] })
        .join(LOG_COVERAGE.out.coverage_with_meta.map { run_name, sample_id, coverage -> [sample_id, coverage] })
        .join(BAM2FASTA_BESTREF.out.fasta.map { run_name, sample_id, fasta, fai -> [sample_id, fasta] })
        .join(BAM2FASTA_BESTREF.out.vcf_with_meta.map { run_name, sample_id, vcf -> [sample_id, vcf] })
        .join(BAM2FASTA_BESTREF.out.vcf_index_with_meta.map { run_name, sample_id, vcf_index -> [sample_id, vcf_index] })
        .join(BAM2FASTA_BESTREF.out.stats_with_meta.map { run_name, sample_id, stats -> [sample_id, stats] })
        .join(ch_best_ref_cram_output.map { run_name, sample_id, cram, crai -> [sample_id, [cram, crai]] })
        .join(CREATE_REPORT_BESTREF.out.report_with_meta.map { run_name, sample_id, report -> [sample_id, report] })
        .join(CREATE_REPORT_BESTREF.out.nucfreq_with_meta.map { run_name, sample_id, nucfreq -> [sample_id, nucfreq] })
        .join(ch_consensus_with_meta.map { run_name, sample_id, fasta -> [sample_id, fasta] })
        .join(ch_iupac_cram_output.map { run_name, sample_id, cram, crai -> [sample_id, [cram, crai]] })
        .join(CREATE_REPORT_IUPAC.out.report_with_meta.map { run_name, sample_id, report -> [sample_id, report] })
        .join(CREATE_REPORT_IUPAC.out.nucfreq_with_meta.map { run_name, sample_id, nucfreq -> [sample_id, nucfreq] })
        .join(ANNOTATE_VADR.out.gff.map { run_name, sample_id, gff -> [sample_id, gff] })
        .join(ANNOTATE_VADR.out.bed.map { run_name, sample_id, bed -> [sample_id, bed] })
        .join(FILTER_VCF.out.filtered_vcfs.map { run_name, sample_id, vcfs -> [sample_id, vcfs] })
        .join(FILTER_VCF.out.indices.map { run_name, sample_id, indices -> [sample_id, indices] })
        .join(FILTER_VCF.out.stats.map { run_name, sample_id, stats -> [sample_id, stats] })
        .join(ch_main_blast_with_meta.map { run_name, sample_id, blast -> [sample_id, blast] })
        .join(ch_iupac_blast_with_meta.map { run_name, sample_id, blast -> [sample_id, blast] })
        .join(ch_pilon_iupac_blast_with_meta.map { run_name, sample_id, blast -> [sample_id, blast] })
        .join(ANNOTATE_RESISTANCE.out.tsv_with_meta.map { run_name, sample_id, tsv -> [sample_id, tsv] })
        .join(ANNOTATE_RESISTANCE.out.bed_with_meta.map { run_name, sample_id, bed -> [sample_id, bed] })
        .join(ANNOTATE_RESISTANCE.out.gff_with_meta.map { run_name, sample_id, gff -> [sample_id, gff] })
        .join(ANNOTATE_RESISTANCE.out.drug_tsv_with_meta.map { run_name, sample_id, drug_tsv -> [sample_id, drug_tsv] })
        .map { sample_id, sample_meta, hostile_json_path, main_fasta_meta, main_cram_meta, coverage_tsv,
                bestref_fasta, bestref_vcf, bestref_vcf_index, bestref_vcf_stats, bestref_cram_meta, bestref_report, bestref_nucfreq,
                iupac_fasta, iupac_cram_meta, iupac_report, iupac_nucfreq, vadr_gff, vadr_bed,
                filtered_vcfs, filtered_indices, filtered_stats, main_blast, iupac_blast, pilon_iupac_blast,
                resistance_tsv, resistance_bed, resistance_gff, resistance_drug_tsv ->
            tuple(
                sample_meta[0],
                sample_id,
                sample_meta[1],
                hostile_json_path,
                main_fasta_meta[0],
                main_fasta_meta[1],
                main_blast,
                main_cram_meta[0],
                main_cram_meta[1],
                iupac_fasta,
                iupac_blast,
                iupac_cram_meta[0],
                iupac_cram_meta[1],
                iupac_report,
                iupac_nucfreq,
                bestref_fasta,
                bestref_vcf,
                bestref_vcf_index,
                bestref_vcf_stats,
                bestref_cram_meta[0],
                bestref_cram_meta[1],
                bestref_report,
                bestref_nucfreq,
                coverage_tsv,
                vadr_gff,
                vadr_bed,
                pilon_iupac_blast,
                resistance_tsv,
                resistance_bed,
                resistance_gff,
                resistance_drug_tsv,
                filtered_vcfs,
                filtered_indices,
                filtered_stats
            )
        }

    FINALIZE_RESULTS(ch_final_results_input)

    def sample_info_json_path = resolved_sample_info_json ? resolved_sample_info_json.toString() : ''

    BUILD_QC_SUMMARY(FINALIZE_RESULTS.out.results_dir_with_meta, sample_info_json_path)

    ch_qc_summary_by_run = BUILD_QC_SUMMARY.out.json_with_meta
        .map { run_name, sample_id, qc_json ->
            tuple(run_name, qc_json)
        }
        .groupTuple(by: 0)

    AGGREGATE_QC_SUMMARY(ch_qc_summary_by_run)
    
    // Output final results
    // ch_consensus_with_meta.view { "Final consensus: $it" }
}
