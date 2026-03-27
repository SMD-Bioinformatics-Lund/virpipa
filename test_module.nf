#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.module = params.module ?: 'help'
params.subsample_reads = params.subsample_reads ?: 1000

include { REMOVE_HOSTILE } from './modules/local/hostile/main'
include { SUBSAMPLE_READS } from './modules/local/subsample/main'
include { BAM2FASTA } from './modules/local/bam2fasta/main'
include { CREATE_CONSENSUS } from './modules/local/consensus/main'
include { FILTER_VCF } from './modules/local/filter_vcf/main'
include { VARIANT_CALLING } from './modules/local/variantcall/main'
include { CREATE_CRAM } from './modules/local/cram/main'
include { LOG_COVERAGE } from './modules/local/coverage/main'
include { SUBTYPE_BLAST } from './modules/local/subtype/main'
include { CREATE_REPORT } from './modules/local/report/main'

workflow {
    def test_input = params.input ?: "${projectDir}/assets/test_samplesheet.csv"
    def test_outdir = params.outdir

    if (params.module == 'help') {
        println """
============================================
VirPipa Module Test Workflow
============================================

Available modules to test locally:
  - hostile    : Remove human reads with hostile
  - subsample  : Subsample FASTQ reads with seqtk
  - bam2fasta  : Build consensus FASTA and VCF from a fixture BAM
  - consensus  : Build 0.15 IUPAC consensus from a fixture VCF
  - filter_vcf : Build min-fraction filtered VCFs from a fixture pilon VCF
  - variantcall: Build the pilon VCF from a fixture BAM and reference
  - cram       : Build the polished CRAM from a fixture BAM and reference
  - coverage   : Build the coverage TSV from a fixture CRAM
  - subtype    : Build the BLAST subtype hits from a fixture consensus FASTA
  - report     : Build the report TSV and nucleotide frequencies from report fixtures

Usage:
  nextflow run test_module.nf -profile local_containers,tiny --module hostile
  nextflow run test_module.nf -profile local,tiny --module subsample --subsample_reads 100
  nextflow run test_module.nf -profile local_containers,tiny --module hostile --outdir test_output_tiny
  nextflow run test_module.nf -profile local --module bam2fasta --outdir test_output_bam2fasta
  nextflow run test_module.nf -profile local --module consensus --outdir test_output_consensus
  nextflow run test_module.nf -profile local_containers --module filter_vcf --outdir test_output_filter_vcf
  nextflow run test_module.nf -profile local_containers --module variantcall --outdir test_output_variantcall
  nextflow run test_module.nf -profile local_containers --module cram --outdir test_output_cram
  nextflow run test_module.nf -profile local --module coverage --outdir test_output_coverage
  nextflow run test_module.nf -profile local_containers --module subtype --outdir test_output_subtype
  nextflow run test_module.nf -profile local --module report --outdir test_output_report

Notes:
  - Sentieon-dependent modules are intentionally not wired here yet.
  - Default input is ${test_input}
  - Default output follows `params.outdir` from config unless `--outdir` is provided.
  - `local` uses tools from `skrotis`; `local_containers` uses Apptainer images.

============================================
"""
        return
    }

    def resolvePath = { raw_path ->
        def input_path = raw_path.toString().trim()
        if (!input_path) {
            return null
        }

        def candidate = file(input_path)
        if (candidate.exists()) {
            return candidate
        }

        if (input_path.startsWith('/fs1/')) {
            def mounted_path = file("/mnt${input_path}")
            if (mounted_path.exists()) {
                return mounted_path
            }
        }

        return candidate
    }

    Channel
        .fromPath(test_input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = (row.clarity_sample_id ?: row.sample ?: row.id ?: '').toString().trim()
            if (!sample_id) {
                error "Samplesheet row is missing sample id"
            }

            def read1 = (row.read1 ?: row.fastq_1 ?: '').toString().trim()
            if (!read1) {
                error "Samplesheet row for sample '${sample_id}' is missing read1/fastq_1"
            }

            def read2 = (row.read2 ?: row.fastq_2 ?: '').toString().trim()
            if (!read2) {
                error "Samplesheet row for sample '${sample_id}' is missing read2/fastq_2"
            }

            def run_name = (row.run_name ?: params.run_name ?: 'test').toString().trim()
            tuple(run_name, sample_id, resolvePath(read1), resolvePath(read2))
        }
        .set { ch_samples }

    if (params.module == 'hostile') {
        REMOVE_HOSTILE(ch_samples, params.hostile_cache_dir ?: '')
    } else if (params.module == 'subsample') {
        SUBSAMPLE_READS(ch_samples, params.subsample_reads)
    } else if (params.module == 'bam2fasta') {
        BAM2FASTA(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/bam2fasta/SAMPLE001-pilon.r11b2L25.bwa.umi.filter.sort.bam"),
                    file("${projectDir}/assets/test_data/bam2fasta/SAMPLE001-pilon.r11b2L25.bwa.umi.filter.sort.bam.bai"),
                    file("${projectDir}/assets/test_data/bam2fasta/SAMPLE001.fasta.org"),
                    '1.0-iupac'
                )
            ),
            1.0
        )
    } else if (params.module == 'consensus') {
        CREATE_CONSENSUS(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/consensus/SAMPLE001-pilon.vcf.gz"),
                    file("${projectDir}/assets/test_data/consensus/SAMPLE001.fasta")
                )
            ),
            '0.15'
        )
    } else if (params.module == 'filter_vcf') {
        FILTER_VCF(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/filter_vcf/SAMPLE001-pilon.vcf.gz"),
                    file("${projectDir}/assets/test_data/filter_vcf/SAMPLE001-pilon.vcf.gz.csi"),
                    'SAMPLE001-pilon'
                )
            )
        )
    } else if (params.module == 'variantcall') {
        VARIANT_CALLING(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/variantcall/SAMPLE001-pilon.r11b2L25.bwa.umi.filter.sort.bam"),
                    file("${projectDir}/assets/test_data/variantcall/SAMPLE001-pilon.r11b2L25.bwa.umi.filter.sort.bam.bai"),
                    file("${projectDir}/assets/test_data/variantcall/SAMPLE001.fasta"),
                    'pilon'
                )
            )
        )
    } else if (params.module == 'cram') {
        CREATE_CRAM(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/cram/SAMPLE001-pilon.r11b2L25.bwa.umi.filter.sort.bam"),
                    file("${projectDir}/assets/test_data/cram/SAMPLE001-pilon.r11b2L25.bwa.umi.filter.sort.bam.bai"),
                    file("${projectDir}/assets/test_data/cram/SAMPLE001.fasta"),
                    'SAMPLE001'
                )
            )
        )
    } else if (params.module == 'coverage') {
        LOG_COVERAGE(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/coverage/SAMPLE001.cram"),
                    file("${projectDir}/assets/test_data/coverage/SAMPLE001.cram.crai"),
                    file("${projectDir}/assets/test_data/coverage/SAMPLE001.fasta")
                )
            )
        )
    } else if (params.module == 'subtype') {
        SUBTYPE_BLAST(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/subtype/SAMPLE001-0.15-iupac.fasta")
                )
            ),
            Channel.value(file('/mnt/fs1/jonas/hcv/refgenomes/hcvglue'))
        )
    } else if (params.module == 'report') {
        CREATE_REPORT(
            Channel.of(
                tuple(
                    'fixture_run',
                    'SAMPLE001',
                    file("${projectDir}/assets/test_data/report/SAMPLE001-0.15-iupac.vcf.gz.stats"),
                    file("${projectDir}/assets/test_data/report/SAMPLE001-0.15-iupac.cram"),
                    file("${projectDir}/assets/test_data/report/SAMPLE001-0.15-iupac.cram.crai"),
                    file("${projectDir}/assets/test_data/report/SAMPLE001-0.15-iupac.fasta"),
                    '3a-D17763'
                )
            )
        )
    } else {
        error "Unsupported module '${params.module}'. Supported modules: hostile, subsample, bam2fasta, consensus, filter_vcf, variantcall, cram, coverage, subtype, report"
    }
}
