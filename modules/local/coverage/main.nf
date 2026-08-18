process LOG_COVERAGE {
    tag { "${sample_id}" }
    label 'process_low'
    
    cpus 2
    memory '4 GB'
    time '30m'
    
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', enabled: params.publish_mode == 'debug'
    
    input:
        tuple val(run_name), val(sample_id), path(cram), path(crai), path(ref_fasta)
    
    output:
        path "${sample_id}-coverage.tsv", emit: coverage_tsv
        path "${sample_id}-coverage-1x.bed", emit: coverage_1x_bed
        tuple val(run_name), val(sample_id), path("${sample_id}-coverage.tsv"), path("${sample_id}-coverage-1x.bed"), emit: coverage_with_meta
    
    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'

    def samtools = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"

    """
    set -euo pipefail

    printf "id\\t1x\\t10x\\t100x\\t1000x\\n${sample_id}\\t" > ${sample_id}-coverage.tsv
    ${samtools} depth -r "${sample_id}:100-9600" ${cram} | \\
        awk 'BEGIN { total=9501; cov1=0; cov10=0; cov100=0; cov1000=0 }
        { if (\$3 >= 1) cov1++; if (\$3 >= 10) cov10++; if (\$3 >= 100) cov100++; if (\$3 >= 1000) cov1000++ }
        END {
            printf "%.2f\\t", (cov1/total)*100;
            printf "%.2f\\t", (cov10/total)*100;
            printf "%.2f\\t", (cov100/total)*100;
            printf "%.2f\\n", (cov1000/total)*100;
        }' >> ${sample_id}-coverage.tsv

    ${samtools} depth -aa --reference ${ref_fasta} ${cram} | \
        awk 'BEGIN { OFS="\\t" }
        \$3 >= 1 {
            start = \$2 - 1
            if (active && \$1 == chrom && start == end) {
                end = \$2
            } else {
                if (active) print chrom, interval_start, end
                chrom = \$1
                interval_start = start
                end = \$2
                active = 1
            }
        }
        \$3 < 1 && active {
            print chrom, interval_start, end
            active = 0
        }
        END { if (active) print chrom, interval_start, end }' > ${sample_id}-coverage-1x.bed
    """
}
