process LOG_COVERAGE {
    tag { "${sample_id}" }
    label 'process_low'
    
    cpus 2
    memory '4 GB'
    time '30m'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(cram), path(crai), path(ref_fasta)
    
    output:
        path "*.tsv", emit: coverage_tsv
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
        
        """
        ${samtools} coverage ${cram} | \\
            awk 'NR==1 {for (i=1;i<=NF;i++) printf "%s\\t", \$i; print ""} NR>1 {for (i=1;i<=NF;i++) printf "%s\\t", \$i; print ""}' | \\
            sed -n '2,\$p' > ${sample_id}-coverage.tsv
        """
    } else {
        """
        echo "LOG_COVERAGE requires container with samtools"
        exit 1
        """
    }
}
