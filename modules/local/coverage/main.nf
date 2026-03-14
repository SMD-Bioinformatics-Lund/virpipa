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
        # Get coverage stats per position using depth
        # Bash pipeline uses region 100-9600 (9501 positions) with fixed total
        ${samtools} depth -r ${sample_id}:100-9600 ${cram} | awk -v sample="${sample_id}" '
        BEGIN { total=9501; cov1=0; cov10=0; cov100=0; cov1000=0 }
        {
            if (\$3 >= 1) cov1++
            if (\$3 >= 10) cov10++
            if (\$3 >= 100) cov100++
            if (\$3 >= 1000) cov1000++
        }
        END {
            pct1 = (cov1/total)*100
            pct10 = (cov10/total)*100
            pct100 = (cov100/total)*100
            pct1000 = (cov1000/total)*100
            print "id\\t1x\\t10x\\t100x\\t1000x"
            printf "%s\\t%.2f\\t%.2f\\t%.2f\\t%.2f\\n", sample, pct1, pct10, pct100, pct1000
        }
        ' > ${sample_id}-coverage.tsv
        """
    } else {
        """
        echo "LOG_COVERAGE requires container with samtools"
        exit 1
        """
    }
}
