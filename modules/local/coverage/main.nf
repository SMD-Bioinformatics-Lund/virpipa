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
        # Bash pipeline uses region 100-9600 (9501 positions)
        ${samtools} depth -r ${sample_id}:100-9600 ${cram} | awk -v sample="${sample_id}" '
        {
            if (\$3 >= 1) c1++
            if (\$3 >= 10) c10++
            if (\$3 >= 100) c100++
            if (\$3 >= 1000) c1000++
            total++
        }
        END {
            if (total > 0) {
                pct1 = (c1/total)*100
                pct10 = (c10/total)*100
                pct100 = (c100/total)*100
                pct1000 = (c1000/total)*100
            } else {
                pct1 = pct10 = pct100 = pct1000 = 0
            }
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
