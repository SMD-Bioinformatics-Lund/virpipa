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
        # Get coverage stats
        ${samtools} coverage -H ${cram} > coverage_header.tsv
        ${samtools} coverage ${cram} > coverage_data.tsv
        
        # Transpose and format
        awk 'NR==1 {for(i=1;i<=NF;i++) h[i]=\$i; next} {for(i=1;i<=NF;i++) d[i]=d[i] (d[i]?FS:"") \$i} END {for(i=1;i<=NF;i++) print h[i] FS d[i]}' coverage_data.tsv > ${sample_id}.coverage.tsv
        
        # Calculate coverage at different thresholds
        # Columns: coverage, bases, fraction
        awk '
        NR>1 {
            cov=\$3
            if (cov >= 1) c1++
            if (cov >= 10) c10++
            if (cov >= 100) c100++
            if (cov >= 1000) c1000++
            total++
        }
        END {
            if (total > 0) {
                print "1x\t" c1/total
                print "10x\t" c10/total  
                print "100x\t" c100/total
                print "1000x\t" c1000/total
            }
        }' coverage_data.tsv > ${sample_id}-coverage.tsv
        
        rm coverage_header.tsv coverage_data.tsv
        """
    } else {
        """
        echo "LOG_COVERAGE requires container with samtools"
        exit 1
        """
    }
}
