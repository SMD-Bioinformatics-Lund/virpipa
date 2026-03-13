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
        awk 'NR==1 {for(i=1;i<=NF;i++) h[i]=\$i; next} {for(i=1;i<=NF;i++) d[i]=d[i] (d[i]?FS:"") \$i} END {for(i=1;i<=NF;i++) print h[i] FS d[i]}' coverage_data.tsv > \${sample_id}.coverage.tsv
        
        # Calculate coverage at different thresholds as percentages
        # Columns in samtools coverage: rname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq
        awk '
        NR>1 {
            covbases=\$5
            total=\$3
            if (covbases >= 1) c1+=covbases
            if (covbases >= 10) c10+=covbases
            if (covbases >= 100) c100+=covbases
            if (covbases >= 1000) c1000+=covbases
            total_len+=\$3
        }
        END {
            if (total_len > 0) {
                printf "id\t1x\t10x\t100x\t1000x\n"
                printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\n", "\${sample_id}", (c1/total_len)*100, (c10/total_len)*100, (c100/total_len)*100, (c1000/total_len)*100
            }
        }' coverage_data.tsv > \${sample_id}-coverage.tsv
        
        rm coverage_header.tsv coverage_data.tsv
        """
    } else {
        """
        echo "LOG_COVERAGE requires container with samtools"
        exit 1
        """
    }
}
