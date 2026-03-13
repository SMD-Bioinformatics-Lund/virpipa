process SELECT_BEST_REFERENCE {
    tag { "${sample_id}" }
    label 'process_low'
    
    cpus 1
    memory '1 GB'
    time '10m'
    
    input:
        tuple val(run_name), val(sample_id), path(stats_files), path(ref_files)
    
    output:
        tuple val(run_name), val(sample_id), val(best_ref_name), path(best_ref_file), emit: best_ref
        path "*.txt", emit: log
    
    script:
    """
    # Copy ref files to work dir
    cp *.fa . 2>/dev/null || true
    
    # Find the reference with lowest error rate from stats files
    # Stats file pattern: sample-refname.bwa.umi.filter.sort.bam.stats
    
    best_err=""
    best_ref=""
    best_file=""
    
    for stats in *.stats; do
        [[ -e "\$stats" ]] || continue
        
        # Extract ref name - the stats file has pattern: SAMPLE001-1a-AF009606.bwa.umi.filter.sort.bam.stats
        refname=\$(basename "\$stats" | sed 's/.*-\\([^.]*\\).bwa.umi.filter.sort.bam.stats/\\1/')
        
        # Get error rate from stats file
        errrate=\$(awk '\$1=="SN" && \$2=="error" && \$3=="rate:" { print \$4; exit }' "\$stats")
        
        [[ -z "\$errrate" ]] && continue
        
        echo "Found \$refname with error rate \$errrate"
        
        if [[ -z "\$best_err" ]] || (( \$(echo "\$errrate < \$best_err" | bc -l) )); then
            best_err=\$errrate
            best_ref=\$refname
        fi
    done
    
    echo "Best reference: \$best_ref with error rate \$best_err" > best_ref.log
    
    # Find the actual reference file
    for f in *.fa; do
        base=\$(basename "\$f" .fa)
        if [[ "\$base" == "\$best_ref" ]]; then
            best_file="\$f"
            break
        fi
    done
    
    if [[ -z "\$best_file" ]]; then
        echo "ERROR: Could not find ref file for \$best_ref"
        exit 1
    fi
    
    echo "Best ref file: \$best_file"
    echo "\$best_ref" > ref_name.txt
    """
}
