process SELECT_BEST_REFERENCE {
    tag { "${sample_id}" }
    label 'process_low'
    
    cpus 1
    memory '1 GB'
    time '10m'
    
    input:
        tuple val(run_name), val(sample_id), path(stats_files), val(ref_dir)
    
    output:
        tuple val(run_name), val(sample_id), val(best_ref_name), path("*.fa"), emit: best_ref
        path "*.txt", emit: log
    
    script:
    def refDir = ref_dir
    """
    # Get absolute path to ref_dir
    REF_DIR=\$(readlink -f "${refDir}")
    
    # Find the reference with lowest error rate from stats files
    
    best_err=""
    best_ref=""
    best_file=""
    
    for stats in *.stats; do
        [[ -e "\$stats" ]] || continue
        
        # Extract ref name - the stats file has pattern: SAMPLE001-3a-D17763.bam.stats
        # We want to extract "3a-D17763" - capture everything between first - and .bam.stats
        refname=\$(echo "\$stats" | sed 's/^[^-]*-\\([^.]*\\).*/\\1/')
        
        # Get error rate from stats file (convert scientific notation to decimal for bc)
        errrate=\$(awk '\$1=="SN" && \$2=="error" && \$3=="rate:" { printf "%.6f", \$4; exit }' "\$stats")
        
        [[ -z "\$errrate" ]] && continue
        
        echo "Found \$refname with error rate \$errrate"
        
        if [[ -z "\$best_err" ]] || (( \$(echo "\$errrate < \$best_err" | bc -l) )); then
            best_err=\$errrate
            best_ref=\$refname
        fi
    done
    
    echo "Best reference: \$best_ref with error rate \$best_err" > best_ref.log
    
    # Find the reference file in ref_dir
    best_file="\${REF_DIR}/\${best_ref}.fa"
    
    if [[ ! -f "\$best_file" ]]; then
        echo "ERROR: Could not find ref file: \$best_file"
        ls -la "\${REF_DIR}/"*.fa | head -5
        exit 1
    fi
    
    echo "Best ref file: \$best_file"
    cp "\$best_file" .
    
    # Write best ref name to file for Nextflow output
    echo "\$best_ref" > best_ref_name.txt
    """
}
