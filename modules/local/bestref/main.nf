process SELECT_BEST_REFERENCE {
    tag { "${sample_id}" }
    label 'process_low'
    
    cpus 1
    memory '1 GB'
    time '10m'
    
    input:
        tuple val(run_name), val(sample_id), path(stats_files), val(ref_dir)
    
    output:
        tuple val(run_name), val(sample_id), path("*.fa"), emit: best_ref
        path "*.txt", emit: log
    
    script:
    def refDir = ref_dir
    """
    set -euo pipefail

    # Get absolute path to ref_dir
    REF_DIR=\$(readlink -f "${refDir}")
    
    # Find the reference with lowest error rate from stats files
    
    best_ref=""
    best_file=""
    > best_ref_candidates.tsv
    
    for stats in *.stats; do
        [[ -e "\$stats" ]] || continue
        
        # Extract the reference name from bash-style stats names like:
        # SAMPLE001-3a-D17763.r11b2L25.bwa.umi.filter.sort.bam.stats
        refname=\$(basename "\$stats")
        refname=\${refname#${sample_id}-}
        refname=\${refname%.bwa.umi.filter.sort.bam.stats}
        refname=\${refname%.*}
        
        # Get error rate from stats file (convert scientific notation to decimal for bc)
        errrate=\$(awk '\$1=="SN" && \$2=="error" && \$3=="rate:" { printf "%.6f", \$4; exit }' "\$stats")
        
        [[ -z "\$errrate" ]] && continue
        
        printf "%s\t%s\n" "\$refname" "\$errrate" | tee -a best_ref_candidates.tsv
    done

    best_ref=\$(sort -k2,2g best_ref_candidates.tsv | head -n 1 | cut -f1)
    best_err=\$(sort -k2,2g best_ref_candidates.tsv | head -n 1 | cut -f2)
    
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
