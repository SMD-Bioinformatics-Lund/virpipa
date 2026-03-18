process CREATE_REPORT {
    tag { "${sample_id}" }
    label 'process_low'
    
    cpus 2
    memory '4 GB'
    time '10m'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf_stats), path(cram), path(crai), path(ref_fasta)
        val(subtype)
    
    output:
        path "*.report.tsv", emit: report
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
        def ref_name = ref_fasta.baseName
        
        """
        # Create report from VCF stats - matching bash pipeline format
        # Extract ref name from fasta
        ref_name=\$(echo \${ref_fasta} | xargs -I{} basename {} .fasta)
        ref_name=\${ref_name#*-}
        
        echo "# VCF stats" > \${sample_id}.report.tsv
        echo "subtype\t${subtype}" >> \${sample_id}.report.tsv
        echo "reference\t\${ref_name}" >> \${sample_id}.report.tsv
        echo "id\t${sample_id}" >> \${sample_id}.report.tsv
        
        # Extract stats from VCF stats file
        snps=\$(grep "^SN.*number of SNPs" \${vcf_stats} | awk '{print \$NF}')
        multiallelic=\$(grep "^SN.*number of multiallelic SNP sites" \${vcf_stats} | awk '{print \$NF}')
        af0=\$(grep "^AF.*0.000000" \${vcf_stats} | awk '{print \$4}')
        af99=\$(grep "^AF.*0.990000" \${vcf_stats} | awk '{print \$4}')
        
        echo "snps\t\${snps:-0}" >> \${sample_id}.report.tsv
        echo "multiallelic_snps\t\${multiallelic:-0}" >> \${sample_id}.report.tsv
        echo "AF0\t\${af0:-0}" >> \${sample_id}.report.tsv
        echo "AF0.99\t\${af99:-0}" >> \${sample_id}.report.tsv
        
        # Add COVERAGE section - matching bash pipeline format
        # Bash uses: samtools coverage | datamash transpose | sed -n '2,\$p'
        echo "# COVERAGE" >> \${sample_id}.report.tsv
        ${samtools} coverage \${cram} | \\
            awk '{
                for (i=1; i<=NF; i++) {
                    rows[i] = (rows[i] ? rows[i] "\\t" : "") \$i
                }
            }
            END {
                for (i=1; i<=NF; i++) {
                    print rows[i]
                }
            }' | \\
            sed -n '2,\$p' >> \${sample_id}.report.tsv
        """
    } else {
        error "CREATE_REPORT requires container_dir"
    }
}
