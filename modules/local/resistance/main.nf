process ANNOTATE_RESISTANCE {
    tag { sample_id }
    label 'process_medium'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf), path(gff), path(fasta)
        val subtype
        path rules_csv
    
    output:
        path "*_resistance.tsv", emit: tsv
        path "*_resistance.bed", emit: bed
        path "*_resistance_by_drug.tsv", emit: drug_tsv
    
    script:
    def scripts_dir = params.scripts_dir ?: '${projectDir}/scripts'
    
    """
    python ${scripts_dir}/annotate_vcf_resistance.py \\
        --vcf ${vcf} \\
        --gff ${gff} \\
        --fasta ${fasta} \\
        --subtype ${subtype} \\
        --sample-name ${sample_id} \\
        --rules ${rules_csv} \\
        --output-dir .
    """
}
