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
        path "*_resistance.tsv", emit: tsv, optional: true
        path "*_resistance.bed", emit: bed, optional: true  
        path "*_resistance_by_drug.tsv", emit: drug_tsv, optional: true
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def scripts_dir = params.scripts_dir ?: '${projectDir}/scripts'
    
    def python = container_dir ?
        "apptainer exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" :
        'python3'
    
    """
    ${python} ${scripts_dir}/annotate_vcf_resistance.py --vcf ${vcf} --gff ${gff} --fasta ${fasta} --subtype ${subtype} --sample-name ${sample_id} --rules ${rules_csv} --output-dir .
    """
}
