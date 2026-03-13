process CREATE_CONSENSUS {
    tag { "${sample_id}:${vcf}" }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf), path(fasta), path(fai)
        val min_freq
    
    output:
        path "*iupac.fasta", emit: consensus
        path "*.fasta", emit: fasta
        path "*.fai", emit: fai
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def scripts_dir = params.scripts_dir ?: '${projectDir}/scripts'
    
    def python_cmd = container_dir ? 
        "apptainer exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" :
        "python3"
    
    """
    ${python_cmd} ${scripts_dir}/consensus_fasta_iupac.py ${fasta} ${vcf} ${min_freq}
    """
}
