process POLISH_PILON {
    tag { "${sample_id}:${assembly}" }
    label 'process_high'
    
    cpus 16
    memory '32 GB'
    time '4h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/pilon", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(bam), path(bai)
        path assembly
        val assembly_name
    
    output:
        path "*.fasta", emit: polished
        path "*.fai", emit: fai
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def pilon = "apptainer exec -B ${bind_paths} ${container_dir}/pilon-1.24.sif pilon"
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_latest.sif samtools"
        
        """
        ${pilon} --genome ${assembly} --bam ${bam} --output ${sample_id}-${assembly_name} --threads ${task.cpus}
        
        ${samtools} faidx ${sample_id}-${assembly_name}.fasta
        """
    } else {
        def mamba_env = System.getenv('CONDA_PREFIX') ?: '/home/jonas/miniforge3/envs/skrotis'
        
        """
        export PATH="${mamba_env}/bin:\$PATH"
        
        pilon --genome ${assembly} --bam ${bam} --output ${sample_id}-${assembly_name} --threads ${task.cpus}
        
        samtools faidx ${sample_id}-${assembly_name}.fasta
        """
    }
}
