process ASSEMBLE_SPADES {
    tag { sample_id }
    label 'process_high'
    
    cpus 16
    memory '64 GB'
    time '8h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/spades", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2)
    
    output:
        path "${sample_id}.spades", emit: contigs
        path "${sample_id}.spades/*.txt", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def spades = "apptainer exec -B ${bind_paths} ${container_dir}/spades_3.15.5.sif spades.py"
        
        """
        ${spades} --rnaviral -1 ${read1} -2 ${read2} -o ${sample_id}.spades -t ${task.cpus}
        
        mv ${sample_id}.spades/contigs.fasta ${sample_id}.spades 2>/dev/null || true
        """
    } else {
        // Direct execution (local testing with mamba)
        def mamba_env = System.getenv('CONDA_PREFIX') ?: '/home/jonas/miniforge3/envs/skrotis'
        
        """
        export PATH="${mamba_env}/bin:\$PATH"
        
        spades.py --rnaviral -1 ${read1} -2 ${read2} -o ${sample_id}.spades -t ${task.cpus}
        
        mv ${sample_id}.spades/contigs.fasta ${sample_id}.spades 2>/dev/null || true
        """
    }
}
