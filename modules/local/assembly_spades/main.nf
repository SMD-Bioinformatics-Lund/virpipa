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
        tuple val(run_name), val(sample_id), path("${sample_id}.spades.fasta"), emit: contigs
        path "${sample_id}.spades/*.txt", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    
    if (container_dir) {
        def spades = "${container_runtime} exec -B ${bind_paths} ${container_dir}/spades_3.15.5.sif spades.py"
        
        """
        ${spades} --rnaviral -1 ${read1} -2 ${read2} -o ${sample_id}.spades -t ${task.cpus}
        
        cp ${sample_id}.spades/contigs.fasta ${sample_id}.spades.fasta
        """
    } else {
        // Direct execution (local testing with mamba)
        def mamba_env = System.getenv('CONDA_PREFIX') ?: '/home/jonas/miniforge3/envs/skrotis'
        
        """
        export PATH="${mamba_env}/bin:\$PATH"
        
        spades.py --rnaviral -1 ${read1} -2 ${read2} -o ${sample_id}.spades -t ${task.cpus}
        
        cp ${sample_id}.spades/contigs.fasta ${sample_id}.spades.fasta
        """
    }
}
