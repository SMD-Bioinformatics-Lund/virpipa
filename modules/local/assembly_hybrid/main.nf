process ASSEMBLE_HYBRID {
    tag { "${sample_id}:${genome_name}" }
    label 'process_high'
    
    cpus 16
    memory '32 GB'
    time '4h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/mummer", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(contigs)
        path ref_genome
        val genome_name
    
    output:
        tuple val(run_name), val(sample_id), path("*.hybrid.fasta"), emit: hybrid_assembly
        path "*.delta", emit: delta
        path "*.tiling", emit: tiling
        path "*.log", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def scripts_dir = params.scripts_dir ?: '${projectDir}/scripts'
    
    if (container_dir) {
        def mummer = "apptainer exec -B ${bind_paths} ${container_dir}/mummer3.23.sif"
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
        def bwa = "apptainer exec -B ${bind_paths} ${container_dir}/bwa-0.7.19.sif bwa"
        def python = "apptainer exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python"
        
        """
        # Align contigs to reference with mummer
        ${mummer} nucmer --maxmatch -p ${sample_id} ${ref_genome} ${contigs}
        ${mummer} delta-filter -q ${sample_id}.delta > ${sample_id}.delta-filter
        ${mummer} show-tiling ${sample_id}.delta-filter > ${sample_id}.tiling
        
        # Build hybrid reference
        ${python} ${scripts_dir}/build_hybrid_reference.py \\
            ${ref_genome} ${contigs} ${sample_id}.tiling ${sample_id} \\
            2> ${sample_id}-hybrid-ref.log
        
        # Index hybrid (bwa index optional - container may not exist on HPC)
        ${samtools} faidx ${sample_id}.hybrid.fasta
        ${bwa} index ${sample_id}.hybrid.fasta 2>/dev/null || true
        """
    } else {
        def mamba_env = System.getenv('CONDA_PREFIX') ?: '/home/jonas/miniforge3/envs/skrotis'
        
        """
        export PATH="${mamba_env}/bin:\$PATH"
        
        # Align contigs to reference with mummer
        nucmer --maxmatch -p ${sample_id} ${ref_genome} ${contigs}
        delta-filter -q ${sample_id}.delta > ${sample_id}.delta-filter
        show-tiling ${sample_id}.delta-filter > ${sample_id}.tiling
        
        # Build hybrid reference
        python ${scripts_dir}/build_hybrid_reference.py \\
            ${ref_genome} ${contigs} ${sample_id}.tiling ${sample_id} \\
            2> ${sample_id}-hybrid-ref.log
        
        # Index hybrid
        samtools faidx ${sample_id}.hybrid.fasta
        bwa index ${sample_id}.hybrid.fasta
        """
    }
}
