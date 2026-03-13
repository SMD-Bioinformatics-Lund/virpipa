process SUBSAMPLE_READS {
    tag { sample_id }
    label 'process_low'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/fastq", mode: 'copy', pattern: '*.fastq.gz'
    
    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2)
        val(nreads)
    
    output:
        tuple val(run_name), val(sample_id), path("*_R1_*.fastq.gz"), path("*_R2_*.fastq.gz"), emit: reads
        path "*.log", emit: logs
    
    script:
    def r1base = read1.baseName
    def r2base = read2 ? read2.baseName : ''
    def r2_arg = read2 ? "${read2}" : ""
    def r2_out = read2 ? "${sample_id}_R2_001.fastq.gz" : ""
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    // Use apptainer if container_dir is set
    // Otherwise use direct commands (assumes tools in PATH)
    if (container_dir) {
        def seqtk = "apptainer exec -B ${bind_paths} ${container_dir}/seqtk_1.3.sif seqtk"
        def pigz = "apptainer exec -B ${bind_paths} ${container_dir}/pigz-2.3.4.sif pigz"
        
        """
        echo "Subsampling to ${nreads} reads" > subsample.log
        echo "Input R1: ${read1}" >> subsample.log
        echo "Input R2: ${read2}" >> subsample.log
        
        ${seqtk} sample -s100 ${read1} ${nreads} | ${pigz} -p 4 -c > ${sample_id}_R1_001.fastq.gz
        ${r2_arg ? "${seqtk} sample -s100 ${read2} ${nreads} | ${pigz} -p 4 -c > ${sample_id}_R2_001.fastq.gz" : "touch ${sample_id}_R2_001.fastq.gz"}
        
        echo "Output R1: ${sample_id}_R1_001.fastq.gz" >> subsample.log
        echo "Output R2: ${sample_id}_R2_001.fastq.gz" >> subsample.log
        """
    } else {
        // Direct execution - use full path to mamba environment
        def mamba_env = System.getenv('CONDA_PREFIX') ?: '/home/jonas/miniforge3/envs/skrotis'
        """
        export PATH="${mamba_env}/bin:\$PATH"
        
        echo "Subsampling to ${nreads} reads" > subsample.log
        echo "Input R1: ${read1}" >> subsample.log
        echo "Input R2: ${read2}" >> subsample.log
        
        seqtk sample -s100 ${read1} ${nreads} | pigz -p 4 -c > ${sample_id}_R1_001.fastq.gz
        ${r2_arg ? "seqtk sample -s100 ${read2} ${nreads} | pigz -p 4 -c > ${sample_id}_R2_001.fastq.gz" : "touch ${sample_id}_R2_001.fastq.gz"}
        
        echo "Output R1: ${sample_id}_R1_001.fastq.gz" >> subsample.log
        echo "Output R2: ${sample_id}_R2_001.fastq.gz" >> subsample.log
        """
    }
}
