process MAP_READS {
    tag { "${sample_id}:${genome_name}" }
    // Use process_medium for local testing, process_high for HPC with sentieon
    
    cpus { params.use_sentieon ? 16 : 8 }
    memory { params.use_sentieon ? '32 GB' : '8 GB' }
    time '8h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bai'
    
    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2), path(genome), val(genome_name)
    
    output:
        tuple val(run_name), val(sample_id), path("*.bam"), path("*.bai"), emit: bams
        path "*.bam.bai", emit: bai
        path "*.stats", emit: stats
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def use_sentieon = params.use_sentieon ?: false
    
    if (use_sentieon && container_dir) {
        // Use sentieon (production)
        def sentieon = "apptainer exec -B ${bind_paths} ${container_dir}/sentieon_202308.03.sif sentieon"
        
        """
        # Use absolute path for genome to ensure bwa index is found
        genome_path=\$(readlink -f ${genome})
        
        # UMI extraction -> bwa mem -> umi consensus
        ${sentieon} umi extract -d 3M2S+T,3M2S+T ${read1} ${read2} | \\
        ${sentieon} bwa mem \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            -k 11 -B 2 -L 25 \\
            -p -C \${genome_path} - | \\
        ${sentieon} umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz
        
        # Align consensus reads
        ${sentieon} bwa mem \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            -k 11 -B 2 -L 25 \\
            -p -C \${genome_path} consensus.fastq.gz | \\
        ${sentieon} util sort -i - --sam2bam --umi_post_process -o ${sample_id}-${genome_name}.bam
        
        # Index
        ${sentieon} util index -a bwa ${sample_id}-${genome_name}.bam
        
        # Stats
        ${sentieon} util stats ${sample_id}-${genome_name}.bam > ${sample_id}-${genome_name}.bam.stats
        """
    } else if (container_dir) {
        // Use bwa + samtools (alternative without sentieon)
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
        def bwa = "apptainer exec -B ${bind_paths} ${container_dir}/bwa-0.7.19.sif bwa"
        def cpus = task.cpus
        def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina"
        
        """
        # Resolve actual genome path (handles both symlinks and regular files)
        actual_path=\$(readlink -f ${genome})
        base=\$(basename \${actual_path})
        
        ${bwa} mem -t ${cpus} -R "${rg}" \${actual_path} ${read1} ${read2} | \\
        ${samtools} view -bS - | \\
        ${samtools} sort -o ${sample_id}-${genome_name}.bam
        
        ${samtools} index ${sample_id}-${genome_name}.bam
        
        ${samtools} stats ${sample_id}-${genome_name}.bam > ${sample_id}-${genome_name}.bam.stats
        """
    } else {
        // Direct execution (local testing)
        """
        echo "MAP_READS requires either sentieon or container. Set params.use_sentieon=true or provide container_dir"
        exit 1
        """
    }
}
