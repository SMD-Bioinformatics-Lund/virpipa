process MAP_READS {
    tag { "${sample_id}:${genome}" }
    label 'process_high'
    
    cpus 16
    memory '32 GB'
    time '8h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bai'
    
    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2)
        path genome
        val genome_name
    
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
        # UMI extraction and consensus
        ${sentieon} umi extract -d 3M2S+T,3M2S+T ${read1} ${read2} -o umi_r1.fq.gz -O umi_r2.fq.gz
        
        # BWA alignment with UMI
        ${sentieon} bwa mem \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            -k 11 -B 2 -L 25 \\
            -p -C ${genome} umi_r1.fq.gz umi_r2.fq.gz | \\
        ${sentieon} util sort -i - --sam2bam -o ${sample_id}-${genome_name}.bam
        
        # Index
        ${sentieon} util index -a bwa ${sample_id}-${genome_name}.bam
        
        # Stats
        ${sentieon} util stats ${sample_id}-${genome_name}.bam > ${sample_id}-${genome_name}.bam.stats
        """
    } else if (container_dir) {
        // Use fgbio + bwa (alternative without sentieon)
        def fgbio = "apptainer exec -B ${bind_paths} ${container_dir}/fgbio_latest.sif fgbio"
        def bwa = "apptainer exec -B ${bind_paths} ${container_dir}/bwa_latest.sif bwa"
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_latest.sif samtools"
        
        """
        # Extract UMI from read names (fgbio)
        ${fgbio} ExtractUmisFromBam -i ${read1} -o unmapped.bam --read-name-in
		
        # Group by UMI
        ${fgbio} GroupReadsByUmi -i unmapped.bam -o grouped.bam --strategy directional
		
        # Call consensus
        ${fgbio} CallConsensus -i grouped.bam -o consensus.bam --min-reads 3
		
        # Align consensus
        ${bwa} mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" ${genome} consensus.bam | \\
        ${samtools} view -bS - | \\
        ${samtools} sort -o ${sample_id}-${genome_name}.bam
		
        # Index
        ${samtools} index ${sample_id}-${genome_name}.bam
		
        # Stats
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
