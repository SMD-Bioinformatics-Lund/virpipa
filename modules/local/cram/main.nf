process CREATE_CRAM {
    tag { "${sample_id}:${bam_base}" }
    label 'process_medium'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(bam), path(bai), path(ref_fasta)
        val(bam_base)
    
    output:
        path "*.cram", emit: crams
        path "*.cram.crai", emit: indices
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
        
        """
        ${samtools} view -O cram,embed_ref -T ${ref_fasta} ${bam} -o ${sample_id}-${bam_base}.cram
        ${samtools} index ${sample_id}-${bam_base}.cram
        """
    } else {
        """
        echo "CREATE_CRAM requires container with samtools"
        exit 1
        """
    }
}
