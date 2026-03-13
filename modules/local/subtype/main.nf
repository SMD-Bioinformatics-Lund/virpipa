process SUBTYPE_BLAST {
    tag { sample_id }
    label 'process_medium'
    
    cpus 8
    memory '8 GB'
    time '1h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(fasta)
        val blast_db
    
    output:
        path "*.blast", emit: blast
        path "*.txt", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    
    if (container_dir) {
        def blast = "apptainer exec -B ${bind_paths} ${container_dir}/blast_2.16.0.sif blastn"
        
        """
        # Copy database to local directory to avoid memory map issues
        mkdir -p blast_db
        cp -r ${blast_db}/* blast_db/
        
        ${blast} -query ${fasta} -db blast_db/hcvgluerefs -outfmt 6 > ${sample_id}.blast
        
        # Also create detailed report
        ${blast} -query ${fasta} -db blast_db/hcvgluerefs -out ${sample_id}_blast_detailed.txt -evalue 0.001
        """
    } else {
        """
        blastn -query ${fasta} -db ${blast_db} -outfmt 6 > ${sample_id}.blast
        """
    }
}
