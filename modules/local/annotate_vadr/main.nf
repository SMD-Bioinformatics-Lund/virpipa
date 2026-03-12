process ANNOTATE_VADR {
    tag { sample_id }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(fasta)
        val vadr_model
    
    output:
        path "*_mod.gff", emit: gff
        path "*.bed", emit: bed
        path "*.txt", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def scripts_dir = params.scripts_dir ?: '${projectDir}/scripts'
    def vadr_model = params.vadr_model ?: 'vadr-models-flavi'
    def vadr_container = params.vadr_container ?: "${container_dir}/vadr_164.sif"
    
    if (container_dir) {
        def vadr = "apptainer exec -B ${bind_paths} ${vadr_container} vadr"
        
        """
        ${vadr} -i ${fasta} -o vadr_out --ncpu ${task.cpus} --fvadr_models ${vadr_model}
        
        mv vadr_out/*_mod.gff ${sample_id}_vadr_mod.gff 2>/dev/null || true
        mv vadr_out/*.bed ${sample_id}_vadr.bed 2>/dev/null || true
        """
    } else {
        """
        echo "VADR requires container. Set container_dir."
        exit 1
        """
    }
}
