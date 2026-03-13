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
        path "*.tbl", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def scripts_dir = params.scripts_dir ?: '${projectDir}/scripts'
    
    def vadr_container = params.vadr_container ?: "${container_dir}/vadr_164.sif"
    
    if (container_dir) {
        """
        # Run VADR using the script from original pipeline
        bash ${scripts_dir}/vadr_annotate.sh ${fasta} . ${sample_id}
        
        # Move all outputs to current directory
        mv vadr/* . 2>/dev/null || true
        """
    } else {
        """
        echo "VADR requires container. Set container_dir."
        exit 1
        """
    }
}
