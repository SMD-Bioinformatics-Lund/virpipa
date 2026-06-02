process ANNOTATE_VADR {
    tag { sample_id }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', enabled: params.publish_mode == 'debug'
    
    input:
        tuple val(run_name), val(sample_id), path(fasta)
        val vadr_model
    
    output:
        tuple val(run_name), val(sample_id), path("*_mod.gff"), emit: gff
        tuple val(run_name), val(sample_id), path("*.bed"), emit: bed

    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    def scripts_dir = params.scripts_dir != 'scripts' ? params.scripts_dir : (active_profiles.contains('hpc') ? '/fs1/jonas/src/virpipa/scripts' : "${projectDir}/scripts")
    def vadr_container = params.vadr_container ?: "${container_dir}/vadr_164.sif"
    def configured_model_dir = params.vadr_model_dir ?: (active_profiles.contains('local') ? '/home/jonas/resources/vadr/vadr-models-flavi' : '/fs1/resources/ref/micro/vadr/vadr-models-flavi')
    def vadr_model_dir = vadr_model.toString().contains('/') ? vadr_model.toString() : (configured_model_dir ?: vadr_model.toString())
    
    if (container_dir) {
        """
        set -euo pipefail

        export VADR_CONTAINER='${vadr_container}'
        export VADR_BIND='${bind_paths}'
        export VADR_CONTAINER_RUNTIME=${container_runtime}
        export VADR_MODELDIR='${vadr_model_dir}'
        export VADR_ANNOTATE_TBL2GFF='${scripts_dir}/annotate-tbl2gff.pl'

        # Run VADR using the original pipeline helper with explicit local paths.
        bash ${scripts_dir}/vadr_annotate.sh ${fasta} . ${sample_id}

        shopt -s nullglob
        for output in vadr/${sample_id}*_mod.gff vadr/${sample_id}*.bed; do
            mv "\$output" .
        done
        shopt -u nullglob
        """
    } else {
        """
        echo "VADR requires container. Set container_dir."
        exit 1
        """
    }
}
