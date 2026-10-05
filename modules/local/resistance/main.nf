process ANNOTATE_RESISTANCE {
    tag { sample_id }
    label 'process_medium'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    input:
        tuple val(run_name), val(sample_id), path(gff), path(fasta), path(cram), path(crai)
        val subtype
        path rules_json
        path h77_fasta
    
    output:
        path "*_resistance.tsv", emit: tsv, optional: true
        path "*_resistance.bed", emit: bed, optional: true
        path "*_resistance.gff", emit: gff, optional: true
        path "*_resistance_by_drug.tsv", emit: drug_tsv, optional: true
        path "*_resistance.json", emit: json, optional: true
        path "*_resistance_sites.gff3", emit: sites_gff, optional: true
        tuple val(run_name), val(sample_id), path("*_resistance.tsv"), emit: tsv_with_meta, optional: true
        tuple val(run_name), val(sample_id), path("*_resistance.bed"), emit: bed_with_meta, optional: true
        tuple val(run_name), val(sample_id), path("*_resistance.gff"), emit: gff_with_meta, optional: true
        tuple val(run_name), val(sample_id), path("*_resistance_by_drug.tsv"), emit: drug_tsv_with_meta, optional: true
        tuple val(run_name), val(sample_id), path("*_resistance.json"), emit: json_with_meta, optional: true
        tuple val(run_name), val(sample_id), path("*_resistance_sites.gff3"), emit: sites_gff_with_meta, optional: true
    
    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    def scripts_dir = params.scripts_dir != 'scripts' ? params.scripts_dir : "${projectDir}/scripts"
    
    def python = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" :
        'python3'
    
    """
    ${python} ${scripts_dir}/annotate_vcf_resistance.py --gff ${gff} --fasta ${fasta} --cram ${cram} --subtype ${subtype} --sample-name ${sample_id} --rules ${rules_json} --h77-fasta ${h77_fasta} --output-dir .
    """
}
