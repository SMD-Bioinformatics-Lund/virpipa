process ANNOTATE_RESISTANCE {
    tag { sample_id }
    label 'process_medium'
    
    cpus 4
    memory '8 GB'
    time '1h'
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', pattern: '*_resistance*', enabled: params.publish_mode == 'debug'
    
    input:
        tuple val(run_name), val(sample_id), path(gff), path(fasta), path(cram), path(crai), path(track_fasta), path(consensus_metadata), path(track_metadata)
        val subtype
        path rules_json
        path h77_fasta
        path reference_panel
    
    output:
        path "*_resistance.tsv", emit: tsv
        path "*_resistance.bed", emit: bed
        path "*_resistance.gff", emit: gff
        path "*_resistance_by_drug.tsv", emit: drug_tsv
        path "*_resistance.json", emit: json
        path "*_resistance_sites.gff3", emit: sites_gff
        tuple val(run_name), val(sample_id), path("*_resistance.tsv"), emit: tsv_with_meta
        tuple val(run_name), val(sample_id), path("*_resistance.bed"), emit: bed_with_meta
        tuple val(run_name), val(sample_id), path("*_resistance.gff"), emit: gff_with_meta
        tuple val(run_name), val(sample_id), path("*_resistance_by_drug.tsv"), emit: drug_tsv_with_meta
        tuple val(run_name), val(sample_id), path("*_resistance.json"), emit: json_with_meta
        tuple val(run_name), val(sample_id), path("*_resistance_sites.gff3"), emit: sites_gff_with_meta
    
    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    def scripts_dir = params.scripts_dir != 'scripts' ? params.scripts_dir : "${projectDir}/scripts"
    
    def python = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" :
        'python3'
    
    def track_arg = track_fasta ? "--track-reference-fasta ${track_fasta}" : ''
    def cram_arg = cram ? "--cram ${cram}" : ''
    def metadata_arg = consensus_metadata ? "--consensus-metadata ${consensus_metadata}" : '--allow-legacy-consensus'
    def track_metadata_arg = track_metadata ? "--track-reference-metadata ${track_metadata}" : ''
    """
    ${python} ${scripts_dir}/annotate_vcf_resistance.py --gff ${gff} --fasta ${fasta} ${cram_arg} --subtype ${subtype} --sample-name ${sample_id} --rules ${rules_json} --h77-fasta ${h77_fasta} --reference-panel ${reference_panel}/panel.json ${track_arg} ${metadata_arg} ${track_metadata_arg} --output-dir .
    """
}
