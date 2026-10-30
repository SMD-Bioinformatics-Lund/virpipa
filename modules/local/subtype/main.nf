process SUBTYPE_BLAST {
    tag { sample_id }
    label 'process_medium'
    
    cpus 8
    memory '8 GB'
    time '1h'
    
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', pattern: '*.blast', enabled: params.publish_mode == 'debug'

    input:
        tuple val(run_name), val(sample_id), path(fasta)
        val blast_db

    output:
        path "*.blast", emit: blast
        tuple val(run_name), val(sample_id), path("*.blast"), emit: blast_with_meta

    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    def fasta_name = fasta.getName()
    def blast_db_path = blast_db.toString()
    def blast_db_prefix = blast_db_path.endsWith('hcvgluerefs') ? blast_db_path : "${blast_db_path}/hcvgluerefs"
    def blast_db_dir = blast_db_prefix.replaceFirst(/\/hcvgluerefs$/, '')
    def blast_header = 'query acc.ver\tsubject acc.ver\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. star\t\t\tt\ts. end\tevalue\tbit score\tnident\n'
    
    if (container_dir) {
        def blast = "${container_runtime} exec -B ${bind_paths} ${container_dir}/blast_2.16.0.sif blastn"
        def blastdbcmd = "${container_runtime} exec -B ${bind_paths} ${container_dir}/blast_2.16.0.sif blastdbcmd"
        
        """
        set -euo pipefail

        # Copy database to local directory to avoid memory map issues
        mkdir -p blast_db
        cp -r ${blast_db_dir}/* blast_db/

        printf '${blast_header}' > ${fasta_name}.blast
        target_count=\$(${blastdbcmd} -db blast_db/hcvgluerefs -entry all -outfmt '%i' | wc -l)
        ${blast} -query ${fasta} -db blast_db/hcvgluerefs -max_target_seqs "\$target_count" -outfmt '6 std nident' >> ${fasta_name}.blast
        """
    } else {
        """
        set -euo pipefail

        printf '${blast_header}' > ${fasta_name}.blast
        target_count=\$(blastdbcmd -db ${blast_db_prefix} -entry all -outfmt '%i' | wc -l)
        blastn -query ${fasta} -db ${blast_db_prefix} -max_target_seqs "\$target_count" -outfmt '6 std nident' >> ${fasta_name}.blast
        """
    }
}

process RANK_SUBTYPE {
    tag { sample_id }
    label 'process_low'
    cpus 1
    memory '2 GB'
    publishDir { "${params.outdir}/${run_name}/${sample_id}" }, mode: 'copy', enabled: params.publish_mode == 'debug'
    input:
        tuple val(run_name), val(sample_id), path(blast), path(fasta)
        val blast_db
    output:
        tuple val(run_name), val(sample_id), path("${sample_id}_typing.json"), emit: typing
    script:
    def profiles = workflow.profile ?: ''
    def containers = params.container_dir ?: (profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def binds = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (profiles.contains('local') ? '/mnt,/home,/tmp' : '/fs1,/fs2,/local,/mnt/beegfs')
    def runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; else echo singularity; fi)'
    def python = containers ? "${runtime} exec -B ${binds} ${containers}/python_hcvpipe.sif python" : 'python3'
    def scripts = params.scripts_dir != 'scripts' ? params.scripts_dir : "${projectDir}/scripts"
    def prefix = blast_db.toString().endsWith('hcvgluerefs') ? blast_db : "${blast_db}/hcvgluerefs"
    """
    set -euo pipefail
    ${python} ${scripts}/rank_blast.py --blast ${blast} --fasta ${fasta} --database-prefix ${prefix} --output ${sample_id}_typing.json
    """
}
