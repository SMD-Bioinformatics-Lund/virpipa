process SUBTYPE_BLAST {
    tag { sample_id }
    label 'process_medium'
    
    cpus 8
    memory '8 GB'
    time '1h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy', pattern: '*.blast'

    input:
        tuple val(run_name), val(sample_id), path(fasta)
        val blast_db

    output:
        path "*.blast", emit: blast
        tuple val(run_name), val(sample_id), path("*.blast"), emit: blast_with_meta

    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def fasta_name = fasta.getName()
    def blast_db_path = blast_db.toString()
    def blast_db_prefix = blast_db_path.endsWith('hcvgluerefs') ? blast_db_path : "${blast_db_path}/hcvgluerefs"
    def blast_db_dir = blast_db_prefix.replaceFirst(/\/hcvgluerefs$/, '')
    def blast_header = 'query acc.ver\tsubject acc.ver\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. star\t\t\tt\ts. end\tevalue\tbit score\n'
    
    if (container_dir) {
        def blast = "apptainer exec -B ${bind_paths} ${container_dir}/blast_2.16.0.sif blastn"
        
        """
        set -euo pipefail

        # Copy database to local directory to avoid memory map issues
        mkdir -p blast_db
        cp -r ${blast_db_dir}/* blast_db/

        printf '${blast_header}' > ${fasta_name}.blast
        ${blast} -query ${fasta} -db blast_db/hcvgluerefs -outfmt 6 >> ${fasta_name}.blast
        """
    } else {
        """
        set -euo pipefail

        printf '${blast_header}' > ${fasta_name}.blast
        blastn -query ${fasta} -db ${blast_db_prefix} -outfmt 6 >> ${fasta_name}.blast
        """
    }
}
