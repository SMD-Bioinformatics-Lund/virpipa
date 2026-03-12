process REMOVE_HOSTILE {
    tag { sample_id }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/fastq", mode: 'copy', pattern: '*.fastq.gz'
    
    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2)
        val cache_dir
    
    output:
        tuple val(run_name), val(sample_id), path("*_hostile_*.fastq.gz"), emit: reads
        path "*.log", emit: logs
    
    script:
    def r1base = read1.baseName
    def r2base = read2 ? read2.baseName : ''
    def hostile_cache = cache_dir && cache_dir.toString().trim() ? "--hostile-cache-dir ${cache_dir}" : ""
    
    """
    hostile -1 ${read1} -2 ${read2} -o . ${hostile_cache} -t ${task.cpus} -p ${sample_id}
    
    mv ${sample_id}_1_hostile-filtered.fastq.gz ${sample_id}_R1_001_hostile.fastq.gz
    mv ${sample_id}_2_hostile-filtered.fastq.gz ${sample_id}_R2_001_hostile.fastq.gz
    
    echo "Hostile filtering complete" > ${sample_id}_hostile.log
    """
}
