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
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def hostile_cache = cache_dir && cache_dir.toString().trim() ? "--hostile-cache-dir ${cache_dir}" : ""
    
    if (container_dir) {
        def hostile = "apptainer exec -B ${bind_paths} ${container_dir}/hostile_1.1.0.sif hostile"
        
        """
        ${hostile} -1 ${read1} -2 ${read2} -o . ${hostile_cache} -t ${task.cpus} -p ${sample_id}
        
        mv ${sample_id}_1_hostile-filtered.fastq.gz ${sample_id}_R1_001_hostile.fastq.gz
        mv ${sample_id}_2_hostile-filtered.fastq.gz ${sample_id}_R2_001_hostile.fastq.gz
        
        echo "Hostile filtering complete" > ${sample_id}_hostile.log
        """
    } else {
        """
        hostile -1 ${read1} -2 ${read2} -o . ${hostile_cache} -t ${task.cpus} -p ${sample_id}
        
        mv ${sample_id}_1_hostile-filtered.fastq.gz ${sample_id}_R1_001_hostile.fastq.gz
        mv ${sample_id}_2_hostile-filtered.fastq.gz ${sample_id}_R2_001_hostile.fastq.gz
        
        echo "Hostile filtering complete" > ${sample_id}_hostile.log
        """
    }
}
