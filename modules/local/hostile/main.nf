process REMOVE_HOSTILE {
    tag { sample_id }
    label 'process_medium'
    
    cpus 8
    memory '16 GB'
    time '2h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy', pattern: 'hostile.json'
    publishDir "${params.outdir}/${run_name}/${sample_id}/fastq", mode: 'copy', pattern: '*.fastq.gz'
    
    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2)
        val cache_dir
    
    output:
        tuple val(run_name), val(sample_id), path("*_R1_*_hostile*.fastq.gz"), path("*_R2_*_hostile*.fastq.gz"), emit: reads
        path "hostile.json", emit: hostile_json
        path "*.log", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def hostile_cache = params.hostile_cache_dir ?: '/fs1/resources/ref/micro/hostile'
    
    if (container_dir) {
        def hostile = "apptainer exec -B ${bind_paths} ${container_dir}/hostile_1.1.0.sif hostile"
        
        """
        export HOSTILE_CACHE_DIR=${hostile_cache}
        
        ${hostile} clean --offline --fastq1 ${read1} --fastq2 ${read2} --out-dir . > hostile.json
        
        mv ${read1.simpleName}.clean_1.fastq.gz ${sample_id}_R1_001_hostile.fastq.gz
        mv ${read2.simpleName}.clean_2.fastq.gz ${sample_id}_R2_001_hostile.fastq.gz
        
        echo "Hostile filtering complete" > ${sample_id}_hostile.log
        """
    } else {
        """
        export HOSTILE_CACHE_DIR=${hostile_cache}
        
        hostile clean --offline --fastq1 ${read1} --fastq2 ${read2} --out-dir . > hostile.json
        
        mv ${read1.simpleName}.clean_1.fastq.gz ${sample_id}_R1_001_hostile.fastq.gz
        mv ${read2.simpleName}.clean_2.fastq.gz ${sample_id}_R2_001_hostile.fastq.gz
        
        echo "Hostile filtering complete" > ${sample_id}_hostile.log
        """
    }
}
