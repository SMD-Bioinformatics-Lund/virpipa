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
        tuple val(run_name), val(sample_id), path("hostile.json"), emit: hostile_json_with_meta
        path "*.log", emit: logs
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    def hostile_cache = cache_dir ?: params.hostile_cache_dir ?: '/fs1/resources/ref/micro/hostile'
    
    if (container_dir) {
        def hostile = "${container_runtime} exec -B ${bind_paths} ${container_dir}/hostile_1.1.0.sif hostile"
        
        """
        set -euo pipefail

        export HOSTILE_CACHE_DIR=${hostile_cache}
        
        ${hostile} clean --offline --fastq1 ${read1} --fastq2 ${read2} --out-dir . > hostile.json
        
        clean_r1=\$(find . -maxdepth 1 -name '*.clean_1.fastq.gz' -print -quit)
        clean_r2=\$(find . -maxdepth 1 -name '*.clean_2.fastq.gz' -print -quit)

        mv "\${clean_r1}" ${sample_id}_R1_001_hostile.fastq.gz
        mv "\${clean_r2}" ${sample_id}_R2_001_hostile.fastq.gz
        
        echo "Hostile filtering complete" > ${sample_id}_hostile.log
        """
    } else {
        """
        set -euo pipefail

        export HOSTILE_CACHE_DIR=${hostile_cache}
        
        hostile clean --offline --fastq1 ${read1} --fastq2 ${read2} --out-dir . > hostile.json
        
        clean_r1=\$(find . -maxdepth 1 -name '*.clean_1.fastq.gz' -print -quit)
        clean_r2=\$(find . -maxdepth 1 -name '*.clean_2.fastq.gz' -print -quit)

        mv "\${clean_r1}" ${sample_id}_R1_001_hostile.fastq.gz
        mv "\${clean_r2}" ${sample_id}_R2_001_hostile.fastq.gz
        
        echo "Hostile filtering complete" > ${sample_id}_hostile.log
        """
    }
}
