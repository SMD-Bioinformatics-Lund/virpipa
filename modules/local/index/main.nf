process INDEX_REFERENCE {
    tag { genome }
    label 'process_low'
    
    cpus 4
    memory '8 GB'
    time '1h'
    
    publishDir "${params.ref_dir}", mode: 'copy', saveAs: { filename -> "indexed/$filename" }
    
    input:
        path(genome)
    
    output:
        path("indexed/${genome}"), emit: genome
        path("indexed/${genome}.*"), emit: index
        path("indexed/*.fai"), emit: fai
    
    script:
    """
    mkdir -p indexed
    cp ${genome} indexed/
    
    bwa indexed/${genome}
    samtools faidx indexed/${genome}
    """
}
