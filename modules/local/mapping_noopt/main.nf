process MAP_READS_NOOPT {
    tag { "${sample_id}:${genome_name}" }

    cpus { params.use_sentieon ? 16 : 8 }
    memory { params.use_sentieon ? '32 GB' : '8 GB' }
    time '8h'

    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bai'
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.stats'

    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2), path(genome), val(genome_name)

    output:
        tuple val(run_name), val(sample_id), path("${sample_id}-${genome_name}.r11b2L25.bwa.umi.sort.bam"), path("${sample_id}-${genome_name}.r11b2L25.bwa.umi.sort.bam.bai"), emit: raw_bams
        tuple val(run_name), val(sample_id), path("${sample_id}-${genome_name}.r11b2L25.bwa.umi.filter.sort.bam"), path("${sample_id}-${genome_name}.r11b2L25.bwa.umi.filter.sort.bam.bai"), emit: bams
        tuple val(run_name), val(sample_id), path("${sample_id}-${genome_name}.r11b2L25.bwa.umi.filter.sort.bam.stats"), emit: stats

    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def use_sentieon = params.use_sentieon ?: false
    def umi_prefix = "${sample_id}-${genome_name}.r11b2L25.bwa.umi"

    if (use_sentieon && container_dir) {
        def sentieon = "apptainer exec -B ${bind_paths} ${container_dir}/sentieon_202308.03.sif sentieon"
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"

        """
        set -euo pipefail

        ref_path=\$(readlink -f ${genome})
        ref_dir=\$(dirname \${ref_path})
        ref_base=\$(basename \${ref_path})
        mkdir -p ref_copy
        cp -L \${ref_path} ref_copy/
        cp \${ref_dir}/\${ref_base}.* ref_copy/ 2>/dev/null || true
        if [[ ! -f ref_copy/\${ref_base}.bwt ]]; then
            ${sentieon} bwa index ref_copy/\${ref_base}
        fi

        ${sentieon} umi extract -d 3M2S+T,3M2S+T ${read1} ${read2} | \\
        ${sentieon} bwa mem \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            -p -C ref_copy/\${ref_base} - | \\
        ${sentieon} umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz

        if [[ ! -f consensus.fastq.gz ]] && compgen -G 'consensus.fastq.gz_*.gz' > /dev/null; then
            cat consensus.fastq.gz_*.gz > consensus.fastq.gz
        fi

        ${sentieon} bwa mem \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            -p -C ref_copy/\${ref_base} consensus.fastq.gz | \\
        ${sentieon} util sort -i - --sam2bam --umi_post_process -o ${umi_prefix}.sort.bam

        ${samtools} view -@ ${task.cpus} -e "sclen < 30" --with-header -b -o ${umi_prefix}.filter.sort.bam ${umi_prefix}.sort.bam
        ${samtools} index ${umi_prefix}.sort.bam
        ${samtools} index ${umi_prefix}.filter.sort.bam
        ${samtools} stats ${umi_prefix}.filter.sort.bam > ${umi_prefix}.filter.sort.bam.stats
        """
    } else {
        """
        echo "MAP_READS_NOOPT requires sentieon container support. Set params.use_sentieon=true and provide container_dir."
        exit 1
        """
    }
}
