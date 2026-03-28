process POLISH_PILON_LOOP {
    tag { "${sample_id}" }
    label 'process_high'

    cpus 16
    memory '32 GB'
    time '24h'

    publishDir "${params.outdir}/${run_name}/${sample_id}/pilon", mode: 'copy', pattern: '*.fasta'
    publishDir "${params.outdir}/${run_name}/${sample_id}/pilon", mode: 'copy', pattern: '*.fasta.*'
    publishDir "${params.outdir}/${run_name}/${sample_id}/pilon", mode: 'copy', pattern: '*.changes'
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bam.bai'
    publishDir "${params.outdir}/${run_name}/${sample_id}/bam", mode: 'copy', pattern: '*.bam.stats'

    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2), path(hybrid_assembly)

    output:
        tuple val(run_name), val(sample_id), path("${sample_id}.fasta"), path("${sample_id}.fasta.fai"), emit: polished
        tuple val(run_name), val(sample_id), path("${sample_id}-pilon.r11b2L25.bwa.umi.filter.sort.bam"), path("${sample_id}-pilon.r11b2L25.bwa.umi.filter.sort.bam.bai"), emit: final_bam_with_index
        tuple val(run_name), val(sample_id), path("${sample_id}-pilon-iupac.fasta"), emit: polished_pilon_iupac
        path("${sample_id}-iupac.fasta")
        path("${sample_id}-iupac.fasta.fai")
        path("${sample_id}-pilon-iupac.fasta")
        path("${sample_id}-pilon-*.fasta*")
        path("${sample_id}-pilon*.changes")
        path("${sample_id}-pilon*.bam")
        path("${sample_id}-pilon*.bam.bai")
        path("${sample_id}-pilon*.bam.stats")

    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def use_sentieon = params.use_sentieon ?: true
    def maxpolish = params.maxpolish ?: 10

    if (!(use_sentieon && container_dir)) {
        error "POLISH_PILON_LOOP requires use_sentieon=true and container_dir"
    }

    def sentieon = "apptainer exec -B ${bind_paths} ${container_dir}/sentieon_202308.03.sif sentieon"
    def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
    def pilon = "apptainer exec -B ${bind_paths} ${container_dir}/pilon-1.24.sif pilon"

    """
    set -euo pipefail

    sample='${sample_id}'
    maxpolish='${maxpolish}'
    hybrid_ref=\$(readlink -f ${hybrid_assembly})

    stage_reference() {
        local source_ref=\$1
        local staged_dir=\$2
        local source_abs
        local source_dir
        local source_name

        source_abs=\$(readlink -f "\${source_ref}")
        source_dir=\$(dirname "\${source_abs}")
        source_name=\$(basename "\${source_abs}")

        mkdir -p "\${staged_dir}"
        cp -L "\${source_abs}" "\${staged_dir}/\${source_name}"
        cp "\${source_dir}/\${source_name}".* "\${staged_dir}/" 2>/dev/null || true

        if [[ ! -f "\${staged_dir}/\${source_name}.bwt" ]]; then
            ${sentieon} bwa index "\${staged_dir}/\${source_name}"
        fi

        printf '%s\n' "\${staged_dir}/\${source_name}"
    }

    run_mapping() {
        local ref=\$1
        local prefix=\$2
        local use_opt=\$3
        local staged_ref
        local bwa_args=()

        staged_ref=\$(stage_reference "\${ref}" "ref_\${prefix}")

        if [[ "\${use_opt}" == "true" ]]; then
            bwa_args=(-k 11 -B 2 -L 25)
        fi

        ${sentieon} umi extract -d 3M2S+T,3M2S+T ${read1} ${read2} | \\
        ${sentieon} bwa mem \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            "\${bwa_args[@]}" \\
            -p -C "\${staged_ref}" - | \\
        ${sentieon} umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz

        if [[ ! -f consensus.fastq.gz ]] && compgen -G 'consensus.fastq.gz_*.gz' > /dev/null; then
            cat consensus.fastq.gz_*.gz > consensus.fastq.gz
        fi

        ${sentieon} bwa mem \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            "\${bwa_args[@]}" \\
            -p -C "\${staged_ref}" consensus.fastq.gz | \\
        ${sentieon} util sort -i - --sam2bam --umi_post_process -o "\${prefix}.sort.bam"

        ${samtools} view -@ ${task.cpus} -e "sclen < 30" --with-header -b -o "\${prefix}.filter.sort.bam" "\${prefix}.sort.bam"
        ${samtools} index "\${prefix}.sort.bam"
        ${samtools} index "\${prefix}.filter.sort.bam"
        ${samtools} stats "\${prefix}.filter.sort.bam" > "\${prefix}.filter.sort.bam.stats"
    }

    run_pilon_round() {
        local ref=\$1
        local round=\$2
        local prefix="\${sample}-pilon-\${round}.r11b2L25.bwa.umi"

        run_mapping "\${ref}" "\${prefix}" true

        ${samtools} view -f 0x2 -b "\${prefix}.filter.sort.bam" > "paired-\${round}.bam"
        ${samtools} view -F 0x2 -b "\${prefix}.filter.sort.bam" > "unpaired-\${round}.bam"
        ${samtools} index "paired-\${round}.bam"
        ${samtools} index "unpaired-\${round}.bam"

        for iupac in "" "--iupac"; do
            local suffix=""
            if [[ -n "\${iupac}" ]]; then
                suffix="-iupac"
            fi

            ${pilon} \${iupac} --fix all --mindepth 5 --changes \\
                --genome "\${ref}" \\
                --frags "paired-\${round}.bam" \\
                --unpaired "unpaired-\${round}.bam" \\
                --outdir . \\
                --output "\${sample}-pilon-\${round}\${suffix}"

            sed -i 's/_pilon//' "\${sample}-pilon-\${round}\${suffix}.fasta"
            ${samtools} faidx "\${sample}-pilon-\${round}\${suffix}.fasta"
            ${sentieon} bwa index "\${sample}-pilon-\${round}\${suffix}.fasta"
        done
    }

    run_pilon_round "\${hybrid_ref}" 1

    for round in \$(seq 2 "\${maxpolish}"); do
        prev_round=\$((round - 1))
        run_pilon_round "\${sample}-pilon-\${prev_round}.fasta" "\${round}"

        curr_changes=\$(wc -l < "\${sample}-pilon-\${round}.changes")
        prev_changes=\$(wc -l < "\${sample}-pilon-\${prev_round}.changes")

        if [[ "\${curr_changes}" -eq "\${prev_changes}" ]] || [[ "\${round}" -eq "\${maxpolish}" ]]; then
            cp "\${sample}-pilon-\${round}.fasta" "\${sample}.fasta"
            cp "\${sample}-pilon-\${round}-iupac.fasta" "\${sample}-iupac.fasta"
            cp "\${sample}-iupac.fasta" "\${sample}-pilon-iupac.fasta"
            ${samtools} faidx "\${sample}.fasta"
            ${samtools} faidx "\${sample}-iupac.fasta"
            ${sentieon} bwa index "\${sample}.fasta"
            ${sentieon} bwa index "\${sample}-iupac.fasta"
            break
        fi
    done

    if [[ ! -f "\${sample}.fasta" ]]; then
        last_round=\$(find . -maxdepth 1 -name "\${sample}-pilon-[0-9]*.changes" ! -name "*-iupac.changes" -printf '%f\n' | sed -E 's/.*-([0-9]+)\\.changes/\\1/' | sort -n | tail -1)
        cp "\${sample}-pilon-\${last_round}.fasta" "\${sample}.fasta"
        cp "\${sample}-pilon-\${last_round}-iupac.fasta" "\${sample}-iupac.fasta"
        cp "\${sample}-iupac.fasta" "\${sample}-pilon-iupac.fasta"
        ${samtools} faidx "\${sample}.fasta"
        ${samtools} faidx "\${sample}-iupac.fasta"
        ${sentieon} bwa index "\${sample}.fasta"
        ${sentieon} bwa index "\${sample}-iupac.fasta"
    fi

    run_mapping "\${sample}.fasta" "\${sample}-pilon.r11b2L25.bwa.umi" false
    """
}
