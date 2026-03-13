process POLISH_PILON_LOOP {
    tag { "${sample_id}" }
    label 'process_high'
    
    cpus 16
    memory '32 GB'
    time '24h'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/pilon", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(read1), path(read2), path(hybrid_assembly)
    
    output:
        tuple val(run_name), val(sample_id), path("*.fasta"), path("*.fai"), emit: polished
        path "*.bam", emit: bams
        path "*.bai", emit: bais
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def use_sentieon = params.use_sentieon ?: true
    
    if (use_sentieon && container_dir) {
        def sentieon = "apptainer exec -B ${bind_paths} ${container_dir}/sentieon_202308.03.sif sentieon"
        def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
        def pilon = "apptainer exec -B ${bind_paths} ${container_dir}/pilon-1.24.sif pilon"
        
        def run_name = run_name
        def sample = sample_id
        def maxpolish = 10
        
        """
        # Write shell script to file
        cat > polish.sh << 'SCRIPT'
        #!/bin/bash
        set -e
        
        hybrid_fa="${hybrid_assembly}"
        sample="${sample}"
        run_name="${run_name}"
        outdir="${params.outdir}/${run_name}/${sample}/pilon"
        maxpolish=${maxpolish}
        
        mkdir -p "\${outdir}"
        
        map_to_ref() {
            local ref=\$1
            local name=\$2
            
            ref_base=\$(basename \${ref})
            mkdir -p ref_\${name}
            cp -L \${ref} ref_\${name}/
            
            ${sentieon} umi extract -d 3M2S+T,3M2S+T ${read1} ${read2} | \\
            ${sentieon} bwa mem \\
                -R "@RG\\tID:\${sample}\\tSM:\${sample}\\tLB:\${sample}\\tPL:illumina" \\
                -t ${task.cpus} \\
                -k 11 -B 2 -L 25 \\
                -p -C ref_\${name}/\${ref_base} - | \\
            ${sentieon} umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz
            
            ${sentieon} bwa mem \\
                -R "@RG\\tID:\${sample}\\tSM:\${sample}\\tLB:\${sample}\\tPL:illumina" \\
                -t ${task.cpus} \\
                -k 11 -B 2 -L 25 \\
                -p -C ref_\${name}/\${ref_base} consensus.fastq.gz | \\
            ${sentieon} util sort -i - --sam2bam --umi_post_process -o \${sample}-\${name}.bam
            
            ${samtools} index \${sample}-\${name}.bam
            echo "\${sample}-\${name}.bam"
        }
        
        run_pilon() {
            local ref=\$1
            local input_bam=\$2
            local round=\$3
            
            ${samtools} view -f 0x2 \${input_bam} > paired-\${round}.bam
            ${samtools} view -F 0x2 \${input_bam} > unpaired-\${round}.bam
            ${samtools} index paired-\${round}.bam
            ${samtools} index unpaired-\${round}.bam
            
            for iupac in "" "--iupac" ; do
                iupac_suffix=\${iupac:1}
                ${pilon} \${iupac} --fix all --mindepth 5 --changes \\
                    --genome \${ref} \\
                    --frags paired-\${round}.bam \\
                    --unpaired unpaired-\${round}.bam \\
                    --outdir \${outdir} \\
                    --output \${sample}-pilon-\${round}\${iupac_suffix}
                
                sed -i 's/_pilon//' \${outdir}/\${sample}-pilon-\${round}\${iupac_suffix}.fasta
                ${samtools} faidx \${outdir}/\${sample}-pilon-\${round}\${iupac_suffix}.fasta
                ${sentieon} bwa index \${outdir}/\${sample}-pilon-\${round}\${iupac_suffix}.fasta
            done
            
            echo "\${outdir}/\${sample}-pilon-\${round}.fasta"
        }
        
        # Round 1
        echo "=== POLISH ROUND 1 ==="
        bam1=\$(map_to_ref "\${hybrid_fa}" "pilon-1")
        ref1=\$(run_pilon "\${hybrid_fa}" "\${bam1}" 1)
        prev_changes=\$(wc -l < \${outdir}/\${sample}-pilon-1.changes)
        
        # Rounds 2-10
        for round in \$(seq 2 \${maxpolish}); do
            prev_ref="\${outdir}/\${sample}-pilon-\$((round-1)).fasta"
            bam=\$(map_to_ref "\${prev_ref}" "pilon-\${round}")
            current_ref=\$(run_pilon "\${prev_ref}" "\${bam}" \${round})
            
            curr_changes=\$(wc -l < \${outdir}/\${sample}-pilon-\${round}.changes)
            prev_changes_file=\$(wc -l < \${outdir}/\${sample}-pilon-\$((round-1)).changes)
            
            if [[ "\${curr_changes}" -eq "\${prev_changes_file}" ]] || [[ \${round} -eq \${maxpolish} ]]; then
                cp \${outdir}/\${sample}-pilon-\${round}.fasta \${outdir}/\${sample}.fasta
                cp \${outdir}/\${sample}-pilon-\${round}-iupac.fasta \${outdir}/\${sample}-iupac.fasta
                ${samtools} faidx \${outdir}/\${sample}.fasta
                break
            fi
        done
        
        # Final mapping
        echo "=== Final mapping ==="
        final_ref="\${outdir}/\${sample}.fasta"
        final_bam=\$(map_to_ref "\${final_ref}" "pilon")
        
        cp \${sample}-pilon.bam \${outdir}/\${sample}-pilon.bam
        cp \${sample}-pilon.bam.bai \${outdir}/\${sample}-pilon.bam.bai
        cp \${outdir}/\${sample}.fasta ./
        cp \${outdir}/\${sample}.fasta.fai ./
        SCRIPT
        
        chmod +x polish.sh
        ./polish.sh
        """
    } else {
        error "POLISH_PILON_LOOP requires use_sentieon=true and container_dir"
    }
}
