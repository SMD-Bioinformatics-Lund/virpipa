process CREATE_REPORT {
    tag { sample_id }
    label 'process_low'
    
    cpus 2
    memory '4 GB'
    time '10m'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf_stats), path(cram), path(crai), path(ref_fasta), val(subtype)
    
    output:
        path "*.report.tsv", emit: report
        path "*.fastanucfreq.tsv", emit: nucfreq

    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def report_id = vcf_stats.getName().replaceFirst(/\.vcf\.gz\.stats$/, '')
    def samtools = container_dir ?
        "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools" :
        "samtools"
    
    """
    set -euo pipefail

    snps=\$(sed -n 's/^SN\t0\tnumber of SNPs:\t\\([0-9]*\\)\$/\\1/p' ${vcf_stats})
    multiallelic=\$(sed -n 's/^SN\t0\tnumber of multiallelic SNP sites:\t\\([0-9]*\\)\$/\\1/p' ${vcf_stats})
    af0=\$(sed -n 's/^AF\t0\t0.000000\t\\([0-9]*\\)\t.*/\\1/p' ${vcf_stats})
    af99=\$(sed -n 's/^AF\t0\t0.990000\t\\([0-9]*\\)\t.*/\\1/p' ${vcf_stats})

    echo '# VCF stats' > ${report_id}.report.tsv
    printf 'subtype\t%s\n' '${subtype.split('-')[0]}' >> ${report_id}.report.tsv
    printf 'reference\t%s\n' '${subtype}' >> ${report_id}.report.tsv
    printf 'id\t%s\n' '${sample_id}' >> ${report_id}.report.tsv
    printf 'snps\t%s\n' "\${snps}" >> ${report_id}.report.tsv
    printf 'multiallelic_snps\t%s\n' "\${multiallelic}" >> ${report_id}.report.tsv
    printf 'AF0\t%s\n' "\${af0}" >> ${report_id}.report.tsv
    printf 'AF0.99\t%s\n' "\${af99}" >> ${report_id}.report.tsv

    echo '# COVERAGE' >> ${report_id}.report.tsv
    ${samtools} coverage ${cram} | awk 'NR==1 { sub(/^#/, "", \$1); for (i=1; i<=NF; i++) headers[i]=\$i; next } NR==2 { for (i=2; i<=NF; i++) print headers[i] "\\t" \$i }' >> ${report_id}.report.tsv

    awk -vFS="" 'NR>1 {for(i=1;i<=NF;i++)w[toupper(\$i)]++}END{for(i in w) print i,w[i]}' ${ref_fasta} | sort -nr -k2 > ${report_id}.fastanucfreq.tsv
    """
}
