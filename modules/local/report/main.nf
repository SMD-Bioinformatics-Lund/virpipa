process CREATE_REPORT {
    tag { sample_id }
    label 'process_low'
    
    cpus 2
    memory '4 GB'
    time '10m'
    
    beforeScript 'source /etc/profile; module load apptainer'
    
    publishDir "${params.outdir}/${run_name}/${sample_id}/results", mode: 'copy'
    
    input:
        tuple val(run_name), val(sample_id), path(vcf_stats), path(cram), path(crai), path(ref_fasta), val(subtype)
    
    output:
        path "*.report.tsv", emit: report
        path "*.fastanucfreq.tsv", emit: nucfreq
    
    script:
    def container_dir = params.container_dir
    def bind_paths = params.bind_paths ?: '/fs1,/fs2,/local'
    def samtools = "apptainer exec -B ${bind_paths} ${container_dir}/samtools_1.21.sif samtools"
    def ref_name = ref_fasta.baseName
    def sub_bash = subtype
    
    """
    echo "# VCF stats" > ${sample_id}-${ref_name}.report.tsv
    printf "subtype\t${subtype.split('-')[0]}\n" >> ${sample_id}-${ref_name}.report.tsv
    printf "reference\t${ref_name}\n" >> ${sample_id}-${ref_name}.report.tsv
    printf "id\t${sample_id}\n" >> ${sample_id}-${ref_name}.report.tsv
    
    snps=\$(grep "^SN.*number of SNPs" ${vcf_stats} | awk '{print \$NF}')
    multiallelic=\$(grep "^SN.*number of multiallelic SNP sites" ${vcf_stats} | awk '{print \$NF}')
    af0=\$(grep "^AF.*0.000000" ${vcf_stats} | awk '{print \$4}')
    af99=\$(grep "^AF.*0.990000" ${vcf_stats} | awk '{print \$4}')
    
    printf "snps\t\${snps:-0}\n" >> ${sample_id}-${ref_name}.report.tsv
    printf "multiallelic_snps\t\${multiallelic:-0}\n" >> ${sample_id}-${ref_name}.report.tsv
    printf "AF0\t\${af0:-0}\n" >> ${sample_id}-${ref_name}.report.tsv
    printf "AF0.99\t\${af99:-0}\n" >> ${sample_id}-${ref_name}.report.tsv
    
    printf "# COVERAGE\n" >> ${sample_id}-${ref_name}.report.tsv
    ${samtools} coverage ${cram} > coverage.tmp
    awk 'NR>1 {print}' coverage.tmp >> ${sample_id}-${ref_name}.report.tsv
    rm coverage.tmp
    
    awk -vFS="" 'NR>1 {for(i=1;i<=NF;i++)w[toupper(\$i)]++}END{for(i in w) print i,w[i]}' ${ref_fasta} | sort -nr -k2 > ${sample_id}-${ref_name}.fastanucfreq.tsv
    """
}
