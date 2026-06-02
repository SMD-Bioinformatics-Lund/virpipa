process FINALIZE_RESULTS {
    tag { sample_id }
    label 'process_medium'

    cpus 2
    memory '4 GB'
    time '30m'

    input:
        tuple val(run_name), val(sample_id), val(lid), val(hostile_json_path),
            path(main_fasta), path(main_fai), path(main_blast), path(main_cram), path(main_crai),
            path(iupac_fasta), path(iupac_fai), path(iupac_blast), path(iupac_cram), path(iupac_crai), path(iupac_report), path(iupac_nucfreq),
            path(bestref_fasta), path(bestref_vcf), path(bestref_vcf_index), path(bestref_vcf_stats), path(bestref_cram), path(bestref_crai), path(bestref_report), path(bestref_nucfreq),
            path(coverage_tsv), path(vadr_gff), path(vadr_bed), path(pilon_iupac_blast),
            path(resistance_tsv), path(resistance_bed), path(resistance_gff), path(resistance_drug_tsv),
            path(filtered_vcfs), path(filtered_indices), path(filtered_stats)

    output:
        path "results", emit: results_dir
        tuple val(run_name), val(sample_id), val(lid), path("results"), emit: results_dir_with_meta

    script:
    def active_profiles = workflow.profile ?: ''
    def container_dir = params.container_dir ?: (active_profiles.contains('local_containers') ? "${projectDir}/assets/containers" : (active_profiles.contains('hpc') ? '/fs1/resources/containers' : ''))
    def bind_paths = params.bind_paths != '/fs1,/fs2,/local' ? params.bind_paths : (active_profiles.contains('local') ? '/mnt,/home,/tmp' : (active_profiles.contains('hpc') ? '/fs1,/fs2,/local,/mnt/beegfs' : params.bind_paths))
    def container_runtime = params.container_runtime ?: '$(if command -v apptainer >/dev/null 2>&1; then echo apptainer; elif command -v singularity >/dev/null 2>&1; then echo singularity; else echo apptainer; fi)'
    def scripts_dir = params.scripts_dir != 'scripts' ? params.scripts_dir : (active_profiles.contains('hpc') ? '/fs1/jonas/src/virpipa/scripts' : "${projectDir}/scripts")
    def mamba_env = env('CONDA_PREFIX') ?: '/home/jonas/miniforge3/envs/skrotis'
    def outdir_root = params.outdir.toString().startsWith('/') ? params.outdir.toString() : "${launchDir}/${params.outdir}"
    def published_results_dir = "${outdir_root}/${run_name}/${sample_id}"
    def bcftools = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/bcftools_1.21.sif bcftools" :
        "bcftools"
    def python = container_dir ?
        "${container_runtime} exec -B ${bind_paths} ${container_dir}/python_hcvpipe.sif python" :
        "python3"

    """
    set -euo pipefail

    mkdir -p results

    cp ${main_fasta} results/${sample_id}.fasta
    cp ${main_fai} results/${sample_id}.fasta.fai
    cp ${main_blast} results/${sample_id}.fasta.blast
    cp ${main_cram} results/${sample_id}.cram
    cp ${main_crai} results/${sample_id}.cram.crai

    cp ${iupac_fasta} results/${iupac_fasta.getName()}
    cp ${iupac_fai} results/${iupac_fai.getName()}
    cp ${iupac_blast} results/${iupac_blast.getName()}
    cp ${iupac_cram} results/${iupac_cram.getName()}
    cp ${iupac_crai} results/${iupac_crai.getName()}
    cp ${iupac_report} results/${iupac_report.getName()}
    cp ${iupac_nucfreq} results/${iupac_nucfreq.getName()}

    cp ${bestref_fasta} results/${bestref_fasta.getName()}
    cp ${bestref_vcf} results/${bestref_vcf.getName()}
    cp ${bestref_vcf_index} results/${bestref_vcf_index.getName()}
    cp ${bestref_vcf_stats} results/${bestref_vcf_stats.getName()}
    cp ${bestref_cram} results/${bestref_cram.getName()}
    cp ${bestref_crai} results/${bestref_crai.getName()}
    cp ${bestref_report} results/${bestref_report.getName()}
    cp ${bestref_nucfreq} results/${bestref_nucfreq.getName()}

    cp ${coverage_tsv} results/${coverage_tsv.getName()}
    cp ${vadr_gff} results/${vadr_gff.getName()}
    cp ${vadr_bed} results/${vadr_bed.getName()}
    cp ${resistance_tsv} results/${resistance_tsv.getName()}
    cp ${resistance_bed} results/${resistance_bed.getName()}
    cp ${resistance_gff} results/${resistance_gff.getName()}
    cp ${resistance_drug_tsv} results/${resistance_drug_tsv.getName()}
    printf 'query acc.ver\tsubject acc.ver\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. star\t\t\tt\ts. end\tevalue\tbit score\n' > results/${pilon_iupac_blast.getName()}

    for file in ${filtered_vcfs}; do
        cp "\${file}" "results/\$(basename "\${file}")"
    done

    for file in ${filtered_indices}; do
        cp "\${file}" "results/\$(basename "\${file}")"
    done

    for file in ${filtered_stats}; do
        cp "\${file}" "results/\$(basename "\${file}")"
    done

    if [[ -n "${hostile_json_path}" ]] && [[ -f "${hostile_json_path}" ]]; then
        cp "${hostile_json_path}" results/hostile.json
    fi

    export PATH="${mamba_env}/bin:\$PATH"

    m01_vcf=\$(find results -maxdepth 1 -name "${sample_id}-pilon-m0.01.vcf.gz" -print -quit)
    ${bcftools} query -f '[%DP\t%AD]\n' "\${m01_vcf}" | \\
        awk '{
            n = split(\$2, ad, ",")
            max_ad = ad[1]
            for (i = 2; i <= n; i++) {
                if (ad[i] > max_ad) max_ad = ad[i]
            }
            ratio = max_ad / \$1
            if (ratio <= 0.96) print ratio
        }' > ${sample_id}.mixin

    plot_title="${sample_id}"
    if [[ -n "${lid}" ]] && [[ "${lid}" != "${sample_id}" ]]; then
        plot_title="${lid} (${sample_id})"
    fi
    ${python} ${scripts_dir}/kderug.py ${sample_id}.mixin "\${plot_title}" "${sample_id}_display_rug_kde_plot.png"
    mv ${sample_id}_display_rug_kde_plot.png results/

    mkdir -p "${published_results_dir}"
    cp -R results/. "${published_results_dir}/"
    """
}
