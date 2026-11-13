#!/bin/bash
# Fresh, isolated Hopper baseline. Run through sbatch, not the login node.
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH -p normal
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch scripts/runme_consensus_benchmark.sh CAMPAIGN_CSV NEW_OUTPUT_DIR" >&2
    exit 2
fi
input=$(readlink -f "$1")
output=$(realpath -m "$2")
source_root=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
source_commit=$(cat "$source_root/.git/HEAD")
[[ "$source_commit" =~ ^[0-9a-f]{40}$ ]] || {
    echo 'Submit from a detached, pinned clone of the repository' >&2
    exit 2
}
[[ -f "$input" && ! -e "$output" ]]
[[ "$output" == /fs1/jonas/virpipa/consensus_benchmark_* ]] || {
    echo 'Benchmark outputs must use a new /fs1/jonas/virpipa/consensus_benchmark_* directory' >&2
    exit 2
}

set +u
source /etc/profile
set -u
module load Java/23.0.2 nextflow/26.04.3 apptainer
export NXF_OFFLINE=true
export NXF_SYNTAX_PARSER=v2

mkdir -p "$output/launch"
cp "$input" "$output/input.csv"
cp "${BASH_SOURCE[0]}" "$output/replay_wrapper.sh"
printf '%s\n' "$source_commit" > "$output/source_commit.txt"
sha256sum "$input" "$output/input.csv" > "$output/samplesheet.sha256"
nextflow -version > "$output/nextflow.version.txt"
nextflow config "$source_root" -profile slurm,hpc,apptainer > "$output/profile-defaults.config"
apptainer exec -B /fs1,/fs2,/local /fs1/resources/containers/python_hcvpipe.sif python \
    "$source_root/scripts/benchmark_campaign.py" --csv "$output/input.csv" \
    --output "$output/campaign.json"

cd "$output/launch"
nextflow run "$source_root/main.nf" -profile slurm,hpc,apptainer \
    --input "$output/input.csv" --dev --outdir "$output/results" \
    --ref_dir "$source_root/refgenomes" --partition normal \
    --completion_log_dir= --virtitta_import_dir= \
    --publish_mode debug --pipeline_info true --pipeline_info_dag false \
    -work-dir "$output/work" -ansi-log false
sha256sum -c "$output/samplesheet.sha256" > "$output/input_verification.log"
touch "$output/BASELINE_COMPLETE"
