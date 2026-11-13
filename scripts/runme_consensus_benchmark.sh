#!/bin/bash
# Fresh, isolated Hopper baseline. Run through sbatch, not the login node.
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH -p normal
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: sbatch scripts/runme_consensus_benchmark.sh INPUT OUTPUT_DIR [baseline|native] [resume]" >&2
    exit 2
fi
mode=${3:-baseline}
recovery=${4:-}
[[ "$mode" == baseline || "$mode" == native ]]
[[ -z "$recovery" || "$recovery" == resume ]]
input=$(readlink -f "$1")
output=$(realpath -m "$2")
source_root=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
source_commit=$(cat "$source_root/.git/HEAD")
[[ "$source_commit" =~ ^[0-9a-f]{40}$ ]] || {
    echo 'Submit from a detached, pinned clone of the repository' >&2
    exit 2
}
[[ -f "$input" ]]
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

resume_args=()
if [[ "$recovery" == resume ]]; then
    [[ -d "$output/launch" && -d "$output/work" ]]
    [[ "$(cat "$output/source_commit.txt")" == "$source_commit" ]]
    cmp "$input" "$output/input.csv"
    [[ ! -e "$output/${mode^^}_COMPLETE" ]] || { echo 'Experiment already completed'; exit 0; }
    resume_args=(-resume)
    cp "${BASH_SOURCE[0]}" "$output/resume_wrapper_${SLURM_JOB_ID:-manual}.sh"
    if [[ "$mode" == baseline ]]; then
        apptainer exec -B /fs1,/fs2,/local /fs1/resources/containers/python_hcvpipe.sif \
            python - "$source_root/scripts" "$output/campaign.json" <<'PY'
import json
import sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, sys.argv[1])
from benchmark_campaign import fingerprint
with open(sys.argv[2]) as stream:
    expected = json.load(stream)['inputs']
with ThreadPoolExecutor(max_workers=4) as pool:
    actual = list(pool.map(fingerprint, [entry['path'] for entry in expected]))
for old, new in zip(expected, actual):
    if old['sha256'] != new['sha256'] or old['size'] != new['size']:
        raise ValueError('Frozen raw input changed: '+old['path'])
print('Verified frozen raw inputs before resume:', len(actual))
PY
    fi
else
    [[ ! -e "$output" ]]
    mkdir -p "$output/launch"
    cp "$input" "$output/input.csv"
    cp "${BASH_SOURCE[0]}" "$output/replay_wrapper.sh"
    printf '%s\n' "$source_commit" > "$output/source_commit.txt"
    sha256sum "$input" "$output/input.csv" > "$output/samplesheet.sha256"
    nextflow -version > "$output/nextflow.version.txt"
    nextflow config "$source_root" -profile slurm,hpc,apptainer > "$output/profile-defaults.config"
fi
entry=main.nf
if [[ "$mode" == baseline ]]; then
    if [[ "$recovery" != resume ]]; then
        apptainer exec -B /fs1,/fs2,/local /fs1/resources/containers/python_hcvpipe.sif python \
            "$source_root/scripts/benchmark_campaign.py" --csv "$output/input.csv" \
            --output "$output/campaign.json"
    fi
else
    entry=benchmark_native.nf
fi

cd "$output/launch"
nextflow run "$source_root/$entry" -profile slurm,hpc,apptainer \
    --input "$output/input.csv" --dev --outdir "$output/results" \
    --ref_dir "$source_root/refgenomes" --partition normal \
    --completion_log_dir= --virtitta_import_dir= \
    --publish_mode debug --pipeline_info true --pipeline_info_dag false \
    -work-dir "$output/work" -ansi-log false "${resume_args[@]}"
sha256sum -c "$output/samplesheet.sha256" > "$output/input_verification.log"
touch "$output/${mode^^}_COMPLETE"
