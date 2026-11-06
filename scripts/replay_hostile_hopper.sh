#!/bin/bash
# Replay only host filtering on unchanged raw inputs; compare paired-read payloads.
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -p normal
#SBATCH --mem=6G
#SBATCH --time=02:00:00
set -euo pipefail
if [[ $# -lt 4 || $# -gt 5 || ( $# -eq 5 && $5 != '--reorder' ) ]]; then
    echo "Usage: sbatch replay_hostile_hopper.sh READ1 READ2 PAIR_COMPARISON_SCRIPT NEW_OUTPUT_DIR [--reorder]" >&2
    exit 2
fi
reorder_flag=${5:-}
read1=$(readlink -f "$1")
read2=$(readlink -f "$2")
helper=$(readlink -f "$3")
output=$(realpath -m "$4")
for input in "$read1" "$read2" "$helper"; do
    [[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 2; }
done
[[ ! -e "$output" ]] || { echo "Output must be a new directory: $output" >&2; exit 2; }
set +u
# shellcheck source=/dev/null
source /etc/profile
set -u
module load Java/23.0.2 singularity/3.8.0
containers=/fs1/resources/containers
export HOSTILE_CACHE_DIR=/fs1/resources/ref/micro/hostile
mkdir -p "$output"
printf 'reorder\t%s\n' "${reorder_flag:-disabled}" > "$output/options.tsv"
sha256sum "$read1" "$read2" > "$output/source_reads.sha256"
sha256sum "$containers/hostile_1.1.0.sif" "$containers/python_hcvpipe.sif" > "$output/containers.sha256"
cp "${BASH_SOURCE[0]}" "$output/replay_wrapper.sh"
cp "$helper" "$output/compare_fastq_pairs.py"
for repeat in 1 2; do
    mkdir "$output/repeat_$repeat"
    (
        cd "$output/repeat_$repeat"
        singularity exec -B /fs1,/fs2,/local "$containers/hostile_1.1.0.sif" hostile \
            clean --offline ${reorder_flag:+--reorder} --fastq1 "$read1" --fastq2 "$read2" --out-dir . > hostile.json
        mapfile -t r1 < <(find . -maxdepth 1 -name '*.clean_1.fastq.gz' -print)
        mapfile -t r2 < <(find . -maxdepth 1 -name '*.clean_2.fastq.gz' -print)
        [[ ${#r1[@]} -eq 1 && ${#r2[@]} -eq 1 ]]
        mv "${r1[0]}" read1.fastq.gz
        mv "${r2[0]}" read2.fastq.gz
    )
done
sha256sum -c "$output/source_reads.sha256" > "$output/source_verification.log"
singularity exec -B /fs1,/fs2,/local "$containers/python_hcvpipe.sif" python \
    "$output/compare_fastq_pairs.py" \
    --left-r1 "$output/repeat_1/read1.fastq.gz" --left-r2 "$output/repeat_1/read2.fastq.gz" \
    --right-r1 "$output/repeat_2/read1.fastq.gz" --right-r2 "$output/repeat_2/read2.fastq.gz" \
    --outdir "$output/comparison"
touch "$output/COMPLETE"
