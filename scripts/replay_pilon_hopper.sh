#!/bin/bash
# Fixed-input diagnostic only: no remapping and no writes to published results.
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p normal
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch replay_pilon_hopper.sh GENOME PAIRED_BAM UNPAIRED_BAM NEW_OUTPUT_DIR" >&2
    exit 2
fi
genome=$(readlink -f "$1")
paired=$(readlink -f "$2")
unpaired=$(readlink -f "$3")
output=$(realpath -m "$4")
for input in "$genome" "$paired" "$unpaired" "${paired}.bai" "${unpaired}.bai"; do
    [[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 2; }
done
[[ ! -e "$output" ]] || { echo "Output must be a new directory: $output" >&2; exit 2; }

set +u
# shellcheck source=/dev/null
source /etc/profile
set -u
module load Java/23.0.2 singularity/3.8.0
container=/fs1/resources/containers/pilon-1.24.sif
[[ -f "$container" ]] || { echo "Missing Pilon container" >&2; exit 2; }

mkdir -p "$output/inputs"
frozen="$output/inputs"
cp "$genome" "$frozen/genome.fasta"
if [[ -f "${genome}.fai" ]]; then
    cp "${genome}.fai" "$frozen/genome.fasta.fai"
fi
cp "$paired" "$frozen/paired.bam"
cp "${paired}.bai" "$frozen/paired.bam.bai"
cp "$unpaired" "$frozen/unpaired.bam"
cp "${unpaired}.bai" "$frozen/unpaired.bam.bai"
printf 'genome\t%s\npaired\t%s\nunpaired\t%s\n' "$genome" "$paired" "$unpaired" > "$output/source_paths.tsv"
sha256sum "$container" > "$output/container.sha256"
sha256sum "$frozen"/* > "$output/frozen_inputs.sha256"
cp "${BASH_SOURCE[0]}" "$output/replay_wrapper.sh"

printf 'repeat\tsequence_sha256\tchange_count\tchanges_sha256\n' > "$output/repeats.tsv"
for repeat in 1 2 3; do
    destination="$output/repeat_$repeat"
    mkdir "$destination"
    singularity exec -B /fs1,/fs2,/local "$container" pilon \
        --fix all --mindepth 5 --changes \
        --genome "$frozen/genome.fasta" \
        --frags "$frozen/paired.bam" --unpaired "$frozen/unpaired.bam" \
        --outdir "$destination" --output polished > "$destination/pilon.log" 2>&1
    sequence_hash=$(awk '!/^>/ {printf "%s", toupper($0)}' "$destination/polished.fasta" | sha256sum | cut -d' ' -f1)
    changes_hash=$(sha256sum "$destination/polished.changes" | cut -d' ' -f1)
    count=$(wc -l < "$destination/polished.changes")
    printf '%s\t%s\t%s\t%s\n' "$repeat" "$sequence_hash" "$count" "$changes_hash" >> "$output/repeats.tsv"
done
sha256sum -c "$output/frozen_inputs.sha256" > "$output/input_verification.log"
cat "$output/repeats.tsv"
touch "$output/COMPLETE"
