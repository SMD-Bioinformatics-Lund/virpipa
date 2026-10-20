#!/usr/bin/env bash
set -euo pipefail

seq="${1:-}"
resultsdir="${2:-}"
id="${3:-}"

if [[ -z "$seq" ]] || [[ -z "$resultsdir" ]] ; then
	echo "Usage: $(basename "$0") <fasta> <results_dir> [sample_id]"
	exit 1
fi

scriptdir=$(cd "$(dirname "$0")" && pwd)
fullseq=$(readlink -f "$seq")
shortseq=$(basename "$seq")
if [[ -z "$id" ]] ; then
	id=${shortseq%.fasta}
fi

vadr_container="${VADR_CONTAINER:-/fs1/resources/containers/vadr_164.sif}"
vadr_bind="${VADR_BIND:-/fs1:/fs1}"
vadr_mdir="${VADR_MODELDIR:-/fs1/resources/ref/micro/vadr/vadr-models-flavi/}"
annotate_tbl2gff="${VADR_ANNOTATE_TBL2GFF:-${scriptdir}/annotate-tbl2gff.pl}"
container_runtime="${VADR_CONTAINER_RUNTIME:-}"
if [[ -z "$container_runtime" ]]; then
    if command -v apptainer >/dev/null 2>&1; then
        container_runtime=apptainer
    elif command -v singularity >/dev/null 2>&1; then
        container_runtime=singularity
    else
        echo "Neither apptainer nor singularity found in PATH" >&2
        exit 127
    fi
fi
vadr="${container_runtime} -q exec -B ${vadr_bind} ${vadr_container}"

if [[ ! -d "$resultsdir" ]] ; then
	echo "$resultsdir does not exist"
	exit 1
fi
if [[ ! -f "$fullseq" ]] ; then
	echo "$fullseq does not exist"
	exit 1
fi
if [[ ! -f "$annotate_tbl2gff" ]] ; then
	echo "annotate-tbl2gff.pl not found at $annotate_tbl2gff"
	exit 1
fi

if [[ ! -f ${resultsdir}/vadr/${id}.vadr.seqstat ]] ; then
	${vadr} v-annotate.pl -f --mdir "$vadr_mdir" --group HCV --mkey flavi --forcegene "$fullseq" "${resultsdir}/vadr/"
fi

shopt -s nullglob
for i in "${resultsdir}"/vadr/vadr.vadr* ; do
	mv "$i" "${i/vadr.vadr/${id}.vadr}"
done

for status in fail pass ; do
	for tblfile in "${resultsdir}"/vadr/"${id}".vadr."${status}".tbl ; do
		echo "$tblfile"
		if [[ -s "$tblfile" ]] ; then
			echo "${vadr} ${annotate_tbl2gff} ${tblfile}"
			${vadr} "$annotate_tbl2gff" "$tblfile" > "${tblfile/.tbl/.gff}"
			awk -f "${scriptdir}/normalize_vadr_gff.awk" "${tblfile/.tbl/.gff}" > "${tblfile/.tbl/_mod.gff}"
		fi
	done
done
shopt -u nullglob

ns3mut=(36 54 55 56 80 122 132 155 156 168 170 175)
ns5amut=(28 30 31 58 93)
ns5bmut=(159 282 316 321 414 448 495 553 554 556)

annotations=("${resultsdir}"/vadr/"${id}"*_mod.gff)
if [[ ! -f "${annotations[0]}" ]]; then
    echo "VADR produced no gene annotation for ${id}" >&2
    exit 1
fi

# Create vadr subdirectory and bed file
mkdir -p "${resultsdir}/vadr"
> "${resultsdir}/vadr/${id}.vadr.bed"

# This legacy track is gene-relative, not H77 numbering. Skip partial starts.
awk 'BEGIN { FS=OFS="\t";
    sites["NS3"]="36 54 55 56 80 122 132 155 156 168 170 175";
    sites["NS5A"]="28 30 31 58 93";
    sites["NS5B"]="159 282 316 321 414 448 495 553 554 556" }
    $3 == "gene" && $9 !~ /trunc5=1/ {
        for (gene in sites) if ($9 ~ "gene=" gene ";") {
            n=split(sites[gene], positions, " ");
            for (i=1; i<=n; i++) {
                start=$4+(positions[i]-1)*3;
                if (start+2 <= $5) print $1,start-1,start+2,gene "-" positions[i]
            }
        }
    }' "${annotations[@]}" > "${resultsdir}/vadr/${id}.vadr.bed"
