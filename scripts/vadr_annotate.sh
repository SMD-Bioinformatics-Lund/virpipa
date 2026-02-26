#!/usr/bin/env bash

seq="$1"
resultsdir="$2"
id="$3"

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
vadr="apptainer -q exec -B ${vadr_bind} ${vadr_container}"

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

mkdir -p "${resultsdir}/vadr"

if [[ ! -f ${resultsdir}/vadr/${id}.vadr.seqstat ]] ; then
	${vadr} v-annotate.pl --mdir "$vadr_mdir" --group HCV --mkey flavi --forcegene "$fullseq" "${resultsdir}/vadr/"
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
			grep "5'UTR" "${tblfile/.tbl/.gff}" | sed 's/^\([^\t]*\)_.[^\t]*/\1/' | tee "${tblfile/.tbl/_mod.gff}"
			for gene in core E1 E2 p7 NS2 NS3 NS4A NS4B NS5A NS5B ; do
				target=($(grep "$gene" "${tblfile/.tbl/.gff}"))
				printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${target[0]%_pilon}" "${target[1]}" "gene" "${target[3]}" "${target[4]}" "${target[5]}" "${target[6]}" "0" "ID:${gene};gene=${gene};product=${gene} protein;" | tee -a "${tblfile/.tbl/_mod.gff}"
			done
			grep "3'UTR" "${tblfile/.tbl/.gff}" | sed 's/^\([^\t]*\)_.[^\t]*/\1/' | tee -a "${tblfile/.tbl/_mod.gff}"
		fi
	done
done
shopt -u nullglob

ns3mut=(36 54 55 56 80 122 132 155 156 168 170 175)
ns5amut=(28 30 31 58 93)
ns5bmut=(159 282 316 321 414 448 495 553 554 556)

ns3startpos=$(grep 'ID:NS3' "${resultsdir}"/vadr/"${id}"*_mod.gff | cut -f4 | head -n 1)
ns5astartpos=$(grep 'ID:NS5A' "${resultsdir}"/vadr/"${id}"*_mod.gff | cut -f4 | head -n 1)
ns5bstartpos=$(grep 'ID:NS5B' "${resultsdir}"/vadr/"${id}"*_mod.gff | cut -f4 | head -n 1)

[[ -f ${resultsdir}/vadr/${id}.vadr.bed ]] && rm "${resultsdir}/vadr/${id}.vadr.bed"

for i in "${ns3mut[@]}" ; do
	pos=$(( i*3+ns3startpos-3 ))
	printf "%s\t%s\t%s\tNS3-%s\n" "$id" "$((pos-1))" "$((pos+2))" "$i" | tee -a "${resultsdir}/vadr/${id}.vadr.bed"
done

for i in "${ns5amut[@]}" ; do
	pos=$(( i*3+ns5astartpos-3 ))
	printf "%s\t%s\t%s\tNS5A-%s\n" "$id" "$((pos-1))" "$((pos+2))" "$i" | tee -a "${resultsdir}/vadr/${id}.vadr.bed"
done

for i in "${ns5bmut[@]}" ; do
	pos=$(( i*3+ns5bstartpos-3 ))
	printf "%s\t%s\t%s\tNS5B-%s\n" "$id" "$((pos-1))" "$((pos+2))" "$i" | tee -a "${resultsdir}/vadr/${id}.vadr.bed"
done
