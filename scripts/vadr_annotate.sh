#!/usr/bin/env bash

seq="$1"
resultsdir="$2"
id="$3"

fullseq=$(readlink -f $seq)
shortseq=$(basename $seq)
id=${shortseq%.fasta}
vadr='apptainer -q exec -B /fs1:/fs1 /fs1/resources/containers/vadr_164.sif'


if [[ ! -d $resultsdir ]] ; then
	echo $resultsdir does not exist
	exit
fi

# [[ ! -d ${resultsdir}/vadr ]] && mkdir $resultsdir/vadr

if [[ ! -f ${resultsdir}/vadr/${id}.vadr.seqstat ]] ; then
	${vadr} v-annotate.pl --mdir /fs1/resources/ref/micro/vadr/vadr-models-flavi/ --group HCV --mkey flavi --forcegene $fullseq ${resultsdir}/vadr/
fi

for i in ${resultsdir}/vadr/vadr.vadr* ; do
	mv $i ${i/vadr.vadr/${id}.vadr}
done

for status in fail pass ; do
	for tblfile in ${resultsdir}/vadr/${id}.vadr.${status}.tbl ; do
echo $tblfile
		if [[ -s $tblfile ]] ; then
#			ls $tblfile
#			echo ${tblfile/.tbl/.gff}
#			echo 		${vadr} /data/bnf/dev/jonas/hcv/VADR/vadr/vadr-vadr-1.7/miniscripts/annotate-tbl2gff.pl $tblfile
			echo ${vadr} /fs1/jonas/src/virpipa/scripts/annotate-tbl2gff.pl $tblfile
			${vadr} /fs1/jonas/src/virpipa/scripts/annotate-tbl2gff.pl $tblfile > ${tblfile/.tbl/.gff}
#			cat ${tblfile/.tbl/.gff}
#			echo start
			grep "5'UTR" ${tblfile/.tbl/.gff} | sed 's/^\([^\t]*\)_.[^\t]*/\1/' | tee ${tblfile/.tbl/_mod.gff}
			for gene in core E1 E2 p7 NS2 NS3 NS4A NS4B NS5A NS5B ; do
				target=($(grep $gene ${tblfile/.tbl/.gff}))
				printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" ${target[0]%_pilon} ${target[1]} "gene" ${target[3]} ${target[4]} ${target[5]} ${target[6]} "0" "ID:${gene};gene=${gene};product=${gene} protein;" | tee -a ${tblfile/.tbl/_mod.gff}
			done
			grep "3'UTR" ${tblfile/.tbl/.gff} | sed 's/^\([^\t]*\)_.[^\t]*/\1/' | tee -a ${tblfile/.tbl/_mod.gff}
		fi
	done
done

# takes a folder as input and output a bed file with the QCMD mutations


ns3mut=(36 54 55 56 80 122 132 155 156 168 170 175)
ns5amut=(28 30 31 58 93)
ns5bmut=(159 282 316 321 414 448 495 553 554 556)

ns3startpos=$(grep 'ID:NS3' ${resultsdir}/vadr/${id}*_mod.gff | cut -f4)
ns5astartpos=$(grep 'ID:NS5A' ${resultsdir}/vadr/${id}*_mod.gff | cut -f4)
ns5bstartpos=$(grep 'ID:NS5B' ${resultsdir}/vadr/${id}*_mod.gff | cut -f4)

[[ -f ${resultsdir}/vadr/${id}.vadr.bed ]] && rm ${resultsdir}/vadr/${id}.vadr.bed

for i in ${ns3mut[@]} ; do
        pos=$(( ${i}*3+${ns3startpos}-3 ))
        printf "%s\t%s\t%s\tNS3-%s\n" $id $(($pos-1)) $(($pos+2)) $i | tee -a ${resultsdir}/vadr/${id}.vadr.bed
done

for i in ${ns5amut[@]} ; do
        pos=$(( ${i}*3+${ns5astartpos}-3 ))
        printf "%s\t%s\t%s\tNS5A-%s\n" $id $(($pos-1)) $(($pos+2)) $i | tee -a ${resultsdir}/vadr/${id}.vadr.bed
done

for i in ${ns5bmut[@]} ; do
        pos=$(( ${i}*3+${ns5bstartpos}-3 ))
        printf "%s\t%s\t%s\tNS5B-%s\n" $id $(($pos-1)) $(($pos+2))  $i | tee -a ${resultsdir}/vadr/${id}.vadr.bed
done

