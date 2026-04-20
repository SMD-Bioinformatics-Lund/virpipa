#!/bin/bash

shopt -s nullglob 

run_base_dir='/fs2/seqdata/NovaSeq/'
if [[ "$2" == 'header' ]] ; then
	printf "clarityID\tKM-LID\tgenotypCMD\taf005\taf01\taf015\taf02\taf03\taf04\tcdx1\tcdx10\tcdx100\tcdx1000\thostileremoved\ttotalreads\tdepth\tcoverage\trun\tct\tlibconc\tlibfrag\n"
fi

inputdir="$1"

for rundir in "$inputdir" ; do
	[[ ! -d "$rundir" ]] && continue
	for runpath in $rundir/* ; do
	[[ ! -d "$runpath" ]] && continue
		clarityid=$(basename $runpath)
		rundir=$(dirname $runpath)
		rundir=$(basename $rundir)

		blast=$(sed -n '2p' ${runpath}/${clarityid}.fasta.blast | cut -f1,2 | sed 's/.*\t//' | sed 's/_.*//')
		af005=$(sed -n 's/^AF\t0\t0.000000\t\([0-9]*\)\t.*/\1/p' ${runpath}/${clarityid}-pilon-m0.05.vcf.gz.stats)
		af01=$(sed -n 's/^AF\t0\t0.000000\t\([0-9]*\)\t.*/\1/p' ${runpath}/${clarityid}-pilon-m0.1.vcf.gz.stats)
		af015=$(sed -n 's/^AF\t0\t0.000000\t\([0-9]*\)\t.*/\1/p' ${runpath}/${clarityid}-pilon-m0.15.vcf.gz.stats)
		af02=$(sed -n 's/^AF\t0\t0.000000\t\([0-9]*\)\t.*/\1/p' ${runpath}/${clarityid}-pilon-m0.2.vcf.gz.stats)
		af03=$(sed -n 's/^AF\t0\t0.000000\t\([0-9]*\)\t.*/\1/p' ${runpath}/${clarityid}-pilon-m0.3.vcf.gz.stats)
		af04=$(sed -n 's/^AF\t0\t0.000000\t\([0-9]*\)\t.*/\1/p' ${runpath}/${clarityid}-pilon-m0.4.vcf.gz.stats)
		hostile=$(sed -n 's/.*"reads_removed_proportion": \(.*\),/\1/p' ${runpath}/hostile.json)
		reads=$(sed -n 's/ *"reads_in": \([0-9]*\).*/\1/p' ${runpath}/hostile.json)
		depth=$(sed -n 's/meandepth\t\(.*\)/\1/p' ${runpath}/${clarityid}-0.15-iupac.report.tsv)
		coverage=$(sed -n 's/coverage\t\(.*\)/\1/p' ${runpath}/${clarityid}-0.15-iupac.report.tsv)
		covdepth=$(sed -n '2p' ${runpath}/${clarityid}-coverage.tsv | cut -f2,3,4,5)
		#[[ $(compgen -G "$runpath/*.lid") != '' ]] && lid=$(echo $"$runpath"/*.lid | sed 's:.*/\(.*\)\.lid:\1:')
		lidfiles=("$runpath"/*.lid)
		if (( ${#lidfiles[@]} > 0 )) ; then
			lid="${lidfiles[0]}"
			lid="${lid##*/}"
			lid="${lid%.lid}"
		else
			lid=''
		fi
		read -r ct libconc libfrag < <(jq -r --arg id "$clarityid" 'to_entries[] | select(.value.clarity_sample_id == $id) | [.value.CT,.value."Library concentration (ng/ul)",.value."Library fragment length (bp)"] | @tsv' "$run_base_dir"/"$inputdir"/clarity_sample_info.json)
		printf "${clarityid}\t${lid}\t${blast}\t${af005}\t${af01}\t${af015}\t${af02}\t${af03}\t${af04}\t${covdepth}\t${hostile}\t${reads}\t${depth}\t${coverage}\t${rundir}\t${ct}\t${libconc}\t${libfrag}\n"
	done
done
