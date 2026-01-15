#!/bin/bash

fulldir=$(readlink -f $1)
runname=$(basename $fulldir)
subsample=1000000
outdir=/fs1/jonas/hcv/results/
logdir=/fs1/jonas/hcv/logs/
partition=low
dry=''
[[ "$2" == 'dry' ]] && dry='echo'

if [[ ! -d $fulldir ]] ; then
	echo $fulldir is not a directory
	exit
fi

cd $logdir

for sample in ${fulldir}/*R1*gz ; do
	jobname=$(basename $sample)
	jobname=${jobname%%_*}
	if [[ -f "$fulldir"/../"$runname".tsv ]] ; then
		lid=$(grep $jobname "$fulldir"/../"$runname".tsv | cut -f2)
		$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -l $lid -o ${outdir}/${runname} $sample
	else
		$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -o ${outdir}/${runname} $sample
	fi
done

cd - > /dev/null
