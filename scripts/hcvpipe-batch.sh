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
	$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -f -s $subsample -o ${outdir}/${runname} $sample
done

cd - > /dev/null
