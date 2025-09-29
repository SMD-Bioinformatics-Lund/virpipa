#!/bin/bash

fulldir=$(readlink -f $1)
runname=$(basename $fulldir)
subsample=1000000
outdir=/fs1/jonas/hcv/results/
logdir=/fs1/jonas/hcv/logs/

if [[ ! -d $fulldir ]] ; then
	echo $fulldir is not a directory
	exit
fi

cd $logdir

for sample in ${fulldir}/*R1*gz ; do
	jobname=$(basename $sample)
	jobname=${jobname%%_*}
	sbatch -J HCV-${jobname} /fs1/jonas/hcv/scripts/hcvpipe.sh -f -s $subsample -o ${outdir}/${runname} $sample
done

cd - > /dev/null
