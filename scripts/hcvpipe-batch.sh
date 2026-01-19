#!/bin/bash

subsample=1000000
outdir=/fs1/jonas/hcv/results/
logdir=/fs1/jonas/hcv/logs/
partition=low
dry=''

showhelp(){
	echo hcvpipe_batch.sh lskdjfsldkfjlskdfj
}

readopts=$(getopt -o hno:l:p:i:c: --long help,dryrun,outdir:,logdir:,partition:,input:,csv: -n 'error' -- "$@")
eval set -- "$readopts"
dry=''

while true ; do
	case "$1" in
		-h|--help)
			showhelp
			exit 1 ;;
		-n|--dryrun)
			dry='echo'
			shift ;;
		-o|--outdir)
			outdir="$2"
			shift 2 ;;
		-l|--logdir)
			logdir="$2"
			shift 2 ;;
		-p|--partition)
			partition="$2"
			shift 2 ;;
		-i|--inputdir)
			inputdir="$2"
			shift 2 ;;
		-c|--csv)
			csv="$2"
			shift 2 ;;
		--) shift ; break ;;
		*) echo "unknown parameter" ; exit 1 ;;
	esac
done


sbatch_sample(){
	local jobname="$1"
	local sample="$2"
	local lid="$3"
	if [[ ! -z "$lid" ]] ; then
		$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -l $lid -o ${outdir}/${runname} $sample
	else
		$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -o ${outdir}/${runname} $sample
	fi	
}

runraw(){
	# Uses a directory with fastq files
	fulldir=$(readlink -f $inputdir)
	runname=$(basename $fulldir)
	for sample in ${inputdir}/*R1*gz ; do
		jobname=$(basename $sample)
		jobname=${jobname%%_*}
		if [[ -f "$inputdir"/../"$runname".tsv ]] ; then
			lid=$(grep $jobname "$inputdir"/../"$runname".tsv | cut -f2)
			$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -l $lid -o ${outdir}/${runname} $sample
		else
			$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -o ${outdir}/${runname} $sample
		fi
		if [[ -f "$fulldir"/../"$runname".tsv ]] ; then
			lid=$(grep $jobname "$fulldir"/../"$runname".tsv | cut -f2)
		fi
		sbatch_sample "$jobname" "$sample" "$lid"
	done
}

runcsv(){
	# Uses a csv file, parses the header and then created the sbatch commands
	declare -A csvdata
	{
		IFS=',' read -r -a keys
		IFS=',' read -r -a values
	} < "$csv"
	for i in "${!keys[@]}"; do
		# Remove potential carriage returns (\r) from Windows-style CSVs
		key=$(echo "${keys[$i]}" | tr -d '\r')
		val=$(echo "${values[$i]}" | tr -d '\r')
		csvdata["$key"]="$val"
	done
	echo CSV contents:
	for k in "${!csvdata[@]}"; do
		echo "$k: ${csvdata[$k]}"
	done
	echo
	# send the paramteres
		sbatch_sample ${csvdata["clarity_sample_id"]} ${csvdata["read1"]} ${csvdata["sample_name"]}
# 	if [[ -z ${csvdata["sample_name"]} ]] ; then
# 		echo run without lid
# 		$dry sbatch -J HCV-${csvdata["clarity_sample_id"]} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -o ${outdir}/${runname} ${csvdata["read1"]}
# 	else
# 		echo run with lid
# 		$dry sbatch -J HCV-${csvdata["clarity_sample_id"]} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -l ${csvdata["sample_name"]} -o ${outdir}/${runname} ${csvdata["read1"]}
# 		sbatch_sample ${csvdata["clarity_sample_id"]} ${csvdata["read1"]} ${csvdata["sample_name"]}
# 	fi
}

### main program

if [[ ! -d $outdir ]] ; then
	echo $outdir is not a directory
	exit
fi

cd $logdir

if [[ ! -z "$csv" ]] ; then
	runcsv
else
	runraw
fi

cd - > /dev/null
